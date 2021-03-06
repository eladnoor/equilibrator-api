#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 13:26:31 2017

@author: noore
"""
import csv
import logging
import numpy as np
import io
from scipy import linalg
from util.SBtab import SBtabTools

from equilibrator_api.reaction import Reaction
from equilibrator_api.bounds import DEFAULT_BOUNDS, Bounds
from equilibrator_api.max_min_driving_force import PathwayMDFData
from equilibrator_api.thermo_models import PathwayThermoModel
from equilibrator_api.component_contribution import ComponentContribution
from equilibrator_api.settings import RT, DEFAULT_PH, DEFAULT_IONIC_STRENGTH

class PathwayParseError(Exception):
    pass


class InvalidReactionFormula(PathwayParseError):
    pass


class UnbalancedReaction(PathwayParseError):
    pass


class ViolatesFirstLaw(PathwayParseError):
    pass


class Pathway(object):
    """A pathway parsed from user input.

    Designed for checking input prior to converting to a stoichiometric model.
    """

    def __init__(self, reactions, fluxes, dG0_r_primes, name_to_cid=None,
                 bounds=None, pH=DEFAULT_PH,
                 ionic_strength=DEFAULT_IONIC_STRENGTH):
        """Initialize.

        Args:
            reactions: a list of gibbs.reaction.Reaction objects.
            fluxes: numpy.array of relative fluxes in same order as reactions.
            dG0_r_primes: reaction energies.
            bounds: bounds on metabolite concentrations.
                Uses default bounds if None provided.
            pH: (optional) specify the pH at which the dG values are calculated
            ionic_strength: (optional) specify the I at which the dG values are calculated
        """
        assert len(reactions) == len(fluxes)
        assert len(reactions) == len(dG0_r_primes)

        self.reactions = reactions
        self.fluxes = np.array(fluxes)
        self.dG0_r_prime = np.array(dG0_r_primes)
        dGm_corr = np.array([r.dGm_correction() for r in self.reactions])
        self.dGm_r_prime = self.dG0_r_prime + dGm_corr
        self.bounds = bounds or DEFAULT_BOUNDS

        self.S, self.compound_kegg_ids = self._build_stoichiometric_matrix()

        if name_to_cid is not None:
            # invert the dictionary and convert the cids back to names
            cid_to_name = dict(zip(name_to_cid.values(), name_to_cid.keys()))
            self.compound_names = list(map(cid_to_name.get,
                                           self.compound_kegg_ids))
        else:
            self.compound_names = list(self.compound_kegg_ids)
            
        nr, nc = self.S.shape

        # dGr should be orthogonal to nullspace of S
        # If not, dGr is not contained in image(S) and then there
        # is no consistent set of dGfs that generates dGr and the
        # first law of thermo is violated by the model.
        S_T = np.matrix(self.S).T
        Spinv = linalg.pinv(S_T)
        null_proj = np.matrix(np.eye(S_T.shape[0])) - S_T*Spinv
        projected = null_proj * np.matrix(self.dG0_r_prime).T
        if not np.all(projected < 1e-8):
            raise ViolatesFirstLaw(
                'Supplied reaction dG values are inconsistent '
                'with the stoichiometric matrix.')

    @classmethod
    def from_csv_file(cls, f, bounds=None, pH=None, ionic_strength=None):
        """Returns a pathway parsed from an input file.

        Caller responsible for closing f.

        Args:
            f: file-like object containing CSV data describing the pathway.
        """
        reactions = []
        fluxes = []

        for row in csv.DictReader(f):
            rxn_formula = row.get('ReactionFormula')

            flux = float(row.get('Flux', 0.0))
            logging.debug('formula = %f x (%s)', flux, rxn_formula)

            rxn = Reaction.parse_formula(rxn_formula)
            rxn.check_full_reaction_balancing()

            reactions.append(rxn)
            fluxes.append(flux)

        equilibrator = ComponentContribution(pH=pH, ionic_strength=ionic_strength)
        dG0_r_primes, dG0_uncertainties = zip(*map(equilibrator.dG0_prime, reactions))
        dG0_r_primes = list(dG0_r_primes)
        return Pathway(reactions, fluxes, dG0_r_primes, bounds=bounds)

    def _get_compounds(self):
        """Returns a dictionary of compounds by KEGG ID."""
        compounds = {}
        for r in self.reactions:
            for cw_coeff in r.kegg_ids():
                c = cw_coeff.compound
                compounds[c.kegg_id] = c
        return compounds

    def _build_stoichiometric_matrix(self):
        """Builds a stoichiometric matrix.

        Returns:
            Two tuple (S, compounds) where compounds is the KEGG IDs of the compounds
            in the order defining the column order of the stoichiometric matrix S.
        """
        compounds = set()
        for r in self.reactions:
            compounds.update(r.kegg_ids())
        compounds = sorted(compounds)
        
        smat = np.matrix(np.zeros((len(compounds), len(self.reactions))))
        for j, r in enumerate(self.reactions):
            for i, c in enumerate(compounds):
                smat[i, j] = r.get_coeff(c)

        return smat, compounds

    def calc_mdf(self):
        dGs = np.matrix(self.dG0_r_prime).T
        model = PathwayThermoModel(self.S, self.fluxes, dGs,
                                   self.compound_kegg_ids,
                                   concentration_bounds=self.bounds)
        mdf = model.mdf_result
        return PathwayMDFData(self, mdf)

    def print_reactions(self):
        for f, r in zip(self.fluxes, self.reactions):
            print('%sx %s' % (f, r.write_formula()))

    @classmethod
    def from_sbtab(self, sbtab):
        """
            read the sbtab file (can be a filename or file handel)
            and use it to initialize the Pathway
        """
        if type(sbtab) == str:
            with open(sbtab, 'r') as sbtabfile:
                sbtabs = SBtabTools.openMultipleSBtabFromFile(sbtabfile)
        elif type(sbtab) == io.TextIOWrapper:
            sbtabs = SBtabTools.openMultipleSBtabFromFile(sbtab)
        tdict = dict([(t.getTableInformation()[1].upper(), t) for t in sbtabs])
        expected_tnames = ['REACTION', 'RELATIVEFLUX', 'CONCENTRATIONCONSTRAINT',
                           'REACTIONCONSTANT']
        assert set(expected_tnames).issubset(tdict.keys())
    
        sbtabs = [tdict[n] for n in expected_tnames]
        return Pathway.from_separate_sbtabs(*sbtabs)

    @classmethod
    def from_separate_sbtabs(self, reaction_sbtab, flux_sbtab,
                             bounds_sbtab, keqs_sbtab):
        """Returns an initialized Pathway."""
        bounds = Bounds.from_sbtab(bounds_sbtab)

        reaction_df = reaction_sbtab.toDataFrame()
        flux_df = flux_sbtab.toDataFrame()
        bounds_df = bounds_sbtab.toDataFrame()
        keqs_df = keqs_sbtab.toDataFrame()

        name_to_cid = dict(
            zip(bounds_df['Compound'],
                bounds_df['Compound:Identifiers:kegg.compound']))

        reactions = []
        for _, row in reaction_df.iterrows():
            rxn = Reaction.parse_formula(row['ReactionFormula'], name_to_cid,
                                         rid=row['ID'])
            rxn.check_full_reaction_balancing()
            reactions.append(rxn)

        reaction_ids = reaction_df['ID']
        fluxes = flux_df[flux_df['QuantityType'] == 'flux']
        reaction_fluxes = dict(zip(fluxes['Reaction'], fluxes['Value']))
        fluxes_ordered = [float(reaction_fluxes[rid]) for rid in reaction_ids]

        # grab rows containing keqs.
        keqs = keqs_df[keqs_df['QuantityType'] == 'equilibrium constant']
        reaction_keqs = dict(zip(keqs['Reaction'], keqs['Value']))
        dgs = [-RT * np.log(float(reaction_keqs[rid]))
               for rid in reaction_ids]

        # Manually set the delta G values on the reaction objects
        for dg, rxn in zip(dgs, reactions):
            rxn._dg0_prime = dg

        pH = keqs_sbtab.getCustomTableInformation('pH')
        ionic_strength = keqs_sbtab.getCustomTableInformation('IonicStrength')
        pp = Pathway(reactions, fluxes_ordered, dgs,
                     name_to_cid=name_to_cid,
                     bounds=bounds, pH=pH, ionic_strength=ionic_strength)
        return pp
