#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 12:59:43 2017

@author: noore
"""
import json
import os
import re
import logging
from numpy import matrix, zeros, log, nan, inf

from equilibrator_api import settings
from equilibrator_api.compound import Compound

COMPOUND_JSON_FNAME = os.path.join(settings.DATA_DIR, 'cc_compounds.json')

class Reaction(object):
    # load formation energies from the JSON file
    COMPOUND_DICT = {}
    for cd in json.load(open(COMPOUND_JSON_FNAME, 'r')):
        kegg_id = cd.get('CID', 'unknown')
        COMPOUND_DICT[kegg_id] = cd

    REACTION_COUNTER = 0

    def __init__(self, kegg_id_to_coeff, rid=None):
        self.kegg_id_to_coeff = kegg_id_to_coeff

        # Create the relevant "Compound" objects and store in a dictionary
        self.kegg_id_to_compound = {}
        for kegg_id in self.kegg_id_to_coeff.keys():
            compound = Compound(Reaction.COMPOUND_DICT[kegg_id])
            self.kegg_id_to_compound[kegg_id] = compound
        
        if rid is not None:
            self.reaction_id = rid
        else:
            self.reaction_id = 'R%05d' % Reaction.REACTION_COUNTER
            Reaction.REACTION_COUNTER += 1
        
    def kegg_ids(self):
        return self.kegg_id_to_coeff.keys()

    def is_empty(self):
        return len(self.kegg_id_to_coeff) == 0

    def get_coeff(self, kegg_id):
        return self.kegg_id_to_coeff.get(kegg_id, 0)

    def get_compound(self, kegg_id):
        return self.kegg_id_to_compound.get(kegg_id, None)

    def dG0_prime(self, pH=settings.DEFAULT_PH, pMg=settings.DEFAULT_PMG,
                  ionic_strength=settings.DEFAULT_IONIC_STRENGTH):
        dG0_r_prime = 0
        for kegg_id in self.kegg_ids():
            coeff = self.get_coeff(kegg_id)
            compound = self.get_compound(kegg_id)
            dG0_r_prime += coeff * compound.dG0_prime(pH, pMg, ionic_strength)
        return dG0_r_prime

    def _GetSumCoeff(self):
        """
            Calculate the sum of all coefficients (excluding water).
            This is useful for shifting the dG'0 to another set of standard
            concentrations (e.g. 1 mM)
        """
        ids = set(self.kegg_ids()) - set(['C00001'])
        sum_coeff = sum(map(self.get_coeff, ids))
        return sum_coeff

    def _GetAbsSumCoeff(self):
        """
            Calculate the sum of all coefficients (excluding water) in
            absolute value.
            This is useful for calculating the reversibility index.
        """
        ids = set(self.kegg_ids()) - set(['C00001'])
        abs_sum_coeff = sum(map(abs, (map(self.get_coeff, ids))))
        return abs_sum_coeff     

    def dGm_correction(self):
        """
            Calculate the dG' in typical physiological concentrations (1 mM)
        """
        return settings.RT * self._GetSumCoeff() * log(1e-3)

    def dGm_prime(self):
        """
            Calculate the dG' in typical physiological concentrations (1 mM)
        """
        return self.dG0_prime() + self.dGm_correction()

    def reversibility_index(self, pH=settings.DEFAULT_PH, pMg=settings.DEFAULT_PMG,
                            ionic_strength=settings.DEFAULT_IONIC_STRENGTH):
        """
            Calculates the reversiblity index according to Noor et al. 2012:
            https://doi.org/10.1093/bioinformatics/bts317

            Returns:
                ln_RI - the natural log of the RI
        """
        dG0_prime = self.dG0_prime(pH, pMg, ionic_strength)
        return self.calculate_reversibility_index_from_dG0_prime(dG0_prime)

    def calculate_reversibility_index_from_dG0_prime(self, dG0_prime):
        """
            Calculates the reversiblity index according to Noor et al. 2012:
            https://doi.org/10.1093/bioinformatics/bts317

            Returns:
                ln_RI - the natural log of the RI
        """
        sum_coeff = self._GetSumCoeff()
        abs_sum_coeff = self._GetAbsSumCoeff()
        if abs_sum_coeff == 0:
            return inf
        dGm_prime = dG0_prime + settings.RT * sum_coeff * log(1e-3)
        ln_RI = (2.0 / abs_sum_coeff) * dGm_prime / settings.RT
        return ln_RI

    @staticmethod
    def parse_formula_side(s):
        """
            Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
            Ignores stoichiometry.

            Returns:
                The set of CIDs.
        """
        if s.strip() == "null":
            return {}

        compound_bag = {}
        for member in re.split('\s+\+\s+', s):
            tokens = member.split(None, 1)  # check for stoichiometric coeff
            if len(tokens) == 0:
                continue
            if len(tokens) == 1:
                amount = 1
                key = member
            else:
                try:
                    amount = float(tokens[0])
                except ValueError:
                    raise ValueError('could not parse the reaction side: %s'
                                     % s)
                key = tokens[1]
            compound_bag[key] = compound_bag.get(key, 0) + amount

        return compound_bag

    @staticmethod
    def parse_formula(formula, name_to_cid=None, rid=None):
        """
            Parses a two-sided formula such as: 2 C00001 = C00002 + C00003
            
            Args:
                formula     - a string representation of the chemical formula
                name_to_cid - (optional) a dictionary mapping names to KEGG IDs

            Return:
                The set of substrates, products and the reaction direction
        """
        tokens = []
        for arrow in settings.POSSIBLE_REACTION_ARROWS:
            if formula.find(arrow) != -1:
                tokens = formula.split(arrow, 2)
                break

        if len(tokens) < 2:
            raise ValueError('Reaction does not contain an allowed arrow sign:'
                             ' %s' % (arrow, formula))

        left = tokens[0].strip()
        right = tokens[1].strip()

        sparse_reaction = {}
        for cid, count in Reaction.parse_formula_side(left).items():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count

        for cid, count in Reaction.parse_formula_side(right).items():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count

        # remove compounds that are balanced out in the reaction,
        # i.e. their coefficient is 0
        sparse_reaction = dict(filter(lambda x: x[1] != 0,
                                      sparse_reaction.items()))
        if name_to_cid is not None:
            # replace the names of the metabolites with their KEGG IDs
            # using this dictionary
            sparse_reaction = \
                dict(zip(map(name_to_cid.get, sparse_reaction.keys()),
                     sparse_reaction.values()))
            
        if 'C00080' in sparse_reaction:
            sparse_reaction.pop('C00080')
        
        return Reaction(sparse_reaction, rid=rid)

    @staticmethod
    def write_compound_and_coeff(compound_id, coeff):
        if coeff == 1:
            return compound_id
        else:
            return "%g %s" % (coeff, compound_id)

    def write_formula(self):
        """String representation."""
        left = []
        right = []
        for kegg_id in self.kegg_ids():
            coeff = self.get_coeff(kegg_id)
            if coeff < 0:
                left.append(Reaction.write_compound_and_coeff(kegg_id, -coeff))
            elif coeff > 0:
                right.append(Reaction.write_compound_and_coeff(kegg_id, coeff))
        return "%s %s %s" % (' + '.join(left), '=', ' + '.join(right))

    def _get_element_matrix(self):
        # gather the "atom bags" of all compounds in a list 'atom_bag_list'
        elements = set()
        atom_bag_list = []
        for kegg_id in self.kegg_ids():
            comp = self.get_compound(kegg_id)
            atom_bag = comp.get_atom_bag()
            if atom_bag is not None:
                elements = elements.union(atom_bag.keys())
            atom_bag_list.append(atom_bag)
        elements = sorted(elements)

        # create the elemental matrix, where each row is a compound and each
        # column is an element (or e-)
        Ematrix = matrix(zeros((len(atom_bag_list), len(elements))))
        for i, atom_bag in enumerate(atom_bag_list):
            if atom_bag is None:
                Ematrix[i, :] = nan
            else:
                for j, elem in enumerate(elements):
                    Ematrix[i, j] = atom_bag.get(elem, 0)
        return elements, Ematrix

    def _get_reaction_atom_balance(self):
        cids = list(self.kegg_ids())
        coeffs = matrix(list(map(self.get_coeff, cids)))

        elements, Ematrix = self._get_element_matrix()
        conserved = coeffs * Ematrix

        atom_balance_dict = dict([(e, c) for (e, c) in
                                  zip(elements, conserved.flat) if (c != 0)])

        return atom_balance_dict

    def check_half_reaction_balancing(self):
        """
            Returns:
                The number of electrons that are 'missing' in the half-reaction
                or None if the reaction is not atomwise-balanced.
        """
        atom_balance_dict = self._get_reaction_atom_balance()
        n_e = atom_balance_dict.pop('e-', 0)
        if not self._check_balancing(atom_balance_dict):
            return None
        else:
            return n_e

    def check_full_reaction_balancing(self):
        """
            Returns:
                True iff the reaction is balanced for all elements
                (excluding H)
        """
        atom_balance_dict = self._get_reaction_atom_balance()
        return self._check_balancing(atom_balance_dict)

    def _check_balancing(self, atom_balance_dict):
        """
            Use for checking if all elements are conserved.

            Returns:
                An atom_bag of the differences between the sides of the
                reaction. E.g. if there is one extra C on the left-hand
                side, the result will be {'C': -1}.
        """
        if nan in atom_balance_dict.values():
            warning_str = 'cannot test reaction balancing because of ' + \
                          'unspecific compound formulas: %s' % \
                          self.write_formula()
            raise ValueError(warning_str)

        # if there are unbalanced elements, write a full report
        if len(atom_balance_dict) == 0:
            return True
        logging.warning('unbalanced reaction: %s' % self.write_formula())
        for elem, c in atom_balance_dict.items():
            if c != 0:
                logging.warning('there are %d more %s atoms on the '
                                'right-hand side' % (c, elem))
        return False

