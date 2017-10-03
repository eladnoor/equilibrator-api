#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 13:28:00 2017

@author: noore
"""
import numpy as np
import matplotlib.pyplot as plt

class ReactionMDFData(object):

    def __init__(self, reaction, flux, dGr, shadow_price):
        """
        Args:
            reaction: kegg reaction object.
                should be set to contain user-defined dG0
            flux: amount of relative flux in pathway.
            dGr: dG in MDF conditions.
            shadow_price: shadow price associated with this rxn.
        """
        self.reaction = reaction
        self.flux = flux
        self.dGr = dGr
        self.shadow_price = shadow_price

class CompoundMDFData(object):
    def __init__(self, compound, concentration_bounds,
                 concentration, shadow_price):
        self.compound = compound
        self.concentration = concentration
        self.shadow_price = shadow_price
        self.lb, self.ub = concentration_bounds

    @property
    def bounds_equal(self):
        return self.lb == self.ub


class PathwayMDFData(object):

    def __init__(self, parsed_pathway, mdf_result):
        self.parsed_pathway = parsed_pathway
        self.mdf_result = mdf_result
        self.model = mdf_result.model

        rxns = parsed_pathway.reactions
        fluxes = parsed_pathway.fluxes
        dGs = self.mdf_result.dG_r_prime_adj.flatten().tolist()[0]
        prices = self.mdf_result.reaction_prices.flatten().tolist()[0]
        self.reaction_data = [
            ReactionMDFData(*t) for t in zip(rxns, fluxes, dGs, prices)]

        compound_ids = parsed_pathway.compound_kegg_ids
        cbounds = [self.model.concentration_bounds.GetBoundTuple(cid)
                   for cid in compound_ids]
        concs = self.mdf_result.concentrations.flatten().tolist()[0]
        prices = self.mdf_result.compound_prices.flatten().tolist()[0]
        self.compound_data = [CompoundMDFData(*t)
                              for t in zip(compound_ids, cbounds, concs, prices)]

    @property
    def mdf(self):
        return self.mdf_result.mdf

    @property
    def min_total_dG(self):
        return self.mdf_result.min_total_dG

    @property
    def max_total_dG(self):
        return self.mdf_result.max_total_dG

    @property
    def max_total_driving_force(self):
        return -self.min_total_dG

    @property
    def min_total_driving_force(self):
        return -self.max_total_dG

    @property
    def conc_plot(self):
        ys = np.arange(0, len(self.compound_data))
        concs = np.array([c.concentration for c in self.compound_data])
        cnames = [str(c.compound) for c in self.compound_data]
        default_lb = self.model.concentration_bounds.default_lb
        default_ub = self.model.concentration_bounds.default_ub

        cids = [str(c.compound) for c in self.compound_data]
        lbs = [self.model.concentration_bounds.GetLowerBound(cid)
               for cid in cids]
        ubs = [self.model.concentration_bounds.GetUpperBound(cid)
               for cid in cids]
        lbs, ubs = np.array(lbs), np.array(ubs)
        bounds_equal = np.where(lbs == ubs)
        ys_equal = ys[bounds_equal]
        concs_equal = concs[bounds_equal]

        # Special color for metabolites with nonzero shadow prices.
        shadow_prices = np.array([c.shadow_price for c in self.compound_data])
        nz_shadow = np.where(shadow_prices != 0)
        ys_nz_shadow = ys[nz_shadow]
        concs_nz_shadow = concs[nz_shadow]

        conc_figure, ax = plt.subplots(1, 1, figsize=(6, 6))
        #ax = plt.axes([0.2, 0.1, 0.7, 0.8])
        ax.axvspan(1e-8, default_lb, color='y', alpha=0.5)
        ax.axvspan(default_ub, 1e3, color='y', alpha=0.5)
        ax.scatter(concs, ys,
                   label='Variable Concentrations')
        ax.scatter(concs_equal, ys_equal, color='y',
                   label='Fixed Concentrations')
        ax.scatter(concs_nz_shadow, ys_nz_shadow,
                   color='r', label='Variable Concentrations')

        ax.set_yticks(ys)
        ax.set_yticklabels(cnames, fontsize=9)
        ax.set_xlabel('Concentration (M)')
        ax.set_xscale('log')

        ax.set_xlim(1e-7, 1.5)
        ax.set_ylim(-1.5, len(self.compound_data) + 0.5)
        
        conc_figure.tight_layout()
        return conc_figure
    
    @property
    def mdf_plot(self):
        dgs = [0] + [r.dGr for r in self.reaction_data]
        cumulative_dgs = np.cumsum(dgs)

        xs = np.arange(0, len(cumulative_dgs))

        mdf_fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        ax.plot(xs, cumulative_dgs, label='MDF-optimized concentrations')
        ax.set_xticks(xs)

        # TODO: Consider using reaction IDs from the file as xticks?
        ax.set_xlabel('After Reaction Step')
        ax.set_ylabel("Cumulative $\Delta_r G'$ (kJ/mol)")
        ax.legend(loc='best')
        ax.set_title('MDF = %.1f' % self.mdf)

        mdf_fig.tight_layout()
        return mdf_fig
