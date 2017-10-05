#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 15:47:09 2017

@author: noore
"""

import logging
from numpy import matrix, load, zeros, sqrt
import os
from equilibrator_api import settings

PREPROCESS_FNAME = os.path.join(settings.DATA_DIR, 'cc_preprocess.npz')

class ComponentContribution(object):

    def __init__(self, pH=settings.DEFAULT_PH, pMg=settings.DEFAULT_PMG,
                 ionic_strength=settings.DEFAULT_IONIC_STRENGTH):
        # load pre-processing matrices (for the uncertainty estimation)
        relpath = os.path.dirname(os.path.realpath(__file__))
        cc_preprocess_fname = os.path.join(relpath, PREPROCESS_FNAME)
        cc_preprocess = load(cc_preprocess_fname)

        self.v_r = matrix(cc_preprocess['v_r'])
        self.v_g = matrix(cc_preprocess['v_g'])
        self.C1 = matrix(cc_preprocess['C1'])
        self.C2 = matrix(cc_preprocess['C2'])
        self.C3 = matrix(cc_preprocess['C3'])
        self.G1 = matrix(cc_preprocess['G1'])
        self.G2 = matrix(cc_preprocess['G2'])
        self.G3 = matrix(cc_preprocess['G3'])
        self.S = matrix(cc_preprocess['S'])
        self.kegg_ids = cc_preprocess['cids']
        self.Nc = self.C1.shape[0]
        self.Ng = self.C3.shape[0]
        assert self.C1.shape[0] == self.C1.shape[1]
        assert self.C1.shape[1] == self.C2.shape[0]
        assert self.C2.shape[1] == self.C3.shape[0]
        assert self.C3.shape[0] == self.C3.shape[1]
        assert self.C3.shape[0] == self.C3.shape[1]
        
        self.pH = pH
        self.ionic_strength = ionic_strength
        self.pMg = pMg

    def reaction_to_vectors(self, reaction):
        # x is the stoichiometric vector of the reaction, only for the
        # compounds that appeared in the original training set for CC
        x = matrix(zeros((self.Nc, 1)))

        # g is the group incidence vector of all the other compounds
        g = matrix(zeros((self.Ng, 1)))

        for kegg_id in reaction.kegg_ids():
            coeff = reaction.get_coeff(kegg_id)
            compound = reaction.get_compound(kegg_id)
            x += coeff * compound.get_stoich_vector(self.Nc)
            g += coeff * compound.get_group_incidence_vector(self.Ng)

        return x, g

    def reactions_to_matrices(self, reactions):
        """
            Arguments:
                reaction - a KeggReaction object

            Returns:
                X        - the stoichiometric matrix of the reactions (only
                           for compounds that appear in the original training
                           set of CC)
                G        - the group incidence matrix (of all other compounds)
        """
        X = matrix(zeros((self.Nc, len(reactions))))
        G = matrix(zeros((self.Ng, len(reactions))))

        for i, reaction in enumerate(reactions):
            x, g = self.reaction_to_vectors(reaction)
            X[:, i] = x
            G[:, i] = g

        return X, G

    def dG0_prime(self, reaction):
        """
            Calculate the dG'0 of a single reaction or a list of reactions
            
            Returns:
                dG0_r_prime     - estimated Gibbs free energy of reaction
                dG0_uncertainty - standard deviation of estimation, multiply
                                  by 1.96 to get a 95% confidence interval
                                  (which is the value shown on eQuilibrator)
        """

        dG0_r_prime = reaction.dG0_prime(pH=self.pH, pMg=self.pMg,
                                         ionic_strength=self.ionic_strength)

        X, G = self.reactions_to_matrices([reaction])
        U = X.T * self.C1 * X + \
            X.T * self.C2 * G + \
            G.T * self.C2.T * X + \
            G.T * self.C3 * G

        dG0_uncertainty = sqrt(U[0, 0])
        return dG0_r_prime, dG0_uncertainty

    def dG0_prime_multi(self, reactions):
        """
            Calculate the dG'0 of a single reaction or a list of reactions

            Returns:
                dG0_r_prime     - estimated Gibbs free energy of reaction
                U               - correlation matrix of the uncertainties
                                  the values on the diagonal are the variances,
                                  i.e. the square of the standard deviations
                                  that one would get using the dG0_prime()
                                  funtion.
        """

        if not isinstance(reactions, list):
            reactions = [reactions]

        dG0_r_prime = map(lambda r: r.dG0_prime(pH=self.pH, pMg=self.pMg,
                                                ionic_strength=self.ionic_strength),
                          reactions)
        dG0_r_prime = matrix(list(dG0_r_prime)).T

        X, G = self.reactions_to_matrices(reactions)
        U = X.T * self.C1 * X + \
            X.T * self.C2 * G + \
            G.T * self.C2.T * X + \
            G.T * self.C3 * G

        return dG0_r_prime, U

    def E0_prime(self, reaction):
        """
            Calculate the E'0 of a single half-reaction
        """
        n_e = reaction.check_half_reaction_balancing()
        if n_e is None:
            raise ValueError('reaction is not chemically balanced')
        if n_e == 0:
            raise ValueError('this is not a half-reaction, '
                             'electrons are balanced')

        dG0_prime, dG0_uncertainty = self.dG0_prime(reaction)

        E0_prime_mV = 1000.0 * -dG0_prime / (n_e*settings.FARADAY)
        E0_uncertainty = 1000.0 * dG0_uncertainty / (n_e*settings.FARADAY)

        return E0_prime_mV, E0_uncertainty

    @staticmethod
    def WriteCompoundAndCoeff(kegg_id, coeff):
        if coeff == 1:
            return kegg_id
        else:
            return "%g %s" % (coeff, kegg_id)

    @staticmethod
    def DictToReactionString(d):
        """String representation."""
        left = []
        right = []
        for kegg_id, coeff in sorted(d.items()):
            _s = ComponentContribution.WriteCompoundAndCoeff(kegg_id, -coeff)
            if coeff < 0:
                left.append(_s)
            elif coeff > 0:
                right.append(_s)
        return "%s %s %s" % (' + '.join(left), '=', ' + '.join(right))

    @staticmethod
    def Analyze(self, x, g):
        weights_rc = x.T * self.G1
        weights_gc = x.T * self.G2 + g.T * self.G3
        weights = weights_rc + weights_gc

        res = []
        for j in range(self.S.shape[1]):
            d = {self.kegg_ids[i]: self.S[i, j]
                 for i in range(self.Nc)
                 if self.S[i, j] != 0}
            r_string = self.DictToReactionString(d)
            res.append({'w': weights[0, j],
                        'w_rc': weights_rc[0, j].round(4),
                        'w_gc': weights_gc[0, j].round(4),
                        'reaction_string': r_string})
        res.sort(key=lambda d: abs(d['w']), reverse=True)
        return res

    def IsUsingGroupContributions(self, x, g):
        weights_gc = x.T * self.G2 + g.T * self.G3
        sum_w_gc = sum(abs(weights_gc).flat)
        logging.info('sum(w_gc) = %.2g' % sum_w_gc)
        return sum_w_gc > 1e-5

