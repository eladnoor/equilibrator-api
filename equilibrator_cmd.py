#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 15:36:47 2017

@author: noore

A stand-alone version of Component Contribution that can calculate the
Delta-Gr'0 of any reaction (with KEGG notation, i.e. whose reactants are
already cached in our database), at a given pH and I.
"""

import argparse
import logging
import sys
from numpy import sqrt
from preprocessing import Preprocessing, Reaction, FARADAY


def MakeParser():
    parser = argparse.ArgumentParser(
        description=('Estimate the Gibbs energy of a reaction. For example,'
                     'the following calculates dGr0 for ATP hydrolysis '
                     'at pH 6: calc_dGr0.py --ph 6 "C00002 + C00001 = '
                     'C00008 + C00009"'))
    parser.add_argument('--ph', type=float, help='pH level', default=7.0)
    parser.add_argument('--i', type=float,
                        help='ionic strength in M',
                        default=0.1)
    parser.add_argument('reaction', type=str, help='reaction in KEGG notation')
    return parser


###############################################################################
parser = MakeParser()
args = parser.parse_args()

logging.getLogger().setLevel(logging.WARNING)

print 'pH = %.1f' % args.ph
print 'I = %.1f M' % args.i
print 'Reaction: ' + args.reaction

# parse the reaction
reaction = Reaction.parse_formula(args.reaction)

prep = Preprocessing()

n_e = reaction.check_half_reaction_balancing()
dG0_prime, U = prep.dG0_prime(reaction, pH=args.ph,
                              ionic_strength=args.i)
dG0_prime = dG0_prime[0, 0]
uncertainty = sqrt(U[0, 0])

if n_e != 0:  # treat as a half-reaction
    E0_prime_mV = 1000.0 * -dG0_prime / (n_e*FARADAY)
    E0_uncertainty = 1000.0 * uncertainty / (n_e*FARADAY)
    print u'E\'\u00B0 = %.1f \u00B1 %.1f mV' % (E0_prime_mV, E0_uncertainty)
else:
    print u'\u0394G\'\u00B0 = %.1f \u00B1 %.1f kJ/mol' % (dG0_prime, uncertainty)

print r'* the range represents the 95% confidence interval due to ' + \
       'Component Contribution estimation uncertainty'