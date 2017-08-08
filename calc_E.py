# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25th 2015

@author: flamholz
"""

from preprocessing import Preprocessing, Reaction, FARADAY

import argparse
import csv
from numpy import arange, sqrt
import logging

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
                                     'Calculate reduction potentials for '
                                     'a number of reactions.')
    parser.add_argument('infile', type=argparse.FileType(),
                        help='path to input file containing a '
                        'list of reactions')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help='path to output file')
    parser.add_argument('--ionic_strength', default=0.2, type=int,
                        help='ionic strength in molar units.')
    parser.add_argument('--pH_min', default=5, type=int,
                        help='lowest pH to produce E0 for.')
    parser.add_argument('--pH_max', default=9, type=int,
                        help='highest pH to produce E0 for.')
    parser.add_argument('--pH_step', default=0.05, type=float,
                        help='pH increment.')
    logging.getLogger().setLevel(logging.WARNING)

    args = parser.parse_args()

    ionic_strength = args.ionic_strength

    pHs = arange(args.pH_min, args.pH_max + args.pH_step, args.pH_step)

    prep = Preprocessing()

    reactions_and_energies = []
    reader = csv.reader(args.infile)
    for row in reader:
        formula = row[0].strip()
        reaction = Reaction.parse_formula(formula)
        n_e = reaction.check_half_reaction_balancing()
        if n_e == 0:
            logging.warning('one of the input reactions is not a half-reaction'
                            ': ' + formula)
            continue
        E0s = []
        for pH in pHs:
            dG0_prime, U = prep.dG0_prime(reaction, pH=pH,
                                          ionic_strength=ionic_strength)
            uncertainty = sqrt(U[0, 0])
            E0_prime = 1000.0 * -dG0_prime / (n_e*FARADAY)  # mV
            E0s.append(E0_prime)

        reactions_and_energies.append((row, E0s))

    header = ['reaction']
    pH_header = ['pH %.1f (mV)' % pH for pH in pHs]
    header += pH_header
    writer = csv.writer(args.outfile)
    writer.writerow(header)
    for rxn_data, pH_E0 in reactions_and_energies:
        energies_fmted = ['%.2f' % e for e in pH_E0]
        writer.writerow(rxn_data + energies_fmted)
