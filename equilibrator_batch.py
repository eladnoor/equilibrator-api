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

        dG0_prime_list = []
        uncertainty_list = []
        for pH in pHs:
            dG0_prime, U = prep.dG0_prime(reaction, pH=pH,
                                          ionic_strength=ionic_strength)
            dG0_prime_list.append(dG0_prime[0, 0])
            uncertainty_list.append(sqrt(U[0, 0]))

        if n_e != 0:  # treat as a half-reaction
            E0_prime_mV = map(lambda g: '%.2f' % (1000.0 * -g / (n_e*FARADAY)),
                              dG0_prime_list)
            reactions_and_energies.append([formula, 'half-reaction',
                                           'E\'0', 'mV'] + E0_prime_mV)
        else:
            dG0_prime_kj_mol = map(lambda g: '%.2f' % g,
                                   dG0_prime_list)
            reactions_and_energies.append([formula, 'reaction', 'dG\'0',
                                           'kJ/mol'] + dG0_prime_kj_mol)

    
    # write all results to the output CSV file
    writer = csv.writer(args.outfile)
    header = ['formula', 'type', 'estimated_value', 'unit']
    header += ['pH %.1f' % pH for pH in pHs]
    writer.writerow(header)
    writer.writerows(reactions_and_energies)
    args.outfile.flush()
