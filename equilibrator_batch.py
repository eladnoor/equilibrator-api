# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25th 2015

@author: flamholz
"""

from equilibrator_api import EquilibratorAPI, Reaction

import argparse
import csv
from numpy import sqrt, nan
import logging
import sys

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Calculate potentials for a number of reactions.')
    parser.add_argument(
        'infile', type=argparse.FileType(),
        help='path to input file containing reactions')
    parser.add_argument(
        'outfile', type=argparse.FileType('w'),
        help='path to output file')
    parser.add_argument('--i', type=float,
                        help='ionic strength in M',
                        default=0.1)
    parser.add_argument('--ph', type=float, help='pH level', default=7.0)
    logging.getLogger().setLevel(logging.WARNING)

    args = parser.parse_args()

    sys.stderr.write('pH = %.1f\n' % args.ph)
    sys.stderr.write('I = %.1f M\n' % args.i)

    infile_lines = filter(None, map(str.strip, args.infile.readlines()))
    reactions = map(Reaction.parse_formula, infile_lines)

    equilibrator = EquilibratorAPI()

    dG0_prime, U = equilibrator.dG0_prime_multi(
            reactions, pH=args.ph, ionic_strength=args.i)

    writer = csv.writer(args.outfile)
    header = ['reaction', 'pH', 'ionic strength [M]', 'dG\'0 [kJ/mol]',
              'uncertainty [kJ/mol]', 'comment']
    writer.writerow(header)
    for i, r in enumerate(reactions):
        if r.check_full_reaction_balancing():
            row = [infile_lines[i], args.ph, args.i, dG0_prime[i, 0],
                   sqrt(U[i, i]), '']
        else:
            row = [infile_lines[i], args.ph, args.i, nan, nan,
                   'reaction is not chemically balanced']
        writer.writerow(row)

    args.outfile.flush()
