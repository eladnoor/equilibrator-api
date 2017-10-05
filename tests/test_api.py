#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 17:22:11 2017

@author: noore
"""

import unittest

from equilibrator_api.reaction import Reaction
from equilibrator_api.component_contribution import ComponentContribution
from equilibrator_api import settings
from equilibrator_api.pathway import Pathway
import os

class TestReactionParsing(unittest.TestCase):
    
    def test_atp_hydrolysis(self):
        formula = ' C00002 + C00001  <= C00008 +   C00009'
        kegg_ids = set(('C00002', 'C00001', 'C00008', 'C00009'))
        try:
            reaction = Reaction.parse_formula(formula)
        except ValueError as e:
            self.fail('unable to parse the formula\n' + str(e))
            
        self.assertSetEqual(set(reaction.kegg_ids()), kegg_ids)
        for kegg_id in kegg_ids:
            self.assertIsNotNone(reaction.get_compound(kegg_id))
            self.assertNotEqual(reaction.get_coeff(kegg_id), 0)
        self.assertIsNone(reaction.get_compound('C00003'))
        self.assertEqual(reaction.get_coeff('C00003'), 0)
        
    def test_reaction_balancing(self):
        kegg_id_to_coeff = {'C00011' : -1, 'C00001' : -1, 'C01353' : 1}
        reaction = Reaction(kegg_id_to_coeff)
        self.assertTrue(reaction.check_full_reaction_balancing())

        kegg_id_to_coeff = {'C00036' : -1, 'C00149' : 1}  # oxaloacetate = malate
        reaction = Reaction(kegg_id_to_coeff)
        self.assertAlmostEqual(reaction.check_half_reaction_balancing(), 2.0)

        kegg_id_to_coeff = {'C00031' : -1, 'C00469' : 2}  # missing two CO2
        reaction = Reaction(kegg_id_to_coeff)
        self.assertDictEqual(reaction._get_reaction_atom_balance(),
                             {'O': -4, 'C': -2, 'e-': -44})


    def test_gibbs_energy(self):
        kegg_id_to_coeff = {'C00002' : -1, 'C00001' : -1,
                            'C00008' :  1, 'C00009' :  1} # ATP + H2O = ADP + Pi
        reaction = Reaction(kegg_id_to_coeff)

        cc = ComponentContribution(pH=7.0, ionic_strength=0.1)
        dG0_prime, dG0_uncertainty = cc.dG0_prime(reaction)
        
        self.assertAlmostEqual(dG0_prime, -26.4, 1)
        self.assertAlmostEqual(dG0_uncertainty, 0.6, 1)

    def test_reduction_potential(self):
        kegg_id_to_coeff = {'C00036' : -1, 'C00149' : 1}  # oxaloacetate = malate
        reaction = Reaction(kegg_id_to_coeff)

        cc = ComponentContribution(pH=7.0, ionic_strength=0.1)
        E0_prime_mV, E0_uncertainty = cc.E0_prime(reaction)
        
        self.assertAlmostEqual(E0_prime_mV, -175.2, 1)
        self.assertAlmostEqual(E0_uncertainty, 5.3, 1)

    def test_mdf(self):
        sbtab_fname = os.path.join(settings.TEST_DIR, 'pathway_ethanol_SBtab.tsv')
        pp = Pathway.from_sbtab(sbtab_fname)
        mdf_res = pp.calc_mdf()

        self.assertAlmostEqual(mdf_res.mdf, 1.69, 2)


if __name__ == '__main__':
    unittest.main()
