# equilibrator-api

A command-line API with minimal dependencies for calculation of standard thermodynamic potentials of biochemical reactions using the data found on [eQuilibrator](http://equilibrator.weizmann.ac.il/). 

# Current Features 

* Example scripts for singleton and bulk calculations.
* Calculation of standard Gibbs potentials of reactions.
* Calculation of standard reduction potentials of half-cells.

# Example Usage

Import the API and create an instance. Creating the EquilibratorAPI class instance reads all the data that is used to calculate thermodynamic potentials of reactions. 

```python
from equilibrator_api import EquilibratorAPI, Reaction
eq_api = EquilibratorAPI()  # loads data
```

You can parse a reaction from a KEGG-style reaction string. The example given is ATP hydrolysis to ADP and inorganic phosphate.

```python
rxn_str = "C00002 + C00001 = C00008 + C00009"
rxn = Reaction.parse_formula(rxn_str)
```

We highly recommend that you check that the reaction is atomically balanced (conserves atoms) and charge balanced (redox neutral). We've found that it's easy to accidentally write unbalanced reactions in this KEGG-style format and so we always check ourselves.

```python
if not rxn.check_full_reaction_balancing():
	print '%s is not balanced', rxn
```

Now we know that the reaction is "kosher" and we can safely proceed to calculate the standard change in Gibbs potential due to this reaction.

```python
# You control the pH and ionic strength!
# ionic strength is in Molar units.
dG0_prime, dG0_uncertainty = eq_api.dG0_prime(
        rxn, pH=6.5, ionic_strength=0.2) 
print u"dG0' = %.1f \u00B1 %.1f kJ/mol\n" % (
        dG0_prime, dG0_uncertainty)
```

You can also calculate the [reversibility index](https://doi.org/10.1093/bioinformatics/bts317) for this reaction. 

```python
ln_RI = rxn.reversibility_index(pH=6.5, ionic_strength=0.2)
print u'ln(Reversibility Index) = %.1f\n' % ln_RI
```

The reversibility index is a measure of the degree of the reversibility of the reaction that is normalized for stoichiometry. If you are interested in assigning reversibility to reactions we recommend this measure because 1:2 reactions are much "easier" to reverse than reactions with 1:1 or 2:2 reactions. You can see the paper linked above for more information.

# dependencies:
- python 2.7
- numpy (preferably >= 1.12.0)

