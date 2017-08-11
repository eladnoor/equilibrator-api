# equilibrator-api

A command-line API with minimal dependencies for calculation of standard thermodynamic potentials of biochemical reactions using the data found on [eQuilibrator](http://equilibrator.weizmann.ac.il/). 

# Current Features 

* Example scripts for singleton and bulk calculations.
* Calculation of standard Gibbs potentials of reactions.
* Calculation of standard reduction potentials of half-cells.

# Example Usage

Import the API and create an instance.  

```python
from equilibrator_api import EquilibratorAPI, Reaction
eq_api = EquilibratorAPI()
```

Now you can parse a reaction from a KEGG-style reaction string. The example given is ATP hydrolysis to ADP and Pi

```python
rxn_str = "C00002 + C00001 = C00008 + C00009"
rxn = Reaction.parse_formula(rxn_str)
```

Finally you can calculate various things for this reaction

```python
if not rxn.check_full_reaction_balancing():
	print '%s is not balanced', rxn

# You control the pH and ionic strength!
# ionic strength is in Molar units.
dG0_prime, dG0_uncertainty = eq_api.dG0_prime(
        rxn, pH=6.5, ionic_strength=0.2) 
print u"dG0' = %.1f \u00B1 %.1f kJ/mol\n" % (
        dG0_prime, dG0_uncertainty)

# reversibility index is a measure of reaction reversibility
# that accounts for stoichiometry. 
# https://doi.org/10.1093/bioinformatics/bts317
ln_RI = rxn.reversibility_index(pH=6.5, ionic_strength=0.2)
print u'ln(Reversibility Index) = %.1f\n' % ln_RI
```

# dependenceis:
- python 2.7
- numpy (preferably >= 1.12.0)

