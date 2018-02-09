![eQuilibrator logo](img/logo.svg)
<img src="img/logo.svg">

A command-line API with minimal dependencies for calculation of standard thermodynamic potentials of biochemical reactions using the data found on [eQuilibrator](http://equilibrator.weizmann.ac.il/).
Does not require any network connections.

## Current Features

* Example scripts for singleton and bulk calculations.
* Calculation of standard Gibbs potentials of reactions (together with confidence intervals).
* Calculation of standard reduction potentials of half-cells.

To access more advanced features, such as adding new compounds that are not available in the KEGG database,
try using our full-blown [Component Contribution](https://github.com/eladnoor/component-contribution)
package.

## Cite us

If you plan to use results from equilibrator-api in a scientific publication,
please cite our paper:

Noor E, Haraldsdóttir HS, Milo R, Fleming RMT. Consistent estimation of Gibbs energy using component contributions. PLoS Comput Biol. 2013;9: e1003098.

## Example Usage

Import the API and create an instance. Creating the EquilibratorAPI class instance reads all the data that is used to calculate thermodynamic potentials of reactions.

```python
from equilibrator_api import ComponentContribution, Reaction
eq_api = ComponentContribution(pH=6.5, ionic_strength=0.2)  # loads data
```

You can parse a reaction from a KEGG-style reaction string. The example given is ATP hydrolysis to ADP and inorganic phosphate.

```python
rxn_str = "C00002 + C00001 = C00008 + C00009"
rxn = Reaction.parse_formula(rxn_str)
```

We highly recommend that you check that the reaction is atomically balanced (conserves atoms) and charge balanced (redox neutral). We've found that it's easy to accidentally write unbalanced reactions in this KEGG-style format and so we always check ourselves.

```python
if not rxn.check_full_reaction_balancing():
	print('%s is not balanced' % rxn)
```

Now we know that the reaction is "kosher" and we can safely proceed to calculate the standard change in Gibbs potential due to this reaction.

```python
# You control the pH and ionic strength!
# ionic strength is in Molar units.
dG0_prime, dG0_uncertainty = eq_api.dG0_prime(rxn)
print("dG0' = %.1f \u00B1 %.1f kJ/mol\n" % (dG0_prime, dG0_uncertainty))
```

You can also calculate the [reversibility index](https://doi.org/10.1093/bioinformatics/bts317) for this reaction.

```python
ln_RI = rxn.reversibility_index()
print('ln(Reversibility Index) = %.1f\n' % ln_RI)
```

The reversibility index is a measure of the degree of the reversibility of the reaction that is normalized for stoichiometry. If you are interested in assigning reversibility to reactions we recommend this measure because 1:2 reactions are much "easier" to reverse than reactions with 1:1 or 2:2 reactions. You can see the paper linked above for more information.

## Dependencies:
- python 2.7 or 3.6
- numpy (preferably >= 1.13.3)
### Only if running MDF analysis:
- scipy (preferably >= 1.0.0)
- PuLP (1.6.8)
- tablib (0.12.1)
- matplotlib (2.0.2)
- pandas (0.20.3)

