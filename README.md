# gMCSpy

gMCSpy is a python package for the calculation of Genetic Minimal Cut sets (GMCS). In simple terms the idea is to take a metabolic model and calculate the genetic vulnerabilities that will render the biomass production impossible. This is done through a Mixed-Integer Linear problem (MILP) formultion and the use of a linear solver.  
The models must come from the [cobrapy](https://opencobra.github.io/cobrapy/) package and a linear solver must be installed. The package has been design to be used with [Gurobi](https://www.gurobi.com/),  [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) and [SCIP](https://scipopt.org/#scipoptsuite).


## Documentation

[Documentation](https://cjrodriguez98.github.io/GMCS.py/)

## Installation

Install gmcspy from pip

```bash
  pip install gmcspy
```
  
## Quick Start

Using HUMAN GEM (v14) from [here](https://github.com/SysBioChalmers/Human-GEM/releases). 

```python
#Read the model 
from pathlib import Path
from cobra.io import load_matlab_model

mini_mat_path = Path(".") / "./data/Human-GEM14.mat"
model = load_matlab_model(str(mini_mat_path.resolve()))
```

To calculate all the GMCS of length 3 or less, with a maximum of 200 GMCS and using Gurobi as the solver.

```python
#Calculate the genetic minimal cut sets
from gMCSpy import calculateGeneMCS

mcs = calculateGeneMCS(
        cobraModel=model,
        maxKOLength=3,
        maxNumberGMCS=200,
        solver='gurobi'
)
```
**Table with the results:**

| Order      | Solution                                                             |
|------------|----------------------------------------------------------------------|
| order1_0   | frozenset({'ENSG00000106105'})                                       |
| order1_1   | frozenset({'ENSG00000084774'})                                       |
| ...        | ...                                                                  |
| order3_159 | frozenset({'ENSG00000156471', 'ENSG00000185813', 'ENSG00000213930'}) |

## Authors

- Carlos Rodriguez
- Naroa Barrena 
- Luis Valcarcel
- Francisco Planes
