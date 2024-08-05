# gMCSpy    
# ![Python versions](https://img.shields.io/pypi/pyversions/gMCSpy?logo=python) ![Downloads](https://img.shields.io/pypi/dm/gMCSpy) ![Tests](https://github.com/PlanesLab/gMCSpy/actions/workflows/tests.yml/badge.svg)


gMCSpy is a python package for the calculation of Genetic Minimal Cut sets (GMCS). In simple terms the idea is to take a metabolic model and calculate the genetic vulnerabilities that will render the biomass production impossible. This is done through a Mixed-Integer Linear problem (MILP) formultion and the use of a linear solver.  
The models must come from the [cobrapy](https://opencobra.github.io/cobrapy/) package and a linear solver must be installed. The package has been design to be used with [Gurobi](https://www.gurobi.com/),  [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) and [SCIP](https://scipopt.org/#scipoptsuite).


## Documentation

[Documentation](https://planeslab.github.io/gMCSpy/)

## Installation

Install gmcspy from pip

```bash
  pip install gmcspy
```
  
## Quick Start

### Calculating Genetic Minimal Interventions

[Documentation](https://planeslab.github.io/gMCSpy/)

### E. coli core

```python
#Read the model 
from pathlib import Path
from cobra.io import load_model
model = load_model("textbook")
```

To calculate all the GMCS of length 3 or less; using gurobi as solver

```python
#Calculate the genetic minimal cut sets
from gMCSpy import calculateGeneMCS

calculateGeneMCS(
        cobraModel=model,
        maxKOLength=3,
        solver='gurobi'
)

### Using CPLEX

calculateGeneMCS(
        cobraModel=model,
        maxKOLength=3,
        solver='cplex'
)
```

### Human-GEM (*Requires License Activation: CPLEX and Gurobi)
Using HUMAN GEM (v16) from [here](https://github.com/SysBioChalmers/Human-GEM/releases). 
```bash
mkdir data
curl -o data/Human-GEM16.mat --location --remote-header-name https://github.com/SysBioChalmers/Human-GEM/raw/v1.16.0/model/Human-GEM.mat

```

```python
#Read the model 
from pathlib import Path
from cobra.io import load_matlab_model

mini_mat_path = Path(".") / "./data/Human-GEM16.mat"
model = load_matlab_model(str(mini_mat_path.resolve()))
```

To calculate all the GMCS of length 3 or less; using gurobi as solver

```python
#Calculate the genetic minimal cut sets
from gMCSpy import calculateGeneMCS

calculateGeneMCS(
        cobraModel=model,
        maxKOLength=3,
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

- Carlos J. Rodriguez
- Naroa Barrena 
- Danel Olaverri-Mendizabal
- Idoia Ochoa
- Luis V. Valcárcel
- Francisco J. Planes
