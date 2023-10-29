from gMCSpy import calculateGeneMCS
from pathlib import Path
from cobra.io import load_matlab_model



mini_mat_path = (Path(".") / "C:/Users/crodr/Escritorio/Projects/GMCs/GMCS.py/data/Human-GEM-1.16.0_CultureMedia.mat")
model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
model_ecoli.id = "e_coli_core"



solutionDict = calculateGeneMCS(
        cobraModel=model_ecoli,
        maxKOLength=3,
        earlyStop=True,
        removeRelated=True,
        solver='gurobi',
        forceLength=False, 
        checkMCS=True,
        )
