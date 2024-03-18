import straindesign as sd
import cobra
import time
import pandas as pd
from pathlib import Path
from cobra.io import read_sbml_model, load_matlab_model
from gMCSpy import calculateGeneMCS
import os

numWorkers = 16

def testModel(model, model_name, solver, iter):
    path = "/scratch/a905383/inigoObj"
    startTime = time.time()
    solutionDict = calculateGeneMCS(
        cobraModel=model,
        maxKOLength=3,
        maxNumberGMCS=1e6,
        earlyStop=True,
        removeRelated=True,
        verbose=0,
        solver=solver,
        numWorkers = numWorkers,
        timeLimit = 1e4,
    )
    endTime = time.time()
    totalTime = endTime - startTime
    interventions = []
    for order, ko in solutionDict.items():
        ko_string = "__".join(list(ko["solution"]))
        interventions.append(ko_string)
    filename = f"GMCs_python_{model_name}_{solver}_{iter}.csv"
    completePath = path + filename
    df = pd.DataFrame({"gMCS": interventions})
    df.to_csv(completePath)
    return [model_name, totalTime, len(solutionDict), filename]

def benchSuite(models, solvers, iterations, logger):
    for z in enumerate(models):
        for solver in solvers:
                for i in range(iterations):
                    model = z[1]
                    result = testModel(model, model.id, solver, i)
                    logger.debug(
                        f"{result[0]},gmcspy,{solver},{result[1]},{result[2]}"
                    )
              

######################
# Load the models
######################
mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/e_coli_core.mat")
model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
model_ecoli.id = "ecoli_core"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/Recon3D.mat")
model_recon3d = load_matlab_model(str(mini_mat_path.resolve()))
model_recon3d.id = "recon3D"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/iJN1463.mat")
model_iJN = load_matlab_model(str(mini_mat_path.resolve()))
model_iJN.id = "iJN1463"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/iML1515.mat")
model_iML = load_matlab_model(str(mini_mat_path.resolve()))
model_iML.id = "iML1515"


mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/Human-GEM-1.16.0_CultureMedia.mat")
model_Human16Media = load_matlab_model(str(mini_mat_path.resolve()))
model_Human16Media.id = "Human-GEM-1.16.0_CultureMedia"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/Human-GEM-1.16.0_raw.mat")
model_Human16 = load_matlab_model(str(mini_mat_path.resolve()))
model_Human16.id = "Human-GEM-1.16.0"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/Human-GEM-1.17.mat")
model_Human17 = load_matlab_model(str(mini_mat_path.resolve()))
model_Human17.id = "Human-GEM-1.17"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/yeast8.mat")
model_Yeast = load_matlab_model(str(mini_mat_path.resolve()))
model_Yeast.id = "Yeast-GEM-8.7"



import logging
import os

os.chdir('/scratch/a905383/inigoObj/')

# Set the location and filename for the log file
# Set the location and filename for the log file
log_location = (
    "/scratch/a905383/inigoObj/results_python_"
    + time.strftime("%Y%m%d-%H%M%S")
    + ".csv"
)

# Create a logger
logger = logging.getLogger("my_logger")
logger.setLevel(logging.DEBUG)

# Create a file handler
file_handler = logging.FileHandler(log_location)
file_handler.setLevel(logging.DEBUG)

# Customize the log message format
formatter = logging.Formatter("%(message)s")  # Only includes the log message

# Set the formatter for the file handler
file_handler.setFormatter(formatter)

# Add the file handler to the logger
logger.addHandler(file_handler)

# Write the first line to the log
logger.debug("model,algirithm,solver,time,size")

models = [
          model_ecoli,
          model_Yeast,
          model_Human16,
          model_Human17,
          model_iML,
          model_iJN,
          model_Human16Media,
          model_recon3d,
         ]

numIter = 10

solvers = [
           'cplex',
           'gurobi'
        ]

res = benchSuite(models, solvers, numIter, logger)


