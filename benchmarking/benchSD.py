import straindesign as sd
import cobra
import time
import pandas as pd
from pathlib import Path
from cobra.io import read_sbml_model, load_matlab_model
from gMCSpy import calculateGeneMCS
import os

numWorkers = 16

def benchmarkModelSD(model, solver, approach, max_cost, iter):
    # get model objective function
    cobra.Configuration().solver = solver
    cobra.Configuration().processes = numWorkers
    print(f"Solver: {solver}")
    print(cobra.Configuration())
    print(f"When the model\'s solver is \'{model.solver.configuration}', StrainDesign selects {sd.select_solver()}.")
    startTime = time.time()
    objectiveVector = [r.objective_coefficient for r in model.reactions]
    indexC = objectiveVector.index(1)
    objectiveName = model.reactions[indexC].id
    
    try:
        # Building suppress module
        threshold = 1e-3
        constraint = f'{objectiveName}>={threshold}'
        module_suppress  = sd.SDModule(model,'suppress', constraints=constraint)
        sols = sd.compute_strain_designs(model,
                                    sd_modules = module_suppress,
                                    solution_approach = approach,
                                    gene_kos = True,
                                    max_cost = max_cost,
                                    )
    except Exception as e:
        # print the error that made the code crash
        print(e)
        return 'error', 'error'
    endTime = time.time()
    timeTaken = endTime - startTime
    #Transforming the solution into a dataframe
    solsList = []
    for gene in sols.gene_sd:
        geneList = list(gene.keys())
        solsList.append(geneList)
    filename = f"{model.id}_{solver}_{approach}_{iter}.tsv"
    with open(filename, "w") as f:
        for gene in solsList:
            f.write("%s ,\n" % gene)
            
    return timeTaken, len(sols.gene_sd)

def benchSuiteSD(models, solvers, approaches, iterations, logger):
    for z in enumerate(models):
        for solver in solvers:
            for approach in approaches:
                for i in range(iterations):
                    model = z[1]
                    result = benchmarkModelSD(model, solver, approach, max_cost=3, iter=i)
                    logger.debug(
                        f"{model.id},straindesign,{solver},{result[0]},{result[1]}"
                    )
         
                    # Getting all memory using os.popen()
                    total_memory, used_memory, free_memory = map(
                        int, os.popen('free -t -m').readlines()[-1].split()[1:])
 
                    # Memory usage
                    print("RAM memory % used:", round((used_memory/total_memory) * 100, 2))
                    
           

######################
# Load the models
######################
mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/e_coli_core.mat")
model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
model_ecoli.id = "ecoli_core"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/iEK1008.mat")
model_iEK = load_matlab_model(str(mini_mat_path.resolve()))
model_iEK.id = "iEK1008"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/iJN1463.mat")
model_iJN = load_matlab_model(str(mini_mat_path.resolve()))
model_iJN.id = "iJN1463"

mini_mat_path = (Path(".") / "/scratch/a905383/gmcspy/data/iML1515.mat")
model_iML = load_matlab_model(str(mini_mat_path.resolve()))
model_iML.id = "iML1515"





import logging
import os

os.chdir('/scratch/a905383/resultsBenchSDGM')

# Set the location and filename for the log file
# Set the location and filename for the log file
log_location = (
    "/scratch/a905383/gmcspy/matlab_codes/results/results_SD_"
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
          model_iML,
          model_iJN,
         ]




numIter = 10

approaches = ['best', 'any', 'populate']
approaches = ['best']


solvers = [
           'gurobi',
        ] 

res = benchSuiteSD(models, solvers, approaches, numIter, logger)



models = [
          model_ecoli,
          model_iML,
          model_iJN,
         ]


solvers = [
           'cplex',
        ]

res = benchSuiteSD(models, solvers, approaches, numIter, logger)



