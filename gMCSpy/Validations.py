from .ProblemDefinitions import cobra
from .ProblemDefinitions import buildDictionaryFBAWithDeletions
from .ProblemDefinitions import buildDictionaryReactionsFBAWithDeletions
from .MinimalCutSetClass import logging
from .MinimalCutSetClass import tqdm
from. MinimalCutSetClass import prepareModelNutrientGeneMCS
from .Utilities import datetime
from .Utilities import os

import time

from joblib import Parallel, delayed


def checkGMCS(
    geneSolution,
    model: cobra.Model,
    isoFormSeparator=None,
):
    problem = buildDictionaryFBAWithDeletions(model, geneSolution, isoFormSeparator)
    problem.setSolver("gurobi")
    gurobi_default_params = {
        "Threads": 1,
        "OutputFlag": 0
        }
    [problemInterpreted, _] = problem.interpretProblem(parameters=gurobi_default_params)
    problemInterpreted.optimize()
    try:
        cutSetValue = problemInterpreted.ObjVal
    except:
        checkDict = {
        "isGMCS": 'infesible',
        "solution": geneSolution,
        "isCutSet": 'infesible',
        "cutSetValue": 'infesible',
        "isMinimal": 'infesible',
        "setGeneTested": 'infesible',
        "objVals": 'infesible',
        }
        return checkDict

    if problemInterpreted.ObjVal > 1e-6:
        isCutSet = False
    else:
        isCutSet = True

    isMinimal = []
    setGeneTested = []
    objVals = []
    for i in range(len(geneSolution)):
        ko = list(geneSolution)

        if len(ko) > 1:
            ko.remove(ko[i])
            problem = buildDictionaryFBAWithDeletions(model, ko, isoFormSeparator)
            problem.setSolver("gurobi")
            [problemInterpreted, interface] = problem.interpretProblem(parameters=gurobi_default_params)
            problemInterpreted.optimize()
            objVals.append(problemInterpreted.ObjVal)
            if problemInterpreted.ObjVal > 0:
                isMinimal.append(True)
            else:
                isMinimal.append(False)

        elif len(ko) == 1 and isCutSet:
            isMinimal.append(True)
            objVals.append(0)

        setGeneTested.append(ko)

    areAllMinimal = all(isMinimal)
    both = [isCutSet, areAllMinimal]
    isGMCS = all(both)
    checkDict = {
        "isGMCS": isGMCS,
        "solution": geneSolution,
        "isCutSet": isCutSet,
        "cutSetValue": cutSetValue,
        "isMinimal": isMinimal,
        "setGeneTested": setGeneTested,
        "objVals": objVals,
    }
    return checkDict

def checkFromReaction(
    genes: frozenset,
    gDict: dict,
    model: cobra.Model
):
    reactions = []

    for key in gDict.keys():
        if key.issubset(genes):
            [reactions.append(reaction) for reaction in gDict[key]]

    reactions = list(set(reactions))
    reactions = [model.reactions[reaction].id for reaction in reactions]
    problem = buildDictionaryReactionsFBAWithDeletions(model, reactions)
    problem.setSolver("gurobi")
    gurobi_default_params = {
        "Threads": 1,
        "OutputFlag": 0
        }
    [problemInterpreted, _] = problem.interpretProblem(parameters=gurobi_default_params)
    problemInterpreted.optimize()
    cutSetValue = problemInterpreted.ObjVal
    if problemInterpreted.ObjVal > 1e-6:
        isCutSet = False
    else:
        isCutSet = True

    isMinimal = []
    setGeneTested = []
    objVals = []
    for i in range(len(genes)):
        ko = list(genes)
        if len(ko) > 1:
            ko.remove(ko[i])
            # get the reactions for the new gene set
            reactions = []
            setKo = frozenset(ko)
            if setKo in gDict:
                [reactions.append(reaction) for reaction in gDict[setKo]]
            else:
                for key in gDict.keys():
                    if key.issubset(setKo):
                        [reactions.append(reaction) for reaction in gDict[key]]
            if len(reactions) == 0:
                continue
            problem = buildDictionaryReactionsFBAWithDeletions(model, reactions)
            problem.setSolver("gurobi")
            [problemInterpreted, interface] = problem.interpretProblem(parameters=gurobi_default_params)
            problemInterpreted.optimize()
            objVals.append(problemInterpreted.ObjVal)
            if problemInterpreted.ObjVal > 0:
                isMinimal.append(True)
            else:
                isMinimal.append(False)

        elif len(ko) == 1 and isCutSet:
            isMinimal.append(True)
            objVals.append(0)

        setGeneTested.append(ko)

    areAllMinimal = all(isMinimal)
    both = [isCutSet, areAllMinimal]
    isGMCS = all(both)
    checkDict = {
        "isGMCS": isGMCS,
        "solution": genes,
        "isCutSet": isCutSet,
        "cutSetValue": cutSetValue,
        "isMinimal": isMinimal,
        "setGeneTested": setGeneTested,
        "objVals": objVals,
    }
    return checkDict

# Build parrallel function to check all solutions in a list

def checkGMCSParallel(knockoutList, model, workers, isoFormSeparator=None, log_location=None, isNutrient=False, handler=None):
    
    if len(knockoutList) == 0:
        print("No solutions to check")
        return
    
    if handler:
        loggerMain = logging.getLogger(handler)
    else:
        loggerMain = None
    
    timeStart = time.time()
    
    if isNutrient:
        model = prepareModelNutrientGeneMCS(model, [])
        
    # Set the location and filename for the log file
    if log_location is None:
        folderLogPath = "checks/check_"
        log_location = folderLogPath + model.id + datetime.now().strftime("%Y_%m_%d-%I_%M_%S") + ".txt"
        os.makedirs(os.path.dirname(log_location), exist_ok=True)
        print('saving check results to: ' + log_location)
    else:
        log_location = log_location + "checkedSolutions.csv"
        os.makedirs(os.path.dirname(log_location), exist_ok=True)
        print('saving check results to: ' + log_location)
    # Create a logger
    logger = logging.getLogger('checks')
    logger.setLevel(logging.DEBUG)

    # Create a file handler
    file_handler = logging.FileHandler(log_location)
    file_handler.setLevel(logging.DEBUG)

    # Customize the log message format
    formatter = logging.Formatter('%(message)s')  # Only includes the log message

    # Set the formatter for the file handler
    file_handler.setFormatter(formatter)

    # Add the file handler to the logger
    logger.addHandler(file_handler)

    # Write the first line to the log
    logger.debug('model,isGMCS,completeDict')
    
    # zip al parameters into a list repeating model and isoFormSeparator
    listModels = [model] * len(knockoutList)
    listIsoFormSeparator = [isoFormSeparator] * len(knockoutList)
    zipppedParameters = zip(knockoutList, listModels, listIsoFormSeparator)

    # Redirect stdout to a null device
    with open(os.devnull, "w") as devnull:
        original_stdout = os.dup(1)  # Save the original stdout
        os.dup2(devnull.fileno(), 1)  # Replace stdout with the null device

    results_iter = Parallel(n_jobs=workers)(delayed(checkGMCSChild)(zippedParameters) for zippedParameters in tqdm(zipppedParameters, total=len(knockoutList), 
                                                                                                              desc="Checking GMCSs", unit="GMCS"))

    # Restore the original stdout
    os.dup2(original_stdout, 1)

    percentageOfTrue = []
    for result in results_iter:
        logger.debug(f"{model.id},{result['isGMCS']},{result}")
        percentageOfTrue.append(result['isGMCS'])
    
    endTime = time.time()
    timeElapsed = endTime - timeStart
    if loggerMain:
        loggerMain.info(f"checks: {timeElapsed}")
    print('Percentage of True GMCS: ' +  str((sum(percentageOfTrue)/len(percentageOfTrue))*100))

    logger.removeHandler(file_handler)
    file_handler.close()
    

def checkGMCSChild(zippedParameters):
    knockout = zippedParameters[0]
    model = zippedParameters[1]
    isoFormSeparator = zippedParameters[2]  
    checkDict = checkGMCS(knockout, model, isoFormSeparator)
    return checkDict

