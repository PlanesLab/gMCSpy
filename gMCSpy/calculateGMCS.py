from .ProblemDefinitions import cobra
from .ProblemInterpreters import np


from .ProblemInterpreters import beartype
from .ProblemDefinitions import List
from .ProblemDefinitions import scipy
from .ProblemDefinitions import buildDictionaryMCSGeneProblem
from .ProblemDefinitions import buildDictionaryMCSGeneTargetedProblem

from .Validations import checkGMCSParallel

from .calculateMCS import calculateMCS

import re
from ast import And, BoolOp, Or, Name, parse

from collections import defaultdict
import warnings
import time

from .MinimalCutSetClass import MinimalCutSetProblem
from .MinimalCutSetClass import uuid
from .MinimalCutSetClass import tqdm
from .MinimalCutSetClass import logging
from .MinimalCutSetClass import prepareModelNutrientGeneMCS

from datetime import datetime
import pandas as pd

from .Utilities import createLog


def calculateGeneMCS(cobraModel, **kwargs):
    '''
    **Calculate the genetic minimal cut sets of a model**
    Calculates the gMCS of a metabolic network using a mixed integer linear programming (MILP) formulation.
    Among all genes or with a given subset of genes, the function finds the minimal set of genes
    to stop biomass production.

    Required arguments:
        - **cobraModel:** A metabolic model in cobrapy format. See `cobra models <https://cobrapy.readthedocs.io/en/latest/io.html#Reading- and- Writing- Models>`_.

    Optional arguments:
        - **Solver:** Solver to be used. Current compatibility includes cplex, gurobi and SCIP. Default value is gurobi.
        - **maxKOLength :** Maximum length of gMCS to be found. (i.e. number of genes in the gMCS). Default value is 3
        - **maxNumberGMCS :** Maximum number of gMCS to be found. Default value is 1e6
        - **verbose:** If 0, the solver will not print any information. Verbose  = 1, log will be printed into the console. Verbose = 2, log will be printed into the console and a log file will be created. Default value is 0.
        - **earlyStop :** If True if optimality returned by the solver the solver will stop without running a second time. Default value is True.
        - **removeRelated:** If True, the solver will remove the gMCS that are a subset of another gMCS from the solution space. Default value is True.
        - **targetKOs:** List of genes that you want to be part of the gMCS solutions. The list is expected to be of gene names. Default value is an empty list, which means that all genes will be considered. The gMCS solutions will contain at least one of the genes in the list.
        - **geneSubset:** List of genes to be considered in the gMCS problem. The list is expected to be of gene names. Default value is an empty list, which means that all genes will be considered. This is useful when you want to reduce the search space of the gMCS problem. e.g. geneSubset = ["gene1", "gene2", "gene3"] will only consider genes gene1, gene2 and gene3 as possible intervention points. Only ["gene1"], ["gene2"], ["gene3"], ["gene1", "gene2"], ["gene1", "gene3"], ["gene2", "gene3"] and ["gene1", "gene2", "gene3"] could be solutions. This not compatible with targetKOs.
        - **isNutrient:** If True, the solver will consider the exchange reactions as nutrients, assigning artificial genes to exchange reactions to evaluate what metabolites and genes are essential. Default value is False.
        - **exchangeReactions:** List of exchange reactions to be considered as nutrients. The list is expected to be of reaction names. Default value is None.
        - **saveGMatrix:** If True, the solver will save the G matrix in a csv file. Default value is False.
        - **checkSolutions:** If True, the solver will check if the solutions are valid. Default value is False.
        - **saveSolutions:** If True, the solver will save the solutions in a csv file. Default value is True.
        - **timelimit:** Maximum time allowed for the solver to find a solution. Default value is 1e75.
        - **forceLength:** If True, the solver will find MCS of length 1, 2, 3, ..., MAX_LENGTH_MCS. If False, the solver will find MCS of length MAX_LENGTH_MCS. Default value is True for cplex and false for gurobi.
        - **numWorkers:** Number of workers to be used for parallel processing. Default value is 0. The solvers cplex and gurobi already use parallel processing.      
        - **isoformSeparator:** If not None, the isoforms of the genes will merge into one gene. Default value is None.
    '''
    kwargs["solver"] = kwargs.get("solver", "gurobi")
    name = (
        "gMCSpy_"
        + cobraModel.id
        + "_"
        + kwargs["solver"]
        + "_"
        + datetime.now().strftime("%Y%m%d-%H%M%S")
    )
    
    path = "logs/" + name + "/process/"
    createLog(path)
    path = "logs/" + name
    name = ("process")
    handler = name + ".log"
    logging = initiateLogger(name + ".log", path + "/process/" + name + ".log")
    allTime = time.time()

    maxKOLength = kwargs.get('maxKOLength', 3)
    maxNumberGMCS = kwargs.get('maxNumberGMCS', 1e6)
    verbose = kwargs.get('verbose', 0)
    earlyStop = kwargs.get("earlyStop", True)
    removeRelated = kwargs.get("removeRelated", False)
    isoformSeparator = kwargs.get("isoformSeparator", None)
    targetKOs = kwargs.get("targetKOs", None)
    geneSubset = kwargs.get("geneSubset", None)
    isNutrient = kwargs.get("isNutrient", False)
    exchangeReactions = kwargs.get("exchangeReactions", None)
    saveGMatrix = kwargs.get("saveGMatrix", False)
    checkSolutions = kwargs.get("checkSolutions", False)
    saveSolutions =  kwargs.get('saveSolutions', True)
    timeLimit = kwargs.get("timeLimit", 1e4)
    solver = kwargs["solver"]
    
    if kwargs.get("forceLength") is None:
        if solver == "gurobi":
            forceLength = False
        elif solver == "cplex":
            forceLength = True
        elif solver == "scip":
            forceLength = True   
    else:
        forceLength = kwargs.get("forceLength", False)
    
    kwargs.pop("maxKOLength", None)
    kwargs.pop("maxNumberGMCS", None)
    kwargs.pop("verbose", None)
    kwargs.pop("earlyStop", None)
    kwargs.pop("removeRelated", None)
    kwargs.pop("isoformSeparator", None)
    kwargs.pop("targetKOs", None)
    kwargs.pop("geneSubset", None)
    kwargs.pop("isNutrient", None)
    kwargs.pop("exchangeReactions", None)
    kwargs.pop("saveGMatrix", None)
    kwargs.pop("checkMCS", None)
    kwargs.pop('saveSolutions', None)
    kwargs.pop("forceLength", None)
    kwargs.pop("timeLimit", None)

    numWorkers = kwargs.get("numWorkers", 0)
       
    if isNutrient:
        cobraModel = prepareModelNutrientGeneMCS(cobraModel, exchangeReactions)

    logging.info(
        f"model: {cobraModel.id}, maxKOLength: {maxKOLength}, maxNumberGMCS: {maxNumberGMCS}, verbose: {verbose}, forceLength: {forceLength}, earlyStop: {earlyStop}, removeRelated: {removeRelated}, solver: {kwargs['solver']}"
    )
    startTime = time.time()
    if isoformSeparator is not None:
        gObject = buildGMatrix(
            cobraModel,
            maxKOLength,
            isoformSeparator=isoformSeparator,
            verbose=0,
            **kwargs,
        )
    else:
        gObject = buildGMatrix(cobraModel, maxKOLength, verbose=0, **kwargs)
    endTime = time.time()
    logging.info("GMatrix: " + str(endTime - startTime))
    
    gDict = gObject["gDict"]
    gMatrix = gObject["gMatrix"]
    relationships = gObject["relationships"]
    numberNewGenesByKO = gObject["numberNewGenesByKO"]
    
    # Filter the G matrix to keep only genes that are in the geneSubset
    if geneSubset:
        #Make a set of the geneSubset
        geneSet = set(geneSubset)
        # Extract all possible interventions
        gDictKeys = list(gDict.keys())
        # Remove the interventions that are not in the geneSubset
        for key in gDictKeys:
            if not key.issubset(geneSet):
                gDict.pop(key)

        gMatrix = createSparseMatrix(gDict, cobraModel.reactions)
        [relationships, numberNewGenesByKO] = relatedRows(gDict, mergeIsforms(isoformSeparator))
   
    if saveGMatrix:        
        gMatrix_Strings = []
        for ko, reactions  in gDict.items():
            string_ko = "@@".join(list(ko))
            string_reactions = "@@".join(cobraModel.reactions[x].id for x in reactions)
            row = string_ko + " -> " + string_reactions
            gMatrix_Strings.append(row)
        filename = f"Gmat_python_{cobraModel.id}_{solver}.csv"
        gPath = path + "/GMatrix/" 
        createLog(gPath)
        completePath = gPath + filename
        df = pd.DataFrame({'gMatrix': gMatrix_Strings})
        df.to_csv(completePath)


    # Initialize the genes to be knocked out, e.g. genes that you want to be part of the gMCS
    if targetKOs is None:
        targetKOs = []

    # Make sure that targetKOs is a list
    if not isinstance(targetKOs, list):
        raise ValueError("targetKOs must be a list")

    # Initialize the subset of genes to be studied , e.g. the space of genes to calculate the gMCS
    if geneSubset is None:
        geneSubset = []

    # Make sure that geneSubset is a list
    if not isinstance(geneSubset, list):
        raise ValueError("geneSubset must be a list")

    # List of possible interventions
    gDictKeys = list(gDict.keys())

    # If there are reactions to be knocked out, then the problem is a targeted problem
    if targetKOs:
        # Transform each element of the targetKOs into a frozenset to be able to compare it with the keys of the gDict
        targetKOs = [frozenset([gene]) for gene in targetKOs]
        # Check that the targetKOs is in the model
        for gene in targetKOs:
            isInModel = []
            for intervention in gDictKeys:
                if gene.issubset(intervention):
                    isInModel.append(intervention)
                    break
            if not isInModel:
                raise ValueError("One or more genes in targetKOs are not in the model")

        # Check that geneSubset is not an empty list
        if geneSubset:
            # Transform each element of the geneSubset into a frozenset to be able to compare it with the keys of the gDict
            geneSubset = [frozenset([gene]) for gene in geneSubset]
            # Check that the geneSubset is in the model
            for gene in geneSubset:
                isInModel = []
                for intervention in gDictKeys:
                    if gene.issubset(intervention):
                        isInModel.append(intervention)
                        break
                if not isInModel:
                    raise ValueError(
                        "One or more genes in geneSubset are not in the model"
                    )
        else:
            geneSubset = []
            for intervention in gDictKeys:
                if len(intervention) == 1:
                    geneSubset.append(intervention)

        geneSubset = list(set(geneSubset).union(set(targetKOs)))

        # Build the problem object
        startTime = time.time()
        problem = buildDictionaryMCSGeneTargetedProblem(
            cobraModel,
            gDict,
            gMatrix,
            relationships=relationships,
            geneSubset=geneSubset,
            targetKOs=targetKOs,
            maxKOLength=maxKOLength,
            forceLength=forceLength,
            numberNewGenesByKO=numberNewGenesByKO,
            verbose=verbose,
        )
        endTime = time.time()
        logging.info(
            "Time taken to build problem: " + str(endTime - startTime) + " seconds"
        )

        gurobi_default_params = {
            "MIPFocus": 1,
            "TimeLimit": max(1, timeLimit),
            "Threads": numWorkers,
            "PoolSolutions": 10,
            "PoolGap": 0.1,
            "PoolSearchMode": 1,
            "Cutoff": maxKOLength,
            "Presolve": 1,
            "PreSOS1Encoding": 2,
            "Cuts": 2,
        }

    else:
        startTime = time.time()
        problem = buildDictionaryMCSGeneProblem(
            cobraModel,
            gDict,
            gMatrix,
            relationships=relationships,
            maxKOLength=maxKOLength,
            forceLength=forceLength,
            numberNewGenesByKO=numberNewGenesByKO,
            verbose=verbose,
            logPath=path
        )
        endTime = time.time()
        logging.info("BuildProblem: " + str(endTime - startTime))

            
        gurobi_default_params = {
            "MIPFocus": 3,
            "TimeLimit": max(1, timeLimit),
            "Threads": numWorkers,
            "PoolSolutions": 10000,
            "PoolGap": 0.1,
            "PoolSearchMode": 1,
            "Cutoff": maxKOLength,
            "Presolve": 1,
            "PreSOS1Encoding": 2,
            "Cuts": 2,
            #"OptimalityTol": 1e-9,
        }

        if forceLength:
            gurobi_default_params['MIPFocus'] = 3
            gurobi_default_params['TimeLimit'] = 500
            
        '''
            gurobi_default_params["MIPFocus"] = 1
        '''
    problem.setCharacteristicParameter(maxKOLength)
    problem.setModelName(cobraModel.id)
    problem.setLogPath(path)
    # Set the solver
    
    problem.setSolver(solver)

    cplex_default_params = {
        "mip.tolerances.integrality": 1e-5,
        "emphasis.mip": 4,
        "timelimit": max(10, timeLimit),
        "threads": numWorkers,
        "preprocessing.presolve": 1,
        "mip.limits.populate": 1000,
        "mip.pool.relgap": 0.1,
        "mip.tolerances.uppercutoff": maxKOLength,
        "mip.pool.intensity": 3,
        "mip.strategy.probe": 2,
        "mip.strategy.variableselect": 4,
        "preprocessing.sos1reform": 1,
    }


    if solver == "gurobi":
        params = gurobi_default_params
    elif solver == "cplex":
        params = cplex_default_params
    else:
        params = {}
    # Interpret the problem, depending on the solver each problem has to be interpreted differently into the solver's interface
    optimization = problem.interpretProblem(verbose=2, parameters=params)

    # Import the common problem methods
    problemMethods = MinimalCutSetProblem()

    # Get the objective variables, they are the variables that we are interested in
    knockouts = list(gDict.keys())
    zpVariables = [
        key for key, data in problem.getVariables().items() if "zp_" in data["name"]
    ]
    solutionVariables = {"indices": zpVariables, "names": knockouts}
    # Set the bounds for the solver, e.g. maximum number of solutions, maximum numbers in a solution, or if we want to force the step by step calculation etc.
    bounds = {
        "MAX_SOLUTION_LENGTH": maxKOLength,
        "MAX_SOLUTIONS": maxNumberGMCS,
        "forceLength": forceLength,
    }

    # Solve the problem iteratively
    solutionDict = problemMethods.iterateOverSolver(
        problem=problem,
        problemInterpreted=optimization,
        solver=solver,
        solutionVariables=solutionVariables,
        bounds=bounds,
        mergeSolutions=True,
        earlyStop=earlyStop,
        removeRelated=removeRelated,
        verbose=2,
        handler=handler,
        iterationProgress=True,
    )
    allEndTime = time.time()

    logging.info("totalTime: " + str(allEndTime - allTime))
    
    if saveSolutions:
        solutionPath = path + '/solutions/'
        solutionFileName = 'solutions.log'
        createLog(solutionPath)
        loggingSolutions = initiateLogger(solutionFileName, solutionPath + solutionFileName)
        loggingSolutions.info('gene,reactions')
        geneSolutionList = []
        
        for order, ko in solutionDict.items():
            reactions_log_list = []
            for gene in ko['solution']:
                reactions = gDict[frozenset([gene])]
                reactions_log_list.append([cobraModel.reactions[x].id for x in reactions])
            reactions = ','.join(list(set().union(*reactions_log_list)))
            ko_log = ','.join(list(ko['solution']))         
            loggingSolutions.info(f'[{ko_log}]:[{reactions}]')

        solutionHandler = loggingSolutions.handlers[0]
        loggingSolutions.removeHandler(solutionHandler)
        solutionHandler.close()
               
    if checkSolutions:
        geneSolutionList = []
        for key, items in solutionDict.items():
            geneSolutionList.append(items['solution'])
        
        checkGMCSParallel(geneSolutionList, cobraModel, 16, isoFormSeparator=isoformSeparator, log_location=path + "/checks/", handler=handler)

    handler = logging.handlers[0]
    logging.removeHandler(handler)
    handler.close()
    return solutionDict

def buildGMatrix(
    cobraModel: cobra.core.Model,
    maxKOLength: int = 1e6,
    name: str = None,
    isoformSeparator: str = None,
    verbose: int = 0,
    **kwargs,
):
    """

    **Build the G matrix, G dict and the relationships between perturbations**

    The G matrix is a matrix with the perturbations of the genes in the model, that render a reaction infeasible.
    Each row of the matrix corresponds to a perturbation, and each column corresponds to a reaction.
    This matrix is used to calculate the genetic minimal cut sets of the model.

    :param cobraModel: The cobra model that represents the metabolism of an organism.
    :param name: The name of the model. If None, the id of the model will be used.
    :param separateIsoform: If not None, the isoforms of the genes will merge into one gene.
    :param verbose: The level of verbosity of the function.

    :return: A List with the G matrix, the G dictionary and the relationships between perturbations.


    """
    
    maxKOLength = kwargs.get('maxKOLength', 3)
    dictManager = mergeIsforms(isoformSeparator)
        
    modelReactions = [reaction.id for reaction in cobraModel.reactions]
    # create a dictionary with the gpr rules associated to each reaction to get a unique list of gprs
    gprDict = getGPRDict(cobraModel)
    gDict = {}
    
    # Analyze the GPRs to find the perturbations that stop the reactions
    for key, value in tqdm(gprDict.items()):
            calculate_GPR_MCS(gDict, value['gpr'], value['reactions'], dictManager)


    # Filter out the interventions that are larger than the maxKOLength
    filteredGDict = gDict.copy()
    for key in gDict.keys():
        if len(key) > maxKOLength:
            filteredGDict.pop(key)
            
    newGDict = transformReactionIdsIntoIndexes(filteredGDict, modelReactions)
    
            
    simpG = simplifyGMatrix(newGDict, dictManager)

    [relationships, numberNewGenesByKO] = relatedRows(simpG, dictManager)
    
    gMatrix = createSparseMatrix(simpG, modelReactions)

    # scipy.sparse.save_npz(cobraModelName + '_gMatrix.npz', gMatrix)
    gObject = {
        "gMatrix": gMatrix,
        "gDict": simpG,
        "relationships": relationships,
        "numberNewGenesByKO": numberNewGenesByKO,
    }
    return gObject

def calculate_GPR_MCS(gDict, gpr, reactions, addToDict):
    # get the genes from the gpr rule
    stringGPR = gpr.to_string()
    genes = [gene for gene in gpr.genes]
    if len(genes) == 1:
        key=frozenset(genes)
        [addToDict(gDict, key, reaction) for reaction in reactions]
    # Only and
    elif bool(re.search(' and ', stringGPR) and not re.search(' or ', stringGPR)):
        # Only one gene needs to be KO to stop the reaction
        for gene in genes:
            key=frozenset([gene])
            [addToDict(gDict, key, reaction) for reaction in reactions]

    # Only or
    elif bool(re.search(' or ', stringGPR) and not re.search(' and ', stringGPR)):
        # All genes need to be KO to stop the reaction
        key=frozenset(genes)
        [addToDict(gDict, key, reaction) for reaction in reactions]  
    
    # And and or
    else:
        treeGenes, model = parseGPRToModel(stringGPR, 'gpr', None)
        solutionDict = calculateMCS(model, MAX_LENGTH_MCS=5, MAX_MCS=1e6, rxnSubset=list(treeGenes))
        for _, data in solutionDict.items():
            key = frozenset(data["solution"])
            [addToDict(gDict, key, reaction) for reaction in reactions]
                   
def getGPRDict(model):
    gpr_dict = {}
    for reaction in model.reactions:
        gpr = reaction.gpr.to_string()
        if gpr in gpr_dict:
            gpr_dict[gpr]['reactions'].append(reaction.id)
        elif gpr == '':
            pass
        else:
            gpr_dict[gpr] = {'reactions': [reaction.id], 'genes': [gene.id for gene in reaction.genes], 'gpr': reaction.gpr}
    return gpr_dict

def simplifyGMatrix(gMatrixDict: defaultdict, __addToDict) -> defaultdict:
    """**Simplify the G matrix**

    Analyze the G matrix to find the rows that are related to each other.
    If a common gene is more perturbations, then the perturbations are related to each other.
    And should be an individual perturbation.

    :param gMatrixDict: The dictionary with the perturbations.

    :return: The simplified G matrix.
    """
    
    # Empty dictionary to store the simplified G matrix
    simpG = gMatrixDict.copy()
    
    # Get the interventions compossed of one gene
    singleGeneInterventions = [gene for gene in gMatrixDict.keys() if len(gene) == 1]
    
    # Get the interventions compossed of more than one gene
    multipleGeneInterventions = [gene for gene in gMatrixDict.keys() if len(gene) > 1]
    
    # Intersect each intervention with the other interventions to find common elements
    for interventionA in multipleGeneInterventions:
        for interventionB in multipleGeneInterventions:
            if interventionA != interventionB:
                intersection = interventionA.intersection(interventionB)
                if intersection:
                    # If the intersection is not empty, then the interventions are related
                    # Check if the intersection is already a single gene intervention
                    if len(intersection) == 1:
                       if intersection not in singleGeneInterventions and intersection not in gMatrixDict.keys():
                           # Intersection must not be a key of the dictionary
                            if intersection not in simpG.keys():
                                 # Add the intersection as a key of the dictionary
                                __addToDict(simpG, intersection, [])
                                simpG[intersection] = []
                    else:
                        for element in intersection:
                            if element not in singleGeneInterventions and element not in gMatrixDict.keys():
                                # Intersection must not be a key of the dictionary
                                if intersection not in simpG.keys():
                                     # Add the intersection as a key of the dictionary
                                    __addToDict(simpG, intersection, [])
                                    simpG[intersection] = []
                                

    return simpG

def createSparseMatrix(gMatrixDict, modelReactions):
    """**Build the G matrix**

    Create sparse matrix from the dictionary with the perturbations.

    :param gMatrixDict: The dictionary with the perturbations.
    :param modelReactions: The reactions of the model.

    :return: The G matrix.
    """
    gMatrix = scipy.sparse.lil_matrix((len(gMatrixDict), len(modelReactions)))
    for i, items in enumerate(gMatrixDict.items()):
        gMatrix[i, list(items[1])] = 1
    gMatrix = gMatrix.tocsc()
    return gMatrix

def relatedRows(gMatrixDict: defaultdict, __addToDict) -> List:
    """**Find the rows that are related to each other**

    If a row is a subset of another row,
    then that means that when the perturbation of the superset is activated, the subset is also affected.
    i. e. gene A and gene B inactivates reaction1, and gene A inactivates reaction2. Knocking out gene A and gene B will affect both reactions.

    :param gMatrixDict: The dictionary with the perturbations.

    :return: A list with the rows that are related to each other.
    """
    relationships = defaultdict(list)
    numberNewGenesByKO = defaultdict(int)
    sortedDict = sorted(gMatrixDict.items(), key=lambda x: len(x[0]))
    for key, value in gMatrixDict.items():
        numberNewGenesByKO[key] = len(key)
        if len(key) < 2:
            continue
        keysConstained = []
        for key2, value2 in sortedDict:
            if len(key2) >= len(key):
                break
            if key2.issubset(key):
                keysConstained.append(key2)
                value.extend(value2)
                __addToDict(relationships, key, key2)
        if len(keysConstained) > 0:
            key2 = frozenset.union(*keysConstained)
            numberNewGenesByKO[key] = len(key) - len(key2)
        gMatrixDict[key] = list(set(value))
    return [relationships, numberNewGenesByKO]

def calculateMCSForGPRs(
    briefing: dict, isoformSeparator: str, verbose: int = 0, **kwargs
):
    """**Calculate the Minimal Cut Sets for GPRs with both AND and OR relationships**

    Given a complex GPR, we parse de GPR to a cobra model and calculate the MCSs for the model.

    :param briefing: The dictionary with the information of the GPRs.
    :param verbose: The verbosity level.

    :return: A List with the List of GPRs and a List with the MCSs.
    """
    gpr = []
    mcs = []
    parmMCS = {
        "MAX_MCS": 10000,
        "targetB": 1e-3,
        "timelimit": 60,
        "forceLength": True,
        "numWorkers": 0,
        "verbose": verbose,
        "useIndicators": True,
        "solver": "gurobi",
    }
    for key, value in kwargs.items():
        if key in parmMCS:
            parmMCS[key] = value
        else:
            warnings.warn("Parameter " + key + " is not defined in calculateMCS.")

    for i, stringGPR in enumerate(briefing["gprsANDOR"]):
        treeGenes, model = parseGPRToModel(stringGPR, i, isoformSeparator)
        if verbose > 0:
            for i in model.reactions:
                print(i)

        solutionDict = calculateMCS(
            model, MAX_LENGTH_MCS=len(treeGenes), rxnSubset=list(treeGenes), **parmMCS
        )
        for key, data in solutionDict.items():
            gpr.append(stringGPR)
            mcs.append(data["solution"])
    return gpr, mcs

@beartype
def parseGPRToModel(stringGPR: str, reactionName, isoformSeparator):
    """**Parse a GPR string to a cobra model**


    The model is used to calculate the MCS to calculate the combinations of genes that are required for the reaction.

    :param stringGPR: The GPR string to parse. (e.g. "gen1 and gen2 or gen3")
    :param reactionName: The name of the reaction to parse. Just to identify the model.

    :return:
        - **treeGenes** - A list of genes that participate in the GPR.
        - **model** - A cobra model with the GPR parsed to a metabolic model.
    """
    if isoformSeparator is not None:
        # Remove the isoform separator and the following string and numbers
        stringGPR = re.sub(f"([{isoformSeparator}][0-9]+)", "", stringGPR)
    [gprTree, treeGenes] = cobra.core.gene.parse_gpr(stringGPR)
    parentMetabolite = cobra.core.Metabolite(f"r{0}", name=f"r{0}")
    model = cobra.core.Model(f"reaction_{reactionName}")
    objectiveReaction = cobra.core.Reaction("objective")
    objectiveReaction.name = "objective"
    objectiveReaction.lower_bound = 0
    objectiveReaction.upper_bound = 1000
    objectiveReaction.add_metabolites({parentMetabolite: -1.0})
    model.add_reactions([objectiveReaction])
    model.objective = "objective"
    level = 0
    model = buildReactionFromBranch(expr=gprTree.body, model=model, level=level)

    return [treeGenes, model]

@beartype
def buildReactionFromBranch(
    expr: BoolOp,
    model: cobra.core.Model,
    level: int = 0,
    parentMetabolite: cobra.core.Metabolite = None,
) -> cobra.core.Model:
    """**Recursive function to generate metabolic model from GPR.**

    Each GPR can be expressed as a metabolic model. To that end, a boolean tree can be build, where a branch can be either an AND, an OR or a gene.
    We move through the tree, and recursively call this function on each node. Until all genes are explored (end nodes).
    As we move through the tree, we add reactions to the model. The reactions are named with a unique identifier, to avoid name conflicts.
    An example building a metabolic model from a GPR is shown below.

    .. figure:: ./figs/binaryTreeExample.png
        :height: 350px
        :width: 450 px
        :align: center

        Figure 1. Example of a binary tree.

    In the example above, the GPR is (g2 and (g5 or g6)) or g7. The function will be called recursively on each node.

    :param expr: A BoolOp object from the ast module to be explored. i.e. a tree.  (ast.BoolOp)
    :param model: A cobra model to build add the tree.
    :param level: The level of the tree, i.e. the depth of the current node.  (int)
    :param parentMetabolite: The metabolite that is the parent of the current node.  (cobra.core.Metabolite)

    :returns:
            - **model**: A cobra model with the tree added as reactions, genes have become metabolites, and the objective function set to the root node, the reaction of interest.
    """
    if isinstance(expr, BoolOp):
        branch = expr.op
        if isinstance(branch, Or):
            if parentMetabolite is None:
                parentMetabolite = cobra.core.Metabolite(
                    f"r{level}", name=f"reaction_{level}"
                )
            reaction = cobra.core.Reaction(str(uuid.uuid4()))
            reaction.name = f"OR_{level}_{parentMetabolite.id}"
            reaction.lower_bound = 0
            reaction.upper_bound = 1000
            reaction.add_metabolites({parentMetabolite: 1.0})

            for i in expr.values:
                if isinstance(i, Name):
                    level = level + 1
                    reaction = cobra.core.Reaction(str(uuid.uuid4()))
                    reaction.name = f"OR_{level}_{parentMetabolite.id}"
                    reaction.lower_bound = 0
                    reaction.upper_bound = 1000
                    reaction.add_metabolites({parentMetabolite: 1.0})
                    childMetabolite = cobra.core.Metabolite(i.id, name=i.id)
                    reaction.add_metabolites({childMetabolite: -1.0})
                    model.add_reactions([reaction])
                    if i.id in model.reactions:
                        continue
                    reaction = cobra.core.Reaction(i.id)
                    reaction.name = i.id
                    reaction.lower_bound = 0
                    reaction.upper_bound = 1000
                    reaction.add_metabolites({childMetabolite: 1.0})
                    model.add_reactions([reaction])
                else:
                    level = level + 1
                    reaction = cobra.core.Reaction(str(uuid.uuid4()))
                    reaction.name = f"OR_{level}_{parentMetabolite.id}"
                    reaction.lower_bound = 0
                    reaction.upper_bound = 1000
                    reaction.add_metabolites({parentMetabolite: 1.0})
                    childMetabolite = cobra.core.Metabolite(
                        f"branch_OR_met{level}_{parentMetabolite.id}",
                        name=f"reaction_{level}",
                    )
                    reaction.add_metabolites({childMetabolite: -1.0})
                    model.add_reactions([reaction])
                    buildReactionFromBranch(i, model, level, childMetabolite)
        elif isinstance(branch, And):
            if parentMetabolite is None:
                parentMetabolite = cobra.core.Metabolite(
                    f"r{level}", name=f"reaction_{level}"
                )
            reaction = cobra.core.Reaction(str(uuid.uuid4()))
            reaction.name = f"AND_{level}_{parentMetabolite.id}"
            reaction.lower_bound = 0
            reaction.upper_bound = 1000
            reaction.add_metabolites({parentMetabolite: 1.0})
            for i in expr.values:
                if isinstance(i, Name):
                    childMetabolite = cobra.core.Metabolite(i.id, name=i.id)
                    reaction.add_metabolites({childMetabolite: -1.0})
                    if i.id in model.reactions:
                        continue
                    newReaction = cobra.core.Reaction(i.id)
                    newReaction.name = i.id
                    newReaction.lower_bound = 0
                    newReaction.upper_bound = 1000
                    newReaction.add_metabolites({childMetabolite: 1.0})
                    model.add_reactions([newReaction])
                else:
                    level = level + 1
                    childMetabolite = cobra.core.Metabolite(
                        f"branch_AND_met{level}_{parentMetabolite.id}",
                        name=f"reaction_{level}",
                    )
                    reaction.add_metabolites({childMetabolite: -1.0})
                    buildReactionFromBranch(i, model, level, childMetabolite)
            model.add_reactions([reaction])

    return model


def mergeIsforms(isoformSeparator):
    if isoformSeparator is None:

        def __addToDict(dict: defaultdict, key: List, value) -> None:
            """**Add to dictionary**

            Internal function to add a value to a dictionary with a key. If the key does not exist, it is created,
            and the value is added to the list of values.
            If the key exists, the value is appended to the list of values.

            Dictionary mut be a `defaultdict <https://docs.python.org/3/library/collections.html#collections.defaultdict>`_.

            :param dict: Dictionary to add gene reaction interaction. (defaultdict)
            :param key: List of genes that render the reaction inactivate. (list), if is only one gene, it must be a list with one element.
            :param value: Any value can be stored in the dictionary, we store the corresponding index of the reactions perturbed. (any)

            """
            #key = frozenset(key)
            if key in dict:
                dict[key].append(value)
            else:
                dict[key] = [value]

    elif isinstance(isoformSeparator, str):

        def __addToDict(dict: defaultdict, key: List, value) -> None:
            key = [name.split(isoformSeparator)[0] for name in key]
            key = frozenset(key)
            if key in dict:
               dict[key].append(value)
            else:
                dict[key] = [value]

    else:
        raise TypeError("isoformSeparator must be a string or None")

    return __addToDict


@beartype
def buildRxnGeneMat(GPRs: List, reactions: List, genes: List):
    """**Build a reaction-gene matrix from a list of GPRs.**

    The binary matrix is build from each reaction GPR,
    each row corresponds to a reaction and each column corresponds to a gene.
    The matrix is sparse and is stored in the CSR format.

    :param GPRs: List of GPRs from a cobra model.
    :param reactions: List of reactions from a cobra model.
    :param genes: List of genes from a cobra model.

    :returns:
        - **rxnGeneMatrix:** A sparse matrix with the number of reactions in rows and the number of genes in columns. (scipy.sparse.csr_matrix)
        - **summary:** A dictionary with the number of reactions classified by the number of genes in their GPR. (dict)
    """
    rxnGeneMatrix = scipy.sparse.lil_matrix((len(reactions), len(genes)))

    # Build reaction-gene matrix
    for position, reaction in enumerate(reactions):
        GPR = reaction.gpr.genes
        for gene in GPR:
            rxnGeneMatrix[position, genes.index(gene)] = 1

    rxnGeneMatrix = rxnGeneMatrix.tocsr()

    rxnZeroGenes = np.where(rxnGeneMatrix.sum(axis=1) == 0)[0]
    rxnOneGene = np.where(rxnGeneMatrix.sum(axis=1) == 1)[0]

    # indices of reactions not in rxnZeroGenes and rxnOneGene
    indexRxnMoreThanOneGene = np.setdiff1d(
        np.arange(len(reactions)), np.concatenate((rxnZeroGenes, rxnOneGene))
    )

    # Keep only GPRs that include two or more genes
    GPRStrings = []
    for gpr in GPRs:
        if len(gpr.genes) >= 2:
            GPRStrings.append(gpr.to_string())

    # identify reactions with "or" and "and" relationships between genes
    rxnOr = [bool(re.search(" or ", i)) for i in GPRStrings]
    rxnAnd = [bool(re.search(" and ", i)) for i in GPRStrings]

    rxnOnlyOr = np.logical_and(rxnOr, np.logical_not(rxnAnd))
    rxnOnlyAnd = np.logical_and(rxnAnd, np.logical_not(rxnOr))
    rxnOrAnd = np.logical_and(rxnOr, rxnAnd)

    dictGPRstring = {
        position: reactions[position].gpr.to_string()
        for position in indexRxnMoreThanOneGene
    }

    gprOnlyOr = np.where(rxnOnlyOr)[0]
    indexRxnOnlyOr = []
    for gpr in gprOnlyOr:
        uniqueGPR = GPRStrings[gpr]
        indexRxnOnlyOr.extend(
            [key for key, value in dictGPRstring.items() if value == uniqueGPR]
        )

    gprRxnOnlyAnd = np.where(rxnOnlyAnd)[0]
    indexRxnOnlyAnd = []
    for gpr in gprRxnOnlyAnd:
        uniqueGPR = GPRStrings[gpr]
        indexRxnOnlyAnd.extend(
            [key for key, value in dictGPRstring.items() if value == uniqueGPR]
        )

    gprRxnOrAnd = np.where(rxnOrAnd)[0]
    indexRxnOrAnd = []
    for gpr in gprRxnOrAnd:
        uniqueGPR = GPRStrings[gpr]
        indexRxnOrAnd.extend(
            [key for key, value in dictGPRstring.items() if value == uniqueGPR]
        )

    # count number of reactions in each category
    nRxnZeroGenes = len(rxnZeroGenes)
    nRxnOneGene = len(rxnOneGene)
    nRxnOnlyOr = len(indexRxnOnlyOr)
    nRxnOnlyAnd = len(indexRxnOnlyAnd)
    nRxnOrAnd = len(indexRxnOrAnd)

    nRxnTotal = nRxnZeroGenes + nRxnOneGene + nRxnOnlyOr + nRxnOnlyAnd + nRxnOrAnd

    # Keep only GPRs that will be used in the analysis of minimal cut sets
    UniqueGPRStringsANDOR = [GPRStrings[i] for i in gprRxnOrAnd]

    # summarize results as a dictionary
    summary = {
        "nRxnZeroGenes": nRxnZeroGenes,
        "indexRxnZeroGenes": rxnZeroGenes,
        "nRxnOneGene": nRxnOneGene,
        "indexRxnOneGene": rxnOneGene,
        "nRxnOnlyOr": nRxnOnlyOr,
        "indexRxnOnlyOr": indexRxnOnlyOr,
        "nRxnOnlyAnd": nRxnOnlyAnd,
        "indexRxnOnlyAnd": indexRxnOnlyAnd,
        "nRxnOrAnd": nRxnOrAnd,
        "gprsANDOR": UniqueGPRStringsANDOR,
        "gprDict": dictGPRstring,
        "nRxnTotal": nRxnTotal,
    }

    return [rxnGeneMatrix, summary]


def initiateLogger(name, log_file, level=logging.INFO):
    handler = logging.FileHandler(log_file)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def transformReactionIdsIntoIndexes(gDict, reactionList):
    newGDict = defaultdict(list)
    for key, value in gDict.items():
        reactionIndexes = []
        for reaction in value:
            reactionIndexes.append(reactionList.index(reaction))
        newGDict[key] = reactionIndexes
    return newGDict



        