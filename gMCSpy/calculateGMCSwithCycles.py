from .ProblemDefinitions import cobra
from .ProblemInterpreters import np

from .OptimizationProblem import OptimizationProblem

from .ProblemInterpreters import beartype
from .ProblemDefinitions import List
from .ProblemDefinitions import scipy
from .ProblemDefinitions import buildDictionaryDossageNetwork, buildDictionaryMCSGeneProblem

from .Validations import checkGMCSParallel

from .calculateMCS import calculateMCS

from .calculateGMCS import buildReactionFromBranch, createSparseMatrix, Name, parse, re, calculateMCS, parseGPRToModel, simplifyGMatrix, relatedRows, mergeIsforms, initiateLogger

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

from bidict import bidict

from bonesis.reprogramming import marker_reprogramming
import mpbn

import sys
import os
from contextlib import contextmanager

@contextmanager
def suppress_output():
    # Save current stdout and stderr
    stdout = sys.stdout
    stderr = sys.stderr
    # Redirect stdout and stderr to devnull (a special file that discards all data written to it)
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')
    try:
        yield
    finally:
        # Restore stdout and stderr
        sys.stdout = stdout
        sys.stderr = stderr

def calculateGMIS(cobraModel, regulatory_dataframe, num_layers=2, **kwargs):
    '''
    '''
    # Name of the analysis
    name = (
        "gMCSpy_SyntheticDosage_"
        + cobraModel.id
        + "_"
        + datetime.now().strftime("%H-%M-%S--%d-%m-%Y")
    )
    
    # Prefered solver is gurobi due to performance experience during benchmarking
    solver = kwargs.get("solver", "gurobi")
    maxKOLength = kwargs.get('maxKOLength', 3)
    timeLimit = kwargs.get("timeLimit", 1e4)
    numWorkers = kwargs.get("numWorkers", 0)
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
    if kwargs.get("forceLength") is None:
        if solver == "gurobi":
            forceLength = False
        elif solver == "cplex":
            forceLength = True
        elif solver == "scip":
            forceLength = True   
    else:
        forceLength = kwargs.get("forceLength", False)
        
        
    # Save all the parameters in a log file    
    path = "logs/" + name + "/"
    createLog(path)
    configName = ("configuration")
    logging = initiateLogger(configName + ".log", path + configName + ".log")
    # Save the parameters in the log file
    logging.info(f"solver: {solver}, maxKOLength: {maxKOLength}, timeLimit: {timeLimit}, numWorkers: {numWorkers}, maxNumberGMCS: {maxNumberGMCS}, num_layers: {num_layers}, verbose: {verbose}, forceLength: {forceLength}, earlyStop: {earlyStop}, removeRelated: {removeRelated}, isoformSeparator: {isoformSeparator}, targetKOs: {targetKOs}, geneSubset: {geneSubset}, isNutrient: {isNutrient}, exchangeReactions: {exchangeReactions}, saveGMatrix: {saveGMatrix}, checkSolutions: {checkSolutions}, saveSolutions: {saveSolutions}, forceLength: {forceLength}")
    logging.info(f"Model: {cobraModel.id}")
    # close logger
    handler = logging.handlers[0]
    handler.close()
    logging.removeHandler(handler)
    
    # Pop the parameters that are not used downstream
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
    path = "logs/" + name + "/process/"
    createLog(path)
    path = "logs/" + name
    name = ("process")
    handler = name + ".log"
    logging = initiateLogger(name + ".log", path + "/process/" + name + ".log")
    allTime = time.time()

    
       
    if isNutrient:
        cobraModel = prepareModelNutrientGeneMCS(cobraModel, exchangeReactions)

    # Prepare the regulatory network
    regulatory_dict = regulatory_dataframe.to_dict('list')
    # Check if the regulatory network has valid sources and targets
    regulatory_dict, referenceDict = checkRegulatoryNetwork(regulatory_dict)
    # Check if the regulatory network has the correct format, should have source_ENSEMBL, target_ENSEMBL, and interaction
    expectedColumns = ['source_ENSEMBL', 'target_ENSEMBL', 'interaction']
    for column in expectedColumns:
        if column not in regulatory_dict:
            raise ValueError(f"Expected column {column} not found in the regulatory network")
    
    startTime = time.time()
    gObject = calculateRegNetGMatrix(cobraModel, regulatory_dict, num_layers, solver, maxKOLength=maxKOLength)
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
        "PoolGap": 0.01,
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
            ko_log_decode = []
            for gene in ko['solution']:
                  if gene.startswith('M_I_'):
                      ko_log_decode.append(decode_string(gene, referenceDict))
                  else:
                      ko_log_decode.append(gene)                
            loggingSolutions.info(f'{ko_log_decode}:[{reactions}]')

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

def interpretExpandedRules(key, expandedRules, genesGPR):
    simpleGPRs = {}
    for source, target, interaction in zip(expandedRules['sources'], expandedRules['targets'], expandedRules['interaction']):
        if interaction == 1:
            if target in simpleGPRs:
                simpleGPRs[target]['rule'] = f'{simpleGPRs[target]["rule"]} | {source}'
                simpleGPRs[target]['genes'].append(source)
            else:
                simpleGPRs[target] = {'rule': f'{target} = {source}', 'genes': [source]}
        elif interaction == -1:
            if target in simpleGPRs:
                simpleGPRs[target]['rule'] = f'{simpleGPRs[target]["rule"]} | !{source}'
                simpleGPRs[target]['genes'].append(source)
            else:
                simpleGPRs[target] = {'rule': f'{target} = !{source}', 'genes': [source]}
                
    sourceSet= set(expandedRules['sources']).union(set(expandedRules['targets']))
    expandedRules = {}
    expandedRules['gpr'] = {'rule': f'gpr = {key}'}
    
    
    for gene in genesGPR:
        if gene not in sourceSet:
            simpleGPRs[gene] = {'rule': 'noRegulatoryInteraction', 'genes': [gene]}
    
    
    for aReaction, boolRule in simpleGPRs.items():
        KO_gene = f'{aReaction}_KO'
        KI_gene = f'{aReaction}_KI'
        if boolRule['rule'] == 'noRegulatoryInteraction':
            expandedRules[aReaction] = {'rule': f'{aReaction} = {KO_gene} | {KI_gene}'}
        elif len(boolRule['genes']) > 1:
            rule = boolRule['rule'].split(' = ')
            expandedRules[aReaction] = {'rule': f"{rule[0]} = (({rule[1]}) | {KI_gene}) & {KO_gene}"}
            for gene in boolRule['genes']:
                if gene in simpleGPRs.keys():
                    continue
                KO_gene = f'{gene}_KO'
                KI_gene = f'{gene}_KI'
                expandedRules[gene] = {'rule': f'{gene} = {KO_gene} | {KI_gene}'}
        else:
            new_aReaction = aReaction.split("!")[-1]
            KO_gene = f'{new_aReaction}_KO'
            KI_gene = f'{new_aReaction}_KI'
            rule = boolRule['rule'].split(' = ')
            expandedRules[new_aReaction] = {'rule': f"{rule[0]} = (({rule[1]}) | {KI_gene}) & {KO_gene}"}
            for gene in boolRule['genes']:
                if gene in simpleGPRs.keys():
                    continue
                KO_gene = f'{gene}_KO'
                KI_gene = f'{gene}_KI'
                expandedRules[gene] = {'rule': f'{gene} = {KO_gene} | {KI_gene}'}
    return expandedRules

def bonesisCalculateCutSets(boolForm, maxKOLength=5):
    f = mpbn.MPBooleanNetwork(boolForm)
    res = marker_reprogramming(f, {"gpr": 0}, maxKOLength, ensure_exists=True)
    #tabulate(list(f.attractors()))
    #print(list(res))
    #print('done')
    return list(res)

def boolTransformer(expandedRules, value):
    newDict = {}
    # transform the grp into a boolean form
    trans_gpr = str(value['gpr'])
    trans_gpr = trans_gpr.replace(' and ', ' & ').replace(' or ', ' | ')
    newDict['gpr'] = trans_gpr

    for source, target, interaction in zip(expandedRules['sources'], expandedRules['targets'], expandedRules['interaction']):
        if interaction == -1:
            source = f'!{source}'
        if target not in newDict.keys():
            newDict[target] = source
        else:
            newDict[target] = f'{newDict[target]} | {source}'
    
    for gene in expandedRules['genes']:
        if gene not in newDict.keys():
            newDict[gene] = gene
            
    return newDict

def calculateRegNetGMatrix(model, regulatory_dict, num_layers=2, solver='gurobi', **kwargs):    
    isoformSeparator = kwargs.get('isoformSeparator', None)
    maxKOLength = kwargs.get('maxKOLength', 3)
    dictManager = mergeIsforms(isoformSeparator)
        
    modelReactions = [reaction.id for reaction in model.reactions]
    # create a dictionary with the gpr rules associated to each reaction to get a unique list of gprs
    gprDict = getGPRDict(model)
    gDict = {}

    
    # go through each gpr rule and expand it to the number of layers, if the rule is not expanded use traditional MCS
    geneSet = set()
    for key, value in tqdm(gprDict.items()):
        for iterat in value['gpr'].genes:
            geneSet.add(iterat)
        expandedRulesWithCycles = stackLayersWithCycles(value['genes'], regulatory_dict, num_layers, solver)
        for iterat in expandedRulesWithCycles['sources']:
            geneSet.add(iterat)
        for iterat in expandedRulesWithCycles['targets']:
            geneSet.add(iterat)
        if len(expandedRulesWithCycles['sources']) == 0:
            calculateTraditionalMCS(gDict, value['gpr'], value['reactions'])
        else:
            boolForm = boolTransformer(expandedRulesWithCycles, value)
            with suppress_output():
                gmcs = bonesisCalculateCutSets(boolForm, maxKOLength = 5)
            if len(gmcs) == 1:
                # This means the only solution found is blocking the original gpr
                continue
            for solution in gmcs:
                if list(solution.keys())[0] == 'gpr':
                    continue
                result_set = []
                for key, val in solution.items():
                    lastName = None
                    if val == 0:
                        lastName = 'KO'
                    elif val == 1:
                        lastName = 'KI'
                    sol = key + '_' + lastName
                    result_set.append(sol)
                #print(f'The gmcs is {solution}')
                #print(f'The result_set is {result_set}')            
                result_set = frozenset(result_set)
                addToDict(gDict, result_set, value['reactions'])

    nameGenesFile = f'Genes_{datetime.now().strftime("%H-%M-%S--%d-%m-%Y")}.txt'
    filename = 'C:/Users/cjrodriguezf/Documents/PhD/Loic/evaluation/Data/Gene_lists_prev/' + nameGenesFile
    # Write set to file, create file if it doesn't exist
    with open(filename, 'w') as file:
        file.write('Genes\n')
        for gene in geneSet:
            file.write(gene + '\n')

    # Make sure that each entry in the gDict is unique list of reactions
    for key, value in gDict.items():
        gDict[key] = list(set(value))
    
    # Filter out the interventions that are larger than the maxKOLength
    filteredGDict = gDict.copy()
    for key in gDict.keys():
        if len(key) > maxKOLength:
            filteredGDict.pop(key)
    
    newGDict = transformReactionIdsIntoIndexes(filteredGDict, modelReactions)
    
    simpG = simplifyGMatrix(newGDict, dictManager)
    
    [relationships, numberNewGenesByKO] = relatedRows(simpG, dictManager)
    
    gMatrix = createSparseMatrix(simpG, modelReactions)
    
    gObject = {
        "gMatrix": gMatrix,
        "gDict": simpG,
        "relationships": relationships,
        "numberNewGenesByKO": numberNewGenesByKO,
    }
    
    
    
    return gObject

def calculateGMCSofGPR_no_cycles(key, expandedRules, value, solver='gurobi'):
    expandedRules = interpretExpandedRules(key, expandedRules, value['genes'])
    duplicatedNetwork = []
    for _, rule in expandedRules.items():
        duplicatedNetwork.extend(parseRule(rule['rule']))
    
    solutionDict = calculateRegNetMCS(solver, duplicatedNetwork)
    for _, data in solutionDict.items():
        geneSet = [gene.split('zp_')[1] for gene in data["solution"]]
        result_set = []
        for solution in geneSet:
            if 'KO' in solution:
                result_set.append(solution.split('_on')[0])
            else:
                result_set.append(solution.split('_off')[0])
        result_set = frozenset(result_set)
    return result_set
                
def calculateTraditionalMCS(gDict, gpr, reactions):
    # get the genes from the gpr rule
    stringGPR = gpr.to_string()
    genes = [gene for gene in gpr.genes]
    if len(genes) == 1:
        # Rename the genes as they are KOs
        genes = [f'{gene}_KO' for gene in genes]
        key=frozenset(genes)
        addToDict(gDict, key, reactions)
    # Only and
    elif bool(re.search(' and ', stringGPR) and not re.search(' or ', stringGPR)):
        # Rename the genes as they are KOs
        genes = [f'{gene}_KO' for gene in genes]
        # Only one gene needs to be KO to stop the reaction
        for gene in genes:
            key=frozenset([gene])
            addToDict(gDict, key, reactions)

    # Only or
    elif bool(re.search(' or ', stringGPR) and not re.search(' and ', stringGPR)):
        # Rename the genes as they are KOs
        genes = [f'{gene}_KO' for gene in genes]
        # All genes need to be KO to stop the reaction
        key=frozenset(genes)
        addToDict(gDict, key, reactions)   # And and or
    else:
        treeGenes, model = parseGPRToModel(stringGPR, 'gpr', None)
        solutionDict = calculateMCS(model, MAX_LENGTH_MCS=5, MAX_MCS=1e6, rxnSubset=list(treeGenes))
        for _, data in solutionDict.items():
            # Rename the genes as they are KOs
            genes = [f'{gene}_KO' for gene in data["solution"]]
            key = frozenset(genes)
            addToDict(gDict, key, reactions)
        
def addToDict(gDict, key, reactions):
    if key in gDict:
        gDict[key] = gDict[key] + (reactions)
    else:
        gDict[key] = reactions
            
def calculateRegNetMCS(solver, duplicatedNetwork, maxKOLength=5):
    model = duplicatedNetworkToModel(duplicatedNetwork)  
    rxnToForce = [reaction.id for reaction in model.reactions if 'redflux' in reaction.id]
    compressedNodes = set()

    for reaction in rxnToForce:
        if '_KI_off' in reaction:
            compressedNodes.add(reaction.split('_KI_off')[0])
        elif '_KO_on' in reaction:
            compressedNodes.add(reaction.split('_KO_on')[0])
    problem = buildDictionaryDossageNetwork(model, forceLength=True, rxnToForce=rxnToForce, compressedNodes=compressedNodes)
        
    timelimit = 500
    numWorkers = 0
        
    gurobi_default_params = {
        "MIPFocus": 3,
        "TimeLimit": max(1, timelimit),
        "Threads": numWorkers,
        "PoolSolutions": 100,
        "PoolGap": 0.01,
        "PoolSearchMode": 1,
        "Presolve": 1,
        }

    cplex_default_params = {
        "mip.tolerances.integrality": 1e-5,
        "mip.strategy.rinsheur": 50,
        "emphasis.mip": 1,
        "timelimit": max(10, timelimit),
        "threads": numWorkers,
        "mip.limits.populate": 20,
        "mip.pool.relgap": 0.001,
        }
    if solver == "gurobi":
        params = gurobi_default_params
    elif solver == "cplex":
        params = cplex_default_params
    else:
        params = {}      
            
    problem.setSolver(solver)
        
        # Interpret the problem, depending on the solver each problem has to be interpreted differently into the solver's interface
    optimization = problem.interpretProblem(params) 
        
        # Import the common problem methods
    problemMethods = MinimalCutSetProblem()

        # Get the objective variables, they are the variables that we are interested in
    obj = problem.getObjective()

        # Set the bounds for the solver, e.g. maximum number of solutions, maximum numbers in a solution, or if we want to force the step by step calculation etc.
    bounds = {
            "MAX_SOLUTION_LENGTH": 7,
            "MAX_SOLUTIONS": 10000,
            "forceLength": True,
        }

    varNames = []
    for i in obj["variables"]:
        varNames.append(problem.getVariable(i)["name"])
        
    obj = obj["variables"]
    solutionVariables = {"indices": obj, "names": varNames}

        # Solve the problem iteratively
    solutionDict = problemMethods.iterateOverSolver(
            problem=problem,
            problemInterpreted=optimization,
            solver=solver,
            solutionVariables=solutionVariables,
            bounds=bounds,
        )
    
    return solutionDict

def duplicatedNetworkToModel(duplicatedNetwork):
    model = cobra.Model('duplicated_model')
    react_list = []
    blueFluxes = []
    for gpr in duplicatedNetwork:
        model, react = reactionsToModel(gpr, model)    
        if react is not None:
            react_list.append(react)
        [product, reactant] = gpr.split('=')
        if '_off' in product:
            gene = product.split('_off')[0]
            blue_gpr = f'{gene}_off = {gene}_KO_off & {gene}_KI_off'
            if blue_gpr == gpr:
                blueFluxes.append(f'{gene}_KO_off')
            
        
    
    reactionListRemove = []
    for reaction in model.reactions:
        if len(reaction.reactants) == 0:
            reactionListRemove.append(reaction.id)
            
    model.remove_reactions(reactionListRemove, remove_orphans=True)            
    
    model.repair()
    
    for metabolite in model.metabolites:
        if '_K' in metabolite.id:
            name = f'{metabolite.id}_redflux'
            redFlux = cobra.core.Reaction(name)
            redFlux.name = name
            redFlux.lower_bound = 0
            redFlux.upper_bound = 10000
            redFlux.add_metabolites({metabolite: 1.0})
            model.add_reactions([redFlux])      
                
    for blue in blueFluxes:  
        metabolite = model.metabolites.get_by_id(blue) 
        name = f'{blue}_blueflux'
        blueflux = cobra.core.Reaction(name)
        blueflux.name = name
        blueflux.lower_bound = 0
        blueflux.upper_bound = 100
        blueflux.add_metabolites({metabolite: 1.0})
        model.add_reactions([blueflux])
        
    gpr_reactions = addOutputReactions(model)
    model.add_reactions(gpr_reactions)
    model.objective = model.reactions.get_by_id("gpr_on")
    return model
    
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
    
def createRegulatoryModel(sources, targets):
    model = cobra.Model('regulatory_model')
    reactions = []
    for reactant, product in zip(sources, targets):
        reaction = cobra.Reaction(f'{reactant}__{product}')   
        reaction.name = f'{reactant}__{product}'
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
        reactant_metabolite = cobra.Metabolite(id=reactant, name=reactant)
        product_metabolite = cobra.Metabolite(id=product, name=product)
        reaction.add_metabolites({reactant_metabolite: -1, product_metabolite: 1})
        reactions.append(reaction)
    model.add_reactions(reactions)
    return model

def solveForCycles(model, solver):
    formulation = OptimizationProblem()
    S = cobra.util.array.create_stoichiometric_matrix(model)
    [rows, cols] = S.shape
    # Append a row of ones at the end of the matrix
    S = np.vstack([S, np.ones(S.shape[1])])
    
    for i in range(cols):
        formulation.addVariable(index=i, prefix=f'v_', vtype='continuous', lb=0, ub='infinity')

    for i in range(rows):
        formulation.addConstraint(name=f'row_{i}', variables=formulation.getVariables(), coefficients=S[i], sense='E', rhs=0)
        
    formulation.addConstraint(name='row_ones', variables=formulation.getVariables(), coefficients=S[-1], sense='G', rhs=1)

    variables = list(formulation.getVariables().keys())
    formulation.addObjective(variables=variables, coefficients=np.ones(cols), sense='minimize')

    if solver == 'gurobi':
        formulation.setSolver('gurobi')
        [problem, interface] = formulation.interpretProblem()
        problem.optimize()
        if problem.status == 3:
            status = False
        else:
            status = True
    
    elif solver == 'cplex':
        formulation.setSolver('cplex')
        [problem, interface] = formulation.interpretProblem()
        problem.solve()
        if 'infeasible' in problem.solution.get_status_string():
            status = False
        else:
            status = True
               
    return status

def checkCycles(sources, targets, solver):
    model = createRegulatoryModel(sources, targets)
    status = solveForCycles(model, solver)
    return status

def constructLayer(genes, regulatory_dict, solver):
    layer_dict = {}
    layer_dict['sources'] = []
    layer_dict['targets'] = []
    layer_dict['interaction'] = []
    layer_dict['genes'] = genes
    for gene in genes:
        if gene in regulatory_dict['target_ENSEMBL']:
            # get the index of all the target genes that match to gene
            target_indices = [i for i, x in enumerate(regulatory_dict['target_ENSEMBL']) if x == gene]
            [layer_dict['sources'].append(regulatory_dict['source_ENSEMBL'][i]) for i in target_indices]
            [layer_dict['targets'].append(regulatory_dict['target_ENSEMBL'][i]) for i in target_indices]
            [layer_dict['interaction'].append(regulatory_dict['interaction'][i]) for i in target_indices]
            layer_dict['genes'] = list(set(layer_dict['genes']).union(set(layer_dict['sources'])))
    hasCycles = checkCycles(layer_dict['sources'], layer_dict['targets'], solver)
    if not hasCycles:
        return layer_dict
        
def constructLayerWithCycles(genes, regulatory_dict, solver):
    layer_dict = {}
    layer_dict['sources'] = []
    layer_dict['targets'] = []
    layer_dict['interaction'] = []
    layer_dict['genes'] = genes
    for gene in genes:
        if gene in regulatory_dict['target_ENSEMBL']:
            # get the index of all the target genes that match to gene
            target_indices = [i for i, x in enumerate(regulatory_dict['target_ENSEMBL']) if x == gene]
            [layer_dict['sources'].append(regulatory_dict['source_ENSEMBL'][i]) for i in target_indices]
            [layer_dict['targets'].append(regulatory_dict['target_ENSEMBL'][i]) for i in target_indices]
            [layer_dict['interaction'].append(regulatory_dict['interaction'][i]) for i in target_indices]
            layer_dict['genes'] = list(set(layer_dict['genes']).union(set(layer_dict['sources'])))
    return layer_dict
    
def stackLayers(genes, regulatory_dict, num_layers, solver):
    if num_layers == 0:
        return {'sources': [], 'targets': [], 'interaction': [], 'genes': []}
    for i in range(num_layers):
        layer = constructLayer(genes, regulatory_dict, solver)
        if layer:
            genes = layer['genes']
            safeLayer = layer.copy()
        else:
            break
    if not layer:
        if i == 0:
            safeLayer = {'sources': [], 'targets': [], 'interaction': [], 'genes': []}
        return safeLayer
    else:
        layer['layer_hasCycles'] = [i + 1, False]
    return layer

def stackLayersWithCycles(genes, regulatory_dict, num_layers, solver):
    if num_layers == 0:
        return {'sources': [], 'targets': [], 'interaction': [], 'genes': []}
    for i in range(num_layers):
        layer = constructLayerWithCycles(genes, regulatory_dict, solver)
        if layer:
            genes = layer['genes']
            safeLayer = layer.copy()
        else:
            break
    if not layer:
        if i == 0:
            safeLayer = {'sources': [], 'targets': [], 'interaction': [], 'genes': []}
        return safeLayer
    else:
        layer['layer_hasCycles'] = [i + 1, False]
    return layer

def parseRule(rule):
    # fisrt step is checking if the rule is written in operands, if not replace words with operands
    parsedRule = rule.replace('and', '&').replace('or', '|')
    # second step is to split the rule into a list of operands and operators
    parsedRule = parsedRule.split(' ')
    shiftedRuleOff = ''
    for each in parsedRule:
        sufix = 'off'
        if each == '=':
            shiftedRuleOff = shiftedRuleOff + each + ' '
        elif each == '&':
            shiftedRuleOff = shiftedRuleOff + '|' + ' '
            continue
        elif each == '|':
            shiftedRuleOff = shiftedRuleOff + '&' + ' '
            continue
        else:
            backParentesisSplit = each.split(')', maxsplit=1)
            # negated operrand check
            negOp = backParentesisSplit[0].split('!')
            if len(negOp) > 1:
                backParentesisSplit[0] = negOp[0] + negOp[1]
                sufix = 'on'
            onString = backParentesisSplit[0] + (f'_{sufix}')
            if len(backParentesisSplit) > 1:
                onString = onString + ')' + backParentesisSplit[1]
            shiftedRuleOff = shiftedRuleOff + onString + ' '    
    shiftedRuleON = ''
    for each in parsedRule:
        sufix = 'on'
        if each == '&' or each == '|' or each == '=':
            shiftedRuleON = shiftedRuleON + each + ' '
            continue
        else:
            backParentesisSplit = each.split(')',  maxsplit=1)
            # negated operrand check
            negOp = backParentesisSplit[0].split('!')
            if len(negOp) > 1:
                backParentesisSplit[0] = negOp[0] + negOp[1]
                sufix = 'off'            
            onString = backParentesisSplit[0] + (f'_{sufix}')
            if len(backParentesisSplit) > 1:
                onString = onString + ')' + backParentesisSplit[1]
            shiftedRuleON = shiftedRuleON + onString + ' '
    shiftedRuleOff =  shiftedRuleOff.rstrip()
    shiftedRuleON = shiftedRuleON.rstrip()
    #if '(' not in shiftedRuleOff and 'KI' in shiftedRuleOff and 'KO' in shiftedRuleOff:
    #    shiftedRuleOff = shiftedRuleOff.replace('&', '|')
    return [shiftedRuleOff, shiftedRuleON]

def parseRuleDossage(rule):
    # fisrt step is checking if the rule is written in operands, if not replace words with operands
    parsedRule = rule.replace('and', '&').replace('or', '|')
    # second step is to split the rule into a list of operands and operators
    parsedRule = parsedRule.split(' ')
    shiftedRuleOff = ''
    for each in parsedRule:
        sufix = 'off'
        if each == '=':
            shiftedRuleOff = shiftedRuleOff + each + ' '
        elif each == '&':
            shiftedRuleOff = shiftedRuleOff + '|' + ' '
            continue
        elif each == '|':
            shiftedRuleOff = shiftedRuleOff + '&' + ' '
            continue
        else:
            backParentesisSplit = each.split(')', maxsplit=1)
            # negated operrand check
            negOp = backParentesisSplit[0].split('!')
            if len(negOp) > 1:
                backParentesisSplit[0] = negOp[0] + negOp[1]
                sufix = 'on'
            onString = backParentesisSplit[0] + (f'_{sufix}')
            if len(backParentesisSplit) > 1:
                onString = onString + ')' + backParentesisSplit[1]
            shiftedRuleOff = shiftedRuleOff + onString + ' '    
    shiftedRuleON = ''
    for each in parsedRule:
        sufix = 'on'
        if each == '&' or each == '|' or each == '=':
            shiftedRuleON = shiftedRuleON + each + ' '
            continue
        else:
            backParentesisSplit = each.split(')',  maxsplit=1)
            # negated operrand check
            negOp = backParentesisSplit[0].split('!')
            if len(negOp) > 1:
                backParentesisSplit[0] = negOp[0] + negOp[1]
                sufix = 'off'            
            onString = backParentesisSplit[0] + (f'_{sufix}')
            if len(backParentesisSplit) > 1:
                onString = onString + ')' + backParentesisSplit[1]
            shiftedRuleON = shiftedRuleON + onString + ' '
    shiftedRuleOff =  shiftedRuleOff.rstrip()
    shiftedRuleON = shiftedRuleON.rstrip()
    return [shiftedRuleOff, shiftedRuleON]



def updateModelFromAssign(model, product, reactant):
    name = reactant
    reaction = cobra.core.Reaction(name)
    reaction.name = name
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    reaction.add_metabolites({cobra.core.Metabolite(product, name=product): 1.0})
    metabolite = cobra.core.Metabolite(reactant, name=reactant)
    reaction.add_metabolites({metabolite: -1.0})
    model.add_reactions([reaction])
    return model

def reactionsToModel(gpr, model):
    [product, reactants] = gpr.split(' = ')
    reactants = reactants.replace('&', 'and').replace('|', 'or')
    boolRule = parse(reactants, mode='eval')
    equality = None
    if isinstance(boolRule.body, Name):
        model = updateModelFromAssign(model, product, reactants)
        equality = reactants
    else:
        model = buildReactionFromBranch(boolRule.body, model, 0, cobra.Metabolite(product))
    return model, equality

def addOutputReactions(model):
    name = 'gpr_on'
    gpr_on = cobra.core.Reaction(name)
    gpr_on.name = name
    gpr_on.lower_bound = 0
    gpr_on.upper_bound = 10000
    gpr_on.add_metabolites({model.metabolites.get_by_id(name): -1.0})
    return [gpr_on]

def checkRegulatoryNetwork(regulatoryDict):
    # Banned symbols
    bannedSymbols = ['*', '+', '/', ',', '~', "'", '-', ':', ';']
    # Implement a bidirectional dictionary as a more efficient way to store the reference
    referenceDict = bidict()
    # Check if the regulatory network has valid sources and targets, raise warning if not and suppress the character
    for i, zipped in enumerate(zip(regulatoryDict['source_ENSEMBL'], regulatoryDict['target_ENSEMBL'])):
        source, target = zipped
        if any(symbol in source for symbol in bannedSymbols) or source[0].isdigit() or 'and' in source or 'or' in source:
            #flag = f'{source} contains banned symbols, removing them, string has been modified'
            #warnings.warn(flag)
            # Encode the string to remove banned symbols
            if source in referenceDict.values():
                identifier = referenceDict.inverse[source]
            else:
                identifier = f'M_I_{str(uuid.uuid4())}'.replace('-', '')
            referenceDict[identifier] = source
            regulatoryDict['source_ENSEMBL'][i] = identifier
        if any(symbol in target for symbol in bannedSymbols) or target[0].isdigit() or 'and' in target or 'or' in target:
            #flag = f'{target} contains banned symbols, removing them, string has been modified'
            #warnings.warn(flag)
            # Encode the string to remove banned symbols
            if target in referenceDict.values():
                identifier = referenceDict.inverse[target]
            else:
                identifier = f'M_I_{str(uuid.uuid4())}'.replace('-', '')
               
            regulatoryDict['target_ENSEMBL'][i] = identifier
        
    return regulatoryDict, referenceDict

def transformReactionIdsIntoIndexes(gDict, reactionList):
    newGDict = defaultdict(list)
    for key, value in gDict.items():
        reactionIndexes = []
        for reaction in value:
            reactionIndexes.append(reactionList.index(reaction))
        newGDict[key] = reactionIndexes
    return newGDict


# Decode a string
def decode_string(string: str, referenceDict: dict):
    list_strings = string.split("_K")
    if len(list_strings) > 1:
        encoded_string = list_strings[0]
        last_name = list_strings[1]
        last_name = f"_K{last_name}"
    else:
        encoded_string = string
        last_name = ""
    decoded_string = referenceDict[encoded_string]
    decoded_string = decoded_string + last_name
    return decoded_string
