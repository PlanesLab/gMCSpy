from .calculateSyntheticDosageGMCS import mergeIsforms, getGPRDict, stackLayers, checkRegulatoryNetwork, interpretExpandedRules, bidict, uuid
from .MinimalCutSetClass import tqdm
from .calculateGMCS import re, initiateLogger
from .ProblemDefinitions import readjustmentProblem
from copy import deepcopy
from datetime import datetime
from .Utilities import createLog

def calculateReadjustment(gMISList, model, regulatory_dataframe, num_layers=2, solver='gurobi', **kwargs):
    
    # Name of the analysis
    name = (
        "gMCSpy_Readjustment_"
        + model.id
        + "_"
        + datetime.now().strftime("%H-%M-%S--%d-%m-%Y")
    )
    
    # Save all the parameters in a log file    
    path = "logs/" + name + "/"
    createLog(path)
    configName = ("configuration")
    logging = initiateLogger(configName + ".log", path + configName + ".log")
    # Save the parameters in the log file
    logging.info(f"solver: {solver}, num_layers: {num_layers}, solver: {solver}")
    logging.info(f"Model: {model.id}")
    logging.info(f"gMIS to check: {gMISList}")
    # close logger
    handler = logging.handlers[0]
    handler.close()
    logging.removeHandler(handler)
    
    
    ### Results logger
    
    path = "logs/" + name + "/results/"
    createLog(path)
    path = "logs/" + name
    name = ("results")
    handler = name + ".log"
    logging = initiateLogger(name + ".log", path + "/results/" + name + ".log")
    
    
    # Prepare the regulatory network
    regulatory_dict = regulatory_dataframe.to_dict('list')
    # Check if the regulatory network has valid sources and targets
    regulatory_dict, referenceDict = checkRegulatoryNetwork(regulatory_dict)
    # create a dictionary with the gpr rules associated to each reaction to get a unique list of gprs
    gprDict = getGPRDict(model)
    
    # Expand all the rules to the desired number of layers
    rulesDict = {}
    i = 0
    for key, value in tqdm(gprDict.items()):
        GPR = key, value
        rules, genes = calculateRules(GPR, regulatory_dict, num_layers, solver, **kwargs)
        rulesDict[i] = {'genes': genes, 'rules': rules}
        i += 1
    
    print('checking gMIS')
    for gMIS in tqdm(gMISList):
        calculateReadjustmentSingleGMIs(gMIS, rulesDict, referenceDict, solver, logging, **kwargs)


def calculateReadjustmentSingleGMIs(gMIS, rulesDict, referenceDict, solver='gurobi', logging=None, **kwargs):        
        # check if the gMIS is in the reference dictionary
        gMISDict = bidict()
        encodedgMIS = []
        for gene in gMIS:
            if '_KO' in gene:
                geneT = gene.split('_KO')[0]
            else:
                geneT = gene.split('_KI')[0]
            if geneT in referenceDict.values():
                encodedgMIS.append(referenceDict.inverse[geneT])
                gMISDict[gene] = referenceDict.inverse[geneT]
            else:
                encodedgMIS.append(geneT)
                gMISDict[gene] = geneT
        
        filteredRules = {}
        for gene in encodedgMIS:
            for key, value in rulesDict.items():
                ruleGenes = value['genes']
                if set([gene]).issubset(set(ruleGenes)):
                    filteredRules[key] = value
        
        #print(gMIS)
        for key, value in filteredRules.items():
            results = calculateReadjustmentInRule(value['rules'], gMISDict, solver)
            if results is not None:
                if results['readjustment']:
                    logging.info(results)
                    break
                    

            
def calculateRules(GPR, regulatory_dict, num_layers=2, solver='gurobi', **kwargs):
    key, value = GPR
    expandedRules = stackLayers(value['genes'], regulatory_dict, num_layers, solver)
    expandedRules = interpretExpandedRules(key, expandedRules, value['genes'])
    genes = list(expandedRules.keys())    
    return expandedRules, genes

def calculateReadjustmentInRule(rules, gMISDict, solver='gurobi'):
    for id, rule in rules.items():
        name, eq = rule['rule'].split(' = ')
        if check_format(eq):
            match_ki = re.search(r'\b\w*_KI\b', eq).group()
            match_ko = re.search(r'\b\w*_KO\b', eq).group()
            name_ki = match_ki.split('_KI')[0]
            transformed = f"({name_ki}_input | {match_ki}) & {match_ko}"
            rules[id]['rule'] = f"{name} = {transformed}"
    rules = ruleTransformation(rules)
    simpGMIs = []
    for key, value in gMISDict.items():
        if '_KO' in key:
            simpGMIs.append(f"{value}_KO")
        elif '_KI' in key:
            simpGMIs.append(f"{value}_KI")
            
    opt = readjustmentProblem(rules, simpGMIs)
    varibleNames = {i.split('--')[0] : int(i.split('--')[1]) for i in opt.getVariableNames()}
    objective = opt.getObjective()
    resultsDict = {i: {'obj': None} for i in simpGMIs}
    
    for intervention in simpGMIs:
        objectiveVariables = deepcopy(objective['variables'])
        objectiveCoefficients = deepcopy(objective['coefficients'])
        typeIntervention = intervention.split('_K')[1]
        readptation = False
        readpCoef = 0
        for index in objectiveVariables:
            name = opt.getVariable(index)['name'].split('--')[0]
            if f"{name}_KI" in simpGMIs and f"{name}_KI" != intervention:
                readpCoef += 1
        
        new_opt = deepcopy(opt)
        for key, value in varibleNames.items():
            if '_KO' in key and key != intervention:
                new_opt.modifyVariable(value, lb=1, ub=1)
            elif '_KI' in key and key != intervention:
                new_opt.modifyVariable(value, lb=0, ub=0)
            elif '_KO' in key and key == intervention:
                new_opt.modifyVariable(value, lb=0, ub=0)
            elif '_KI' in key and key == intervention:
                new_opt.modifyVariable(value, lb=1, ub=1)
        
        gMIS = intervention.split('_K')[0]
        variableIndex = varibleNames.get(gMIS)
        if variableIndex is not None:
            objtIndex = objectiveVariables.index(variableIndex)
            objectiveVariables.pop(objtIndex)
            objectiveCoefficients.pop(objtIndex)
        new_opt.modifyObjective(variables=objectiveVariables, coefficients=objectiveCoefficients)
        new_opt.setSolver(solver)
        problem, interface = new_opt.interpretProblem() 
        if solver == 'gurobi':
            problem.optimize()
            objVal = problem.ObjVal
        elif solver == 'cplex':
            problem.solve()
            objVal = problem.solution.get_objective_value()
        elif solver == 'scip':
            problem.optimize()
            objVal = problem.getObjVal
            
        if typeIntervention == 'O':
            if objVal != -readpCoef:
                readptation = True
        elif typeIntervention == 'I':
            if objVal != -readpCoef:
                readptation = True
        resultsDict[intervention]['obj'] = objVal
        resultsDict[intervention]['readjustment'] = readptation
        resultsDict['readjustment'] = readptation
        if readptation:
            return resultsDict    

def check_format(input_string):
    pattern = r'.*_KO \| .*_KI$'
    return bool(re.match(pattern, input_string))

def ruleTransformation(expandedRules):
    intermediatesDict = bidict()
    newDict = {}
    for key, value in expandedRules.items():
        if bool(re.search(' or ', value['rule']) or re.search(' and ', value['rule'])):
            newRule = value['rule'].replace(' or ', ' | ').replace(' and ', ' & ')
        else: 
            newRule = value['rule']
            
        if bool(re.search(' \| ', newRule) and re.search(' & ', newRule)):
            while bool(re.search(' \| ', newRule) and re.search(' & ', newRule)):    
                matches = extract_parentheses(newRule)
                for match in matches:
                    hasBooleanOp = bool(re.search(' \| ', match) or re.search(' & ', match))
                    if not hasBooleanOp:
                        intermediate = match.replace('(', '').replace(')', '')
                        idPartial = intermediate.replace('!', '')
                        listIntermdiates = frozenset([idPartial] + ['!'])
                        if listIntermdiates not in intermediatesDict:
                            id = str(uuid.uuid4())
                            intermediatesDict[listIntermdiates] = id
                        else:
                            id = intermediatesDict[listIntermdiates]
                        
                    else:
                        if bool(re.search(' \| ', match)):
                            intermediate = match.replace('(', '').replace(')', '')
                            listIntermdiates = frozenset(intermediate.split(' \| ') + ['|'])
                            if listIntermdiates not in intermediatesDict:
                                id = str(uuid.uuid4())
                                intermediatesDict[listIntermdiates] = id
                            else:
                                id = intermediatesDict[listIntermdiates]
                        elif bool(re.search(' & ', match)):
                            intermediate = match.replace('(', '').replace(')', '')
                            listIntermdiates = frozenset(intermediate.split(' & ') + ['&'])
                            if listIntermdiates not in intermediatesDict:
                                id = str(uuid.uuid4())
                                intermediatesDict[listIntermdiates] = id
                            else:
                                id = intermediatesDict[listIntermdiates]
                        
                    dictEntry = match.replace('(', '').replace(')', '')
                    newRule = newRule.replace(match, id)
                    #print(intermediate)
                    #print(dictEntry)
                    dictEntry = f"{id} = {dictEntry}"
                    newDict[id] = {'rule': dictEntry}
            
                newDict[key] = {'rule': newRule}      
        else:
            newDict[key] = value
    return newDict
            

def extract_parentheses(expression):
    # Use regex to find all substrings enclosed in parentheses
    pattern = r'\([^()]*\)'
    matches = re.findall(pattern, expression)
    return matches