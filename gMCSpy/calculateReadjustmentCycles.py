from .calculateGMCSwithCycles import getGPRDict, stackLayersWithCycles, checkRegulatoryNetwork, bidict
from .MinimalCutSetClass import tqdm
from .calculateGMCS import re, initiateLogger
from .ProblemDefinitions import readjustmentProblem
from copy import deepcopy
from datetime import datetime
from .Utilities import createLog
from .calculateGMCSwithCycles import boolTransformer, mpbn

def calculateReadjustment(gMISList, model, regulatory_dataframe, num_layers=2, solver='gurobi', **kwargs):
    # Name of the analysis
    name = (
        "gMCSpy_Readjustment_"
        + model.id
        + "_"
        + datetime.now().strftime("%H-%M-%S--%d-%m-%Y")
    )
    
    # Save all the parameters in a log file    
    pathN = "logs/" + name + "/"
    path = kwargs.get('path', pathN)
    
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
    
    if kwargs.get('path', None) is None:
        path = "logs/" + name
        name = ("results")
        handler = name + ".log"
        createLog(path + "/results/" + name + ".log")
        logging = initiateLogger(name + ".log", path + "/results/" + name + ".log")
    else:
        name = ("results.log")
        createLog(path + '/' + name)
        logging = initiateLogger(name, path + name)
        
    kwargs.pop('path', None)
    
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
        if len(gMIS) == 1:
            continue
        calculateReadjustmentSingleGMIs(gMIS, rulesDict, referenceDict, solver, logging, **kwargs)


def calculateReadjustmentSingleGMIs(gMIS, rulesDict, referenceDict, solver='gurobi', logging=None, **kwargs):        
        precomputedDict = {}
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
            results = calculateReadjustmentInRule(key, value['rules'], gMISDict, precomputedDict)
            if results is not None:
                if results['readjustment']:
                    logging.info(results)
                    break
                    

            
def calculateRules(GPR, regulatory_dict, num_layers=2, solver='gurobi', **kwargs):
    key, value = GPR
    expandedRules = stackLayersWithCycles(value['genes'], regulatory_dict, num_layers, solver)
    boolform = boolTransformer(expandedRules, value)
    genes = list(boolform.keys())
    return boolform, genes

def calculateReadjustmentInRule(keyGPR, rules, gMISDict, precomputedDict):
    network = mpbn.MPBooleanNetwork(rules)
    
    for key, items in gMISDict.items():
        interaction = 0 if key.endswith('_KO') else 1
        provitionalList = list(gMISDict.keys())
        
        #print(provitionalList)
        
        #minimum = sum([-1 for i in provitionalList if i.endswith('_KI')])
        #print(minimum)
        provitionalList.remove(key)
        
        '''
        identifier = frozenset([keyGPR, items, interaction])
        if identifier in precomputedDict:
            act = precomputedDict[identifier]
        else:
            act = network.attractors(reachable_from={items: interaction}, limit=10000)
            precomputedDict[identifier] = act
        '''
        
        total = 0
        constraints = {items: interaction}
        for gene in provitionalList:
            if gene.endswith('_KI'):
                gene_without_lastname = gene.split('_KI')[0]
                tInt = 1
            else:
                gene_without_lastname = gene.split('_KO')[0]
                tInt = 0
            constraints[gene_without_lastname] = tInt
        
        names = bidict()
        resultsDict = {}
        for name in gMISDict.keys():
            gene_without_lastname = name.split('_KI')[0] if name.endswith('_KI') else name.split('_KO')[0]
            names[name] = gene_without_lastname
            resultsDict[name] = []
        break
        # Check if configuration exists
    act = network.attractors(constraints=constraints, limit=1)
    try:
        attractors = list(act)[0]
    except:
        gMIS = list(gMISDict.keys())
        return {'readjustment': True, 'gMIS': gMIS, 'constraints': constraints}
    
    '''
    for key, _ in constraints.items():
        name = names.inverse[key]
        resultsDict[name].append(attractors.get(key))	
    
    total = 0
    for key, value in resultsDict.items():
        val = None
        if 'KO' in key:
            val = min(value)
        else:
            if max(value) is not None:
                val = -1*max(value)
        if val is not None:
            total += val            
    
    if total != minimum:
        return {'readjustment': True, 'gMIS': gMISDict, 'total': total, 'minimum': minimum, 'resultsDict': resultsDict}
    '''    
         
    
