import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from scipy.sparse import eye, hstack

from .ProblemDefinitions import cobra
from .ProblemDefinitions import List

import logging
import time
import uuid
from tqdm import tqdm


"""
The `MinimalCutSetProblem` class is a Python class that represents a minimal cut set problem. 
It provides methods for preparing constraint matrices, building indicators into the problem,
iterating over a solver, and getting corresponding variables.

Example Usage:
    # Create an instance of the MinimalCutSetProblem class
    problem = MinimalCutSetProblem()

    # Prepare constraint matrices
    cobraModel = ...
    constraintMatrices = problem.prepareConstraintMatrices(cobraModel)

    # Build indicators into the problem
    prefix = ...
    variables = ...
    indicatorVariables = ...
    sense = ...
    constant = ...
    activeWhen = ...
    problem.buildIndicatorsIntoProblem(prefix, variables, indicatorVariables, sense, constant, activeWhen)

    # Iterate over a solver
    problemInterpreted = ...
    solver = ...
    solutionVariables = ...
    bounds = ...
    solutionDict = problem.iterateOverSolver(problemInterpreted, solver, solutionVariables, bounds)

    # Get corresponding variables
    realVariables = ...
    interpretedSolutions = problem.getCorrespondingVariables(solutionDict, realVariables)

Main functionalities:
- Preparing constraint matrices for a given cobra model
- Building indicators into the problem
- Iterating over a solver to find solutions
- Getting corresponding variables for the solutions

Methods:
- prepareConstraintMatrices(cobraModel, companionType=None, companionMatrix=None): Prepares constraint matrices for a given cobra model.
- buildIndicatorsIntoProblem(problem, prefix, variables, indicatorVariables, sense, constant, activeWhen): Builds indicators into the problem.
- iterateOverSolver(problem, problemInterpreted, solver, solutionVariables, bounds, mergeSolutions=False, earlyStop=False, removeRelated=False, verbose=0, handler=None, iterationProgress=False): Iterates over a solver to find solutions.
- getCorrespondingVariables(problem, solutions, realVariables, mergeSolutions=False): Gets corresponding variables for the solutions.

Fields:
- forceLength: A boolean field indicating whether to force a specific solution length.
- verbose: An integer field indicating the verbosity level.
- solver: A string field indicating the solver to use.
- constantB: A float field representing a constant value.
- constantC: A float field representing a constant value.
- cobraModel: A field representing a cobra model.
"""
class MinimalCutSetProblem:
    def __init__(self, forceLength=False, verbose=0, solver=None, B=1e-3, C=1e-3) -> None:
        self.forceLength = forceLength
        self.verbose = verbose
        self.solver = solver
        self.constantB = B
        self.constantC = C
        self.cobraModel = None

    @staticmethod
    def prepareConstraintMatrices(cobraModel, companionType=None, companionMatrix=None) -> dict:
        if cobraModel is None:
            raise Exception("No cobra model provided")

        if companionMatrix is None and companionType is None:
            raise Exception("No companion matrix provided, no companion type provided, to build K use companionType='K', to use G provide companionMatrix")

        matrixS = create_stoichiometric_matrix(cobraModel, array_type="lil")

        [_, reactions] = matrixS.shape

        lb, ub = zip(*[(r._lower_bound, r._upper_bound) for r in cobraModel.reactions])
        lb = np.array(lb)
        ub = np.array(ub)

        # Splitting index
        splittingIndex = np.where(lb < 0)[0]

        # Creating companion matrix if K is requested
        if companionType == 'K':
            companionMatrix = eye(reactions, format="lil")

        negativeLbColumns = matrixS[:, splittingIndex]
        matrixS = hstack((matrixS, -negativeLbColumns))

        negativeLbColumns = companionMatrix[:, splittingIndex]
        companionMatrix = hstack((companionMatrix, negativeLbColumns))
        companionMatrix = companionMatrix.tolil()

        constraintMs = {}
        constraintMs["S"] = matrixS
        constraintMs["Companion"] = companionMatrix
        return constraintMs

    @staticmethod
    def buildIndicatorsIntoProblem(problem, prefix, variables, indicatorVariables, sense, constant, activeWhen):
        for j, item in enumerate(indicatorVariables):
            name = f"{prefix}{j}"
            problem.addIndicator(
                name= name,
                indicatorVariable= item,
                activeWhen= activeWhen,
                variable= variables[j],
                coefficients= 1,
                sense= sense,
                b= constant
            )

    @staticmethod
    def iterateOverSolver(problem, problemInterpreted, solver, solutionVariables, bounds,
                          mergeSolutions=False, earlyStop=False, removeRelated=False, verbose=0, handler=None,
                          iterationProgress=False):
        TIMELIMIT = problem.getTimeLimitParameter()
        logger = logging.getLogger(handler)
        interface = problemInterpreted[1]
        problemInterpreted = problemInterpreted[0]
        solutionLength = 0
        solutionCountPerLength = 0
        solutionCount = 0
        solutionDict = {}
        rhsForceLength = 0
        MAX_SOLUTION_LENGTH = bounds['MAX_SOLUTION_LENGTH']
        MAX_SOLUTIONS = bounds['MAX_SOLUTIONS']
        forceLength = bounds['forceLength']
        solutionVariableIndices = solutionVariables['indices']
        solutionVariableNames = solutionVariables['names']

        if solver == "gurobi":
            OPTIMALITY = 2
            if forceLength:
                rhsForceLength = problemInterpreted.getConstrByName('forceLength').RHS
            while(solutionLength <= MAX_SOLUTION_LENGTH and solutionCountPerLength <= MAX_SOLUTIONS and rhsForceLength <= MAX_SOLUTION_LENGTH):
                startTime = time.time()
                problemInterpreted.optimize()
                endTime = time.time()
                numberOfSolutions = problemInterpreted.SolCount
                if numberOfSolutions > 0:
                    solutions = {}
                    for j in range(numberOfSolutions):
                        problemInterpreted.Params.SolutionNumber = j
                        solutions[j] = np.array(problemInterpreted.getAttr("Xn", problemInterpreted.getVars()))

                    for j in range(numberOfSolutions):
                        currentSolution =[]
                        solutionIndex = list(i for i in np.where(solutions[j][solutionVariableIndices] > 0.5)[0])
                        currentSolution = [
                            problemInterpreted.getVars()[solutionVariableIndices[i]]
                            for i in solutionIndex
                        ]
                        names = [var.VarName for var in currentSolution]
                        # To correctly prevent the same solution from being found again, we add a constraint that the sum of the variables in the solution must be less than the total number of variables of order 1 in the solution
                        if mergeSolutions:
                            geneSolutions = [solutionVariableNames[i] for i in solutionIndex]
                            mergedSolutions = frozenset().union(*geneSolutions)
                            solutionLength = len(mergedSolutions)
                            if len(solutionIndex) <= solutionLength:
                                pass
                            else:
                                auxGeneList = []
                                for gene in geneSolutions:
                                    if len(gene) == 1:
                                        auxGeneList.append(gene)
                                newSolutionIndex = [solutionIndex[geneSolutions.index(gene)] for gene in auxGeneList]
                                if len(auxGeneList) == solutionLength:
                                    currentSolution = [problemInterpreted.getVars()[solutionVariableIndices[i]] for i in newSolutionIndex]
                                else: 
                                    pass
                                    #print("Possible minimality occurance")
                                    #print("Individual solutions: ", geneSolutions)
                                    #print("Individual genes found: ", auxGeneList)
                                    #print('##################################')
                            
                            
                        newConExpr = interface.LinExpr(np.ones(len(currentSolution)), currentSolution)
                        problemInterpreted.addLConstr(
                            newConExpr,
                            interface.GRB.LESS_EQUAL,
                            len(currentSolution) - 1,
                            f"{str(uuid.uuid4())}_MultipleSol",
                        )

                        solutionCountPerLength = solutionCountPerLength + 1
                        currentSolution = [solutionVariableNames[i] for i in solutionIndex]
                        solutionInformation = {}

                        solutionInformation["solution"] = currentSolution
                        if mergeSolutions:
                            solutionInformation["solution"] = frozenset().union(*solutionInformation["solution"])
                        if removeRelated:
                            variablesToBound = problem.getContainedVariables(solutionInformation["solution"])
                            afectedVars = 0
                            for name in variablesToBound:
                                problemInterpreted.getVarByName(name).setAttr('UB', 0)
                                afectedVars += 1
                        solLength = len(solutionInformation["solution"])
                        keyName = "order" + str(solLength) + "_" + str(solutionCountPerLength)
                        solutionInformation["order"] = keyName
                        solutionDict[solutionCount] = solutionInformation
                        if  solLength > solutionLength:
                            solutionLength = len(solutionInformation["solution"])
                        solutionCount += 1

                    timeTaken = endTime - startTime
                    if verbose > 0:
                        print(f"Iteration yielded {numberOfSolutions} solutions of length: {solLength}")
                        logger.info(f"timeOrder{solLength}: {timeTaken}")

                    if earlyStop:
                        if (solutionLength == MAX_SOLUTION_LENGTH) & (abs(timeTaken - TIMELIMIT) > 15 ):
                            logger.info("Optimality reached, early stop")
                            break

                else:
                    if forceLength:
                        rhsForceLength = rhsForceLength + 1
                        problemInterpreted.getConstrByName('forceLength').setAttr("RHS", rhsForceLength)
                    else:
                        break
        elif solver == "cplex":
            OPTIMALITY = 0
            if forceLength:
                rhsForceLength = problemInterpreted.linear_constraints.get_rhs('forceLength')
            while(solutionLength <= MAX_SOLUTION_LENGTH and solutionCountPerLength <= MAX_SOLUTIONS and rhsForceLength <= MAX_SOLUTION_LENGTH):
                startTime = time.time()
                problemInterpreted.populate_solution_pool()
                endTime = time.time()
                numberOfSolutions = problemInterpreted.solution.pool.get_num()
                if numberOfSolutions > 0:
                    solutions = {}
                    for j in range(numberOfSolutions):
                        solutions[j] = np.array(problemInterpreted.solution.pool.get_values(j))

                    for j in range(numberOfSolutions):
                        currentSolution =[]
                        solutionIndex = list(i for i in np.where(solutions[j][solutionVariableIndices] > 0.5)[0])
                        currentSolution = [solutionVariableIndices[i] for i in solutionIndex]
                        names = problemInterpreted.variables.get_names(currentSolution)
                        keyName = "order" + str(len(solutionIndex)) + "_" + str(solutionCountPerLength)
                        solutionCountPerLength = solutionCountPerLength + 1
                        
                         # To correctly prevent the same solution from being found again, we add a constraint that the sum of the variables in the solution must be less than the total number of variables of order 1 in the solution
                        if mergeSolutions:
                            geneSolutions = [solutionVariableNames[i] for i in solutionIndex]
                            mergedSolutions = frozenset().union(*geneSolutions)
                            solutionLength = len(mergedSolutions)
                            if len(solutionIndex) <= solutionLength:
                                pass
                            else:
                                auxGeneList = []
                                for gene in geneSolutions:
                                    if len(gene) == 1:
                                        auxGeneList.append(gene)
                                newSolutionIndex = [solutionIndex[geneSolutions.index(gene)] for gene in auxGeneList]
                                if len(auxGeneList) == solutionLength:
                                    currentSolution = [solutionVariableIndices[i] for i in newSolutionIndex]
                                else: 
                                    print("Possible minimality occurance")
                                    print("Individual solutions: ", geneSolutions)
                                    print("Individual genes found: ", auxGeneList)
                                    print('##################################')
                        
                        currentSolutionNames = [solutionVariableNames[i] for i in solutionIndex]
                        solutionInformation = {}
                        solutionInformation["solution"] = currentSolutionNames
                        constraintLength = len(currentSolution) - 1 
                        
                        if mergeSolutions:
                            solutionInformation["solution"] = frozenset().union(*solutionInformation["solution"])
                            
                        problemInterpreted.linear_constraints.add(
                        lin_expr=[
                            interface.SparsePair(
                                ind=currentSolution, val=np.ones(len(currentSolution))
                            )
                        ],
                        rhs=[constraintLength],
                        senses=["L"],
                        names=[f"{'-'.join(names)}_MultipleSol"])
                        if removeRelated:
                            variablesToBound = problem.getContainedVariables(solutionInformation["solution"])
                            for name in variablesToBound:
                                problemInterpreted.variables.set_upper_bounds(name, 0)
                        keyName = "order" + str(len(solutionInformation["solution"])) + "_" + str(solutionCountPerLength)
                        solutionInformation["order"] = keyName
                        solutionDict[solutionCountPerLength] = solutionInformation
                        if len(solutionInformation["solution"]) > solutionLength:
                            solutionLength = len(solutionInformation["solution"])
                        solutionCount += 1


                    if verbose > 0:
                        print(f"Iteration yielded {numberOfSolutions} solutions of length: {len(solutionInformation['solution'])}")
                        timeTaken = endTime - startTime
                        logger.info(f"timeOrder{len(solutionInformation['solution'])}: {timeTaken}")
                    
                    if earlyStop:
                        if (solutionLength == MAX_SOLUTION_LENGTH) & (abs(timeTaken - TIMELIMIT) > 15 ):
                            logger.info("Optimality reached, early stop")
                            break
                        
                else:
                    if forceLength:
                        rhsForceLength = rhsForceLength + 1
                        problemInterpreted.linear_constraints.set_rhs('forceLength', rhsForceLength)
                    else:
                        break

        elif solver == "scip":
            if forceLength:
                rhsForceLength = problem.getConstraintByName('forceLength')['rhs']

            while(solutionLength <= MAX_SOLUTION_LENGTH and solutionCountPerLength <= MAX_SOLUTIONS and rhsForceLength <= MAX_SOLUTION_LENGTH):
                notSkip = True
                startTime = time.time()
                try:
                    problemInterpreted.optimize()
                except:
                    print("Error in SCIP")
                    notSkip = False
                endTime = time.time()

     
                if problemInterpreted.getStatus() == "optimal" and notSkip:
                    solution = []
                    variables = []
                    for j in problemInterpreted.getVars():
                        if problemInterpreted.getVal(j) > 0.5 and 'zp_' in j.name:
                            pos = int(j.name.split('_')[1])
                            variables.append(pos)
                            solution.append(solutionVariableNames[solutionVariableIndices.index(pos)])
                            coefficients = [1 for i in variables]

                    problem.addConstraint(f"{str(uuid.uuid4())}_MultipleSol", variables, coefficients, "L", len(variables) - 1)
                    [problemInterpreted, interface] = problem.interpretProblem(verbose=0)
                    #problemInterpreted.writeProblem(f"model_scip.lp")
                    keyName = "order_" + str(len(solution)) + "_" + str(solutionCountPerLength)
                    solutionCountPerLength = solutionCountPerLength + 1
                    solutionInformation = {}
                    solutionInformation["order"] = keyName
                    solutionInformation["solution"] = solution
                    if mergeSolutions:
                        solutionInformation["solution"] = frozenset().union(*solutionInformation["solution"])
                    solLength = len(solutionInformation["solution"])
                    solutionDict[solutionCountPerLength] = solutionInformation
                    if  solLength > solutionLength:
                        solutionLength = len(solutionInformation["solution"])

                    if verbose > 0:
                        timeTaken = endTime - startTime
                        print(f"Solution of length {solLength} found in {timeTaken} seconds")
                        logger.info(f"timeOrder{solLength}: {timeTaken}")

                else:
                    if forceLength:
                        rhsForceLength = rhsForceLength + 1
                        problem.getConstraintByName('forceLength')['rhs'] = rhsForceLength
                        [problemInterpreted, interface] = problem.interpretProblem(verbose=0)
                    else:
                        break
        return solutionDict

    @staticmethod
    def getCorrespondingVariables(problem, solutions, realVariables, mergeSolutions=False):
        variables = list(problem.getVariableNames().keys())

        zpVariables = [variables[var] for var in problem.getObjective()['variables']]


        #The variables that correspond to the problem variables, e.g the real variables, reactions, genes, Knockouts, etc. 
        interpretedSolutions = {}
        #Match the position of the zp_i variables in the solution
        for j, key in enumerate(solutions):
            solution= solutions[key]
            indices = [zpVariables.index(variable) for variable in solution]
            #Get the corresponding real variables
            solutionInformation = {}
            solutionInformation["order"] = key
            solutionInformation["solution"] = [realVariables[index] for index in indices]
            if mergeSolutions:
                solutionInformation["solution"] = frozenset().union(*solutionInformation["solution"])
            solutionInformation["length"] = len(solution)
            interpretedSolutions[j] = solutionInformation


        return interpretedSolutions
    
    

def prepareModelNutrientGeneMCS(cobraModel: cobra.Model, exchangeReactions: List):
            
    # Perform a reaction splitting, we are only interested in the exchange reactions that bring metabolites into the system (inputs)
    
    # Check what reactions are reversible
    if exchangeReactions is None:
        listOfGeneRules = []
        for reaction in cobraModel.exchanges:
            if len(reaction.metabolites) > 1:
                print("WARNING: Some exchange reactions exchange more than one metabolite.")
            name  = list(reaction.metabolites.keys())[0].id + "_ag"
            if reaction.reversibility:
                # Split the reaction into two reactions, one for the forward direction and one for the reverse direction
                newReaction = cobra.Reaction(reaction.id + "_reverse")
                newReaction.name = reaction.id + "_reverse"
                newReaction.lower_bound = 0
                newReaction.upper_bound = 1000
                for metabolite, value in reaction.metabolites.items():
                    newReaction.add_metabolites({metabolite: -value})     
                
                # Check if reaction is input, if from nothing yields a product, then it is an input
                if len(newReaction.reactants) == 0 and len(newReaction.products) == 1:
                    newReaction.gene_reaction_rule = name
                    listOfGeneRules.append(newReaction.id)
                
                # Last step add reaction to model 
                cobraModel.add_reactions([newReaction])
                
                reaction.lower_bound = 0
                if len(reaction.reactants) == 0 and len(reaction.products) == 1:
                    reaction.gene_reaction_rule = name
                    listOfGeneRules.append(reaction.id)
            else:
                if len(reaction.reactants) == 0 and len(reaction.products) == 0:
                    reaction.gene_reaction_rule = name
                    listOfGeneRules.append(reaction.id)
                    



    print('===================================================================')
    
    
    return cobraModel
