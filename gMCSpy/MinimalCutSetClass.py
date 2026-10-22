import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from scipy.sparse import eye, hstack

from .ProblemDefinitions import cobra
from .ProblemDefinitions import List

import logging
import time
import uuid
import warnings
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
def _extractGurobiSolution(subVector, solutionVariableIndices, solutionVariableNames,
                           allVars, mergeSolutions):
    """Single source of truth for turning one Gurobi integer solution into (a) the gene
    set to report and (b) the no-good cut to exclude it and its supersets.

    Used by BOTH the classic ``iterateOverSolver`` loop and the staged-lazy path, so the
    cut/merge semantics cannot diverge -- the only difference between the two paths is how
    the returned cut is registered (``addLConstr`` vs ``cbLazy``), never what it is. This
    is a verbatim extraction of the classic Gurobi branch's per-solution logic, including
    the multi-gene-row reduction: when more KO rows are selected than the merged gene
    count, the rows overlap, and the cut is reduced to the single-gene rows when they
    alone explain the merged set.

    :param subVector: solution values at the solutionVariableIndices positions (zp vars).
    :returns: (solutionIndex, geneSet, cutVars, cutRhs) where geneSet is a
        ``frozenset`` of genes when mergeSolutions else a list of KO-row keys, and the
        no-good cut is ``sum(cutVars) <= cutRhs``.
    """
    solutionIndex = list(i for i in np.where(subVector > 0.5)[0])
    currentSolution = [allVars[solutionVariableIndices[i]] for i in solutionIndex]
    if mergeSolutions:
        geneSolutions = [solutionVariableNames[i] for i in solutionIndex]
        mergedSolutions = frozenset().union(*geneSolutions)
        solutionLength = len(mergedSolutions)
        if len(solutionIndex) > solutionLength:
            auxGeneList = [gene for gene in geneSolutions if len(gene) == 1]
            newSolutionIndex = [solutionIndex[geneSolutions.index(gene)]
                                for gene in auxGeneList]
            if len(auxGeneList) == solutionLength:
                currentSolution = [allVars[solutionVariableIndices[i]]
                                   for i in newSolutionIndex]
            # else: keep the full selection as the cut (possible minimality occurrence)
        geneSet = mergedSolutions
    else:
        geneSet = [solutionVariableNames[i] for i in solutionIndex]
    return solutionIndex, geneSet, currentSolution, len(currentSolution) - 1


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
        else:
            companionMatrix = companionMatrix.tolil(copy=True)

        if np.any(lb > 0):
            warnings.warn(
                "Reactions with a positive lower bound are treated as irreversible "
                "(lb = 0); the forced flux is not part of the gMCS problem."
            )

        negativeLbColumns = matrixS[:, splittingIndex]
        companionNegativeLbColumns = companionMatrix[:, splittingIndex]

        # A reaction that cannot run forward (ub <= 0, e.g. closed lb = ub = 0 or
        # reverse-only) has no forward direction. Zeroing its forward column makes that
        # dual row vacuous, which is the same as removing the primal variable, while
        # keeping column j = reaction j for the target index used downstream.
        noForward = np.where(ub <= 0)[0]
        matrixS = matrixS.tolil(copy=True)
        matrixS[:, noForward] = 0
        companionMatrix[:, noForward] = 0

        matrixS = hstack((matrixS, -negativeLbColumns))
        companionMatrix = hstack((companionMatrix, companionNegativeLbColumns))
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
        # Pool solutions that differ only in optional zero-weight rows map to the same gene
        # set; each gene set is reported once.
        seenSolutions = set()

        def warnTimeLimit():
            message = ("A solve hit the time limit: the gMCS found so far are a LOWER BOUND, "
                       "not a complete enumeration.")
            logger.warning(message)
            warnings.warn(message)

        if solver == "gurobi":
            OPTIMALITY = 2
            if forceLength:
                rhsForceLength = problemInterpreted.getConstrByName('forceLength').RHS
            while(solutionLength <= MAX_SOLUTION_LENGTH and solutionCountPerLength <= MAX_SOLUTIONS and rhsForceLength <= MAX_SOLUTION_LENGTH):
                startTime = time.time()
                problemInterpreted.optimize()
                endTime = time.time()
                if problemInterpreted.Status == interface.GRB.TIME_LIMIT:
                    warnTimeLimit()
                numberOfSolutions = problemInterpreted.SolCount
                if numberOfSolutions > 0:
                    solutions = {}
                    for j in range(numberOfSolutions):
                        problemInterpreted.Params.SolutionNumber = j
                        solutions[j] = np.array(problemInterpreted.getAttr("Xn", problemInterpreted.getVars()))

                    for j in range(numberOfSolutions):
                        solutionIndex, geneSet, cutVars, cutRhs = _extractGurobiSolution(
                            solutions[j][solutionVariableIndices],
                            solutionVariableIndices, solutionVariableNames,
                            problemInterpreted.getVars(), mergeSolutions)
                        # Preserve the classic side effect: in merge mode the loop-control
                        # solutionLength tracks the current merged order (was inline before
                        # the cut was added).
                        if mergeSolutions:
                            solutionLength = len(geneSet)
                        # Prevent this solution (and its supersets) from being found again.
                        newConExpr = interface.LinExpr(np.ones(len(cutVars)), cutVars)
                        problemInterpreted.addLConstr(
                            newConExpr,
                            interface.GRB.LESS_EQUAL,
                            cutRhs,
                            f"{str(uuid.uuid4())}_MultipleSol",
                        )

                        if frozenset(geneSet) in seenSolutions:
                            continue
                        seenSolutions.add(frozenset(geneSet))

                        solutionCountPerLength = solutionCountPerLength + 1
                        solutionInformation = {}
                        solutionInformation["solution"] = geneSet
                        if removeRelated:
                            variablesToBound = problem.getContainedVariables(solutionInformation["solution"])
                            for name in variablesToBound:
                                problemInterpreted.getVarByName(name).setAttr('UB', 0)
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
                if problemInterpreted.solution.get_status() in (
                    problemInterpreted.solution.status.MIP_time_limit_feasible,
                    problemInterpreted.solution.status.MIP_time_limit_infeasible,
                ):
                    warnTimeLimit()
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
                        if frozenset(solutionInformation["solution"]) in seenSolutions:
                            continue
                        seenSolutions.add(frozenset(solutionInformation["solution"]))
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
                if notSkip and problemInterpreted.getStatus() == "timelimit":
                    warnTimeLimit()

     
                if problemInterpreted.getStatus() == "optimal" and notSkip:
                    solution = []
                    variables = []
                    for j in problemInterpreted.getVars():
                        if problemInterpreted.getVal(j) > 0.5 and 'zp_' in j.name:
                            pos = int(j.name.split('_')[1])
                            variables.append(pos)
                            solution.append(solutionVariableNames[solutionVariableIndices.index(pos)])

                    # Cut on the single-gene rows when they explain the whole gene set, as the
                    # gurobi/cplex branches do, so that gene supersets are excluded too.
                    cutVariables = variables
                    if mergeSolutions:
                        singleGeneVariables = [v for v, row in zip(variables, solution) if len(row) == 1]
                        if len(singleGeneVariables) == len(frozenset().union(*solution)):
                            cutVariables = singleGeneVariables
                    coefficients = [1 for _ in cutVariables]
                    problem.addConstraint(f"{str(uuid.uuid4())}_MultipleSol", cutVariables, coefficients, "L", len(cutVariables) - 1)
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
    def iterateOverSolverStagedLazy(problem, problemInterpreted, solver, solutionVariables,
                                    bounds, mergeSolutions=True, verbose=0, handler=None,
                                    minKOLength=1, exact=False, progressEvery=100):
        """Staged lazy enumeration of gMCS (Gurobi only).

        One lazy branch-and-bound tree per gene-count order ``k = minKOLength ..
        MAX_SOLUTION_LENGTH``, with a cardinality constraint on ``sum(objective)``. Two
        consequences:

          * The objective is bounded to the current order, so the solver spends little
            effort proving the bound -- the dominant cost of the classic per-solution
            re-solve loop on hard models. The search is driven by a ``MIPSOL`` callback
            that records each solution and rejects it (and its supersets) with ``cbLazy``.
          * Enumerating orders 1..k in turn with prior stages' no-good cuts retained makes
            every reported solution minimal by construction; no post-filter is required.

        The per-solution cut/merge logic is the shared :func:`_extractGurobiSolution`, so
        the gene sets and no-good cuts are identical to the classic ``iterateOverSolver``.

        :param exact: If False (default), stage k is bounded by ``sum(objective) <= k``.
            This can never over-constrain: any gMCS of order <= k stays feasible, so
            completeness is guaranteed. If True, stage k is fixed to ``sum(objective) == k``
            which explores a tighter tree but relies on the weighted objective equalling
            the merged gene count for every solution -- faster on some models, but the
            loose bound is the safe default when exactness is paramount.
        :returns: ``(solutionDict, completed)`` with the same solutionDict shape as
            ``iterateOverSolver``. ``completed`` is ``False`` if any stage hit the solver
            time limit -- in that case the set is a lower bound, never a full enumeration.
        """
        if solver != "gurobi":
            raise ValueError("iterateOverSolverStagedLazy is implemented for gurobi only")

        logger = logging.getLogger(handler)
        optModel, interface = problemInterpreted[0], problemInterpreted[1]
        GRB = interface.GRB
        MAX_SOLUTION_LENGTH = bounds["MAX_SOLUTION_LENGTH"]
        MAX_SOLUTIONS = bounds["MAX_SOLUTIONS"]
        solutionVariableIndices = solutionVariables["indices"]
        solutionVariableNames = solutionVariables["names"]

        # cbLazy requires this; set it defensively so the method is correct even if the
        # caller did not put LazyConstraints in the solver parameters.
        optModel.Params.LazyConstraints = 1
        optModel.update()
        allVars = optModel.getVars()
        zpVars = [allVars[i] for i in solutionVariableIndices]

        # Cardinality expression == the objective (gene-count weighting), reusing exactly
        # what the forceLength feature builds. The upper bound caps the current stage's
        # order; the lower bound (exact mode only) additionally pins it. Loose (<=k) can
        # never exclude a valid gMCS of order <= k, so it is the safe default.
        objective = problem.getObjective()
        cardExpr = interface.LinExpr([float(c) for c in objective["coefficients"]],
                                     [allVars[i] for i in objective["variables"]])
        upper = optModel.addConstr(cardExpr <= minKOLength, name="stagedLazyUpper")
        lower = optModel.addConstr(cardExpr >= minKOLength, name="stagedLazyLower") if exact else None
        optModel.update()

        solutionDict = {}
        solutionCount = 0
        completed = True
        runStart = time.time()
        # progress state shared with the callback: last report time and count, for the
        # running rate. A single line every progressEvery solutions gives a tail-able
        # heartbeat during the long order-k stages.
        progress = {"lastTime": runStart, "lastCount": 0}

        for k in range(minKOLength, MAX_SOLUTION_LENGTH + 1):
            upper.setAttr("RHS", k)
            if lower is not None:
                lower.setAttr("RHS", k)
            stageSolutions = []
            baseCount = solutionCount

            def callback(model, where, _stage=stageSolutions, _base=baseCount, _k=k):
                if where != GRB.Callback.MIPSOL:
                    return
                subVector = np.array(model.cbGetSolution(zpVars))
                _, geneSet, cutVars, cutRhs = _extractGurobiSolution(
                    subVector, solutionVariableIndices, solutionVariableNames,
                    allVars, mergeSolutions)
                _stage.append((geneSet, cutVars, cutRhs))
                model.cbLazy(interface.quicksum(cutVars) <= cutRhs)

                running = _base + len(_stage)
                if progressEvery and running % progressEvery == 0:
                    now = time.time()
                    rate = (running - progress["lastCount"]) / max(1e-9, now - progress["lastTime"])
                    line = (f"{running:>8} gMCS  {(now - runStart) / 60:7.1f} min  "
                            f"({rate:.2f}/s)  order={_k}")
                    logger.info(line)
                    if verbose > 0:
                        print(line, flush=True)
                    progress["lastTime"] = now
                    progress["lastCount"] = running

            startTime = time.time()
            optModel.optimize(callback)
            elapsed = time.time() - startTime
            status = optModel.Status

            # Make this stage's cuts permanent so later stages still exclude these
            # solutions and their supersets, independent of cbLazy carrying across solves.
            for geneSet, cutVars, cutRhs in stageSolutions:
                optModel.addLConstr(
                    interface.LinExpr(np.ones(len(cutVars)), cutVars),
                    GRB.LESS_EQUAL, cutRhs, f"{str(uuid.uuid4())}_MultipleSol")
            optModel.update()

            for geneSet, cutVars, cutRhs in stageSolutions:
                solLength = len(geneSet)
                solutionDict[solutionCount] = {
                    "solution": geneSet,
                    "order": "order" + str(solLength) + "_" + str(solutionCount),
                }
                solutionCount += 1

            stageLine = (f"  order {k}: {len(stageSolutions)} gMCS  "
                         f"{(time.time() - runStart) / 60:.1f} min cumulative  "
                         f"status={status}")
            logger.info(stageLine)
            logger.info(f"stagedLazyOrder{k}: {elapsed} count {len(stageSolutions)} "
                        f"status {status}")
            if verbose > 0:
                print(stageLine, flush=True)

            if status in (GRB.TIME_LIMIT, GRB.INTERRUPTED):
                completed = False
                logger.info(f"stagedLazy INCOMPLETE: order {k} hit the time limit")
                break
            if solutionCount >= MAX_SOLUTIONS:
                break

        return solutionDict, completed

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
