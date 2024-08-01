from .ProblemInterpreters import np
from .Utilities import beartype
import scipy
import cobra
from beartype.typing import List
import math

from .MinimalCutSetClass import MinimalCutSetProblem
from .OptimizationProblem import OptimizationProblem


def buildDictionaryMCSProblem(
    cobraModel: cobra.Model,
    forceLength: bool,
    rxnToForce: List[str],
    verbose: int = 0,
    targetB: float = 1e-3,
    cConstant: float = 1e-3,
):
    """
    Builds the dictionary that defines the MILP problem to be solved by the MCS algorithm.

    :param modelStruct: Cobra Model, the metabolic model to be operated on. Only Cobra models are supported.
    :param forceLength: boolean. If true, the length of the solution is forced to be off a given length each iteration, otherwise the length is not forced.
    :param targetB: float. Desired activity level of the metabolic task to be disrupted.
    :param cConstant: float.  Used to activate w variable
    :param rxnToForce: list of strings. Reactions to be forced to be part of the solution. Empty by default.

    If the bounds are left as None they are ment to be infinity (None = Infinity), so each solver can interpret it as it is defined.

    :return: problemDict. Dictionary that defines the MILP problem to be solved by any solver defined in gMCSpy.
    """

    problemMethods = MinimalCutSetProblem()
    problem = OptimizationProblem()

    # Alpha. Used to relate the lower bound of v variables with z variables
    ALPHA = 1

    matrices = problemMethods.prepareConstraintMatrices(cobraModel, "K")

    matrixS = matrices["S"]
    matrixK = matrices["Companion"]

    # S matrix changed size
    [metabolites, reactions] = matrixS.shape

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = objectiveVector.index(1)

    # Target vector
    t = np.zeros(reactions)
    t[indexC] = 1

    if rxnToForce:
        if verbose > 0:
            print("Forcing reactions to be part of the solution: " + str(rxnToForce))
        # Find the index of the reactions in the model
        rxnIndex = [
            cobraModel.reactions.index(cobraModel.reactions.get_by_id(rxn))
            for rxn in rxnToForce
        ]
        matrixK = matrixK[rxnIndex, :]
        reactionsK = matrixK.shape[0]

    # Add variables to the problem

    vIndex = {}
    vIndex["u"] = list(range(0, metabolites))
    vIndex["vp"] = list(range(vIndex["u"][-1] + 1, vIndex["u"][-1] + reactionsK + 1))
    vIndex["w"] = [vIndex["vp"][-1] + 1]
    vIndex["zp"] = list(range(vIndex["w"][-1] + 1, vIndex["w"][-1] + reactionsK + 1))
    vIndex["zw"] = [vIndex["zp"][-1] + 1]
    name = "u_"
    for u in vIndex["u"]:
        problem.addVariable(
            index=u, prefix=name, lb="infinity", ub="infinity", vtype="continuous"
        )

    name = "vp_"
    for vp in vIndex["vp"]:
        problem.addVariable(
            index=vp, prefix=name, lb=0, ub="infinity", vtype="continuous"
        )

    name = "w_"
    w = vIndex["w"][0]
    problem.addVariable(index=w, prefix=name, lb=0, ub="infinity", vtype="continuous")

    name = "zp_"
    count = 0
    for zp in vIndex["zp"]:
        problem.addVariable(
            index=zp, prefix=name, lb=0, ub=1, vtype="binary", relatedTo=count
        )
        count += 1

    name = "zw_"
    zw = vIndex["zw"][0]
    problem.addVariable(index=zw, prefix=name, lb=0, ub=1, vtype="binary")

    vGroups = {}
    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]

    cIndex = {}

    cIndex["Ndual"] = list(range(0, reactions))
    cIndex["forceBioCons"] = [cIndex["Ndual"][-1] + 1]
    cIndex["forceLength"] = [cIndex["forceBioCons"][-1] + 1]

    matrixS = matrixS.transpose().tocsr()
    matrixK = matrixK.transpose().tocsr()
    tIndexConstraint = np.where(t != 0)[0]
    for i in cIndex["Ndual"]:
        constraintName = f"Ndual_constraint{i}"
        variableIndicesU = matrixS[i].indices
        vars = []
        vars = [vIndex["u"][j] for j in variableIndicesU]
        variableIndicesVp = matrixK[i].indices
        vars.extend([vIndex["vp"][j] for j in variableIndicesVp])
        coefs = np.concatenate((matrixS[i].data, matrixK[i].data))
        problem.addConstraint(
            name=constraintName, variables=vars, coefficients=coefs, sense="G", rhs=0
        )

    constraintName = f"Ndual_constraint{tIndexConstraint[0]}"
    constraint = problem.getConstraint(constraintName)
    constraint["variables"].extend([w])
    coefs = np.concatenate((constraint["coefficients"], -t[tIndexConstraint]))
    problem.modifyConstraint(
        constraintName=constraintName,
        variables=constraint["variables"],
        coefficients=coefs,
        sense=constraint["sense"],
        rhs=constraint["rhs"],
    )

    problem.addConstraint(
        name="forceBioCons",
        variables=[w],
        coefficients=[-targetB],
        sense="E",
        rhs=-cConstant,
    )

    if forceLength:
        variables = vIndex["zp"]
        problem.addConstraint(
            name="forceLength",
            variables=variables,
            coefficients=np.ones(len(variables)),
            sense="E",
            rhs=1,
        )

    problem.addObjective(vIndex["zp"], np.ones(len(vIndex["zp"])), sense="minimize")

    # Indicators
    # z = 1  -->  v >= alpha
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOn_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=1,
        sense="G",
        constant=ALPHA,
    )

    # z = 0  -->  v <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOff_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=0,
        sense="L",
        constant=0,
    )

    return problem


def buildDictionaryMCSTargetedProblem(
    cobraModel: cobra.Model,
    forceLength: bool,
    rxnSubset: List[str],
    rxnKnockOut: List[str],
    verbose: int = 0,
    targetB: float = 1e-3,
    cConstant: float = 1e-3,
):
    """
    Builds the dictionary that defines the MILP problem to be solved by the MCS algorithm.

    :param modelStruct: Cobra Model, the metabolic model to be operated on. Only Cobra models are supported.
    :param forceLength: boolean. If true, the length of the solution is forced to be off a given length each iteration, otherwise the length is not forced.
    :param targetB: float. Desired activity level of the metabolic task to be disrupted.
    :param cConstant: float.  Used to activate w variable
    :param rxnToForce: list of strings. Reactions to be forced to be part of the solution. Empty by default.

    If the bounds are left as None they are ment to be infinity (None = Infinity), so each solver can interpret it as it is defined in each solver.

    :return: problemDict. Dictionary that defines the MILP problem to be solved by any solver defined in gMCSpy.
    """

    problemMethods = MinimalCutSetProblem()
    problem = OptimizationProblem()

    # Alpha. Used to relate the lower bound of v variables with z variables
    ALPHA = 1e-3

    matrices = problemMethods.prepareConstraintMatrices(cobraModel, "K")

    matrixS = matrices["S"]
    matrixK = matrices["Companion"]

    # S matrix changed size
    [mUpdated, nUpdated] = matrixS.shape

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = objectiveVector.index(1)

    # Target vector
    t = np.zeros(nUpdated)
    t[indexC] = 1

    if rxnSubset:
        if verbose > 0:
            print("Forcing reactions to be part of the solution: " + str(rxnSubset))
        # Find the index of the reactions in the model
        rxnIndex = [
            cobraModel.reactions.index(cobraModel.reactions.get_by_id(rxn))
            for rxn in rxnSubset
        ]
        # Order the reactions to be forced by the index of the reactions in the model
        rxnIndex = sorted(rxnIndex)
        matrixK = matrixK[rxnIndex, :]
        n = matrixK.shape[0]

    rxnKnockOutIndex = np.zeros(len(cobraModel.reactions))
    rxnKnockOutIndex[
        [
            cobraModel.reactions.index(cobraModel.reactions.get_by_id(rxn))
            for rxn in rxnKnockOut
        ]
    ] = 1

    vIndex = {}
    vGroups = {}
    vIndex["u"] = list(range(0, mUpdated))
    vIndex["vp"] = list(range(vIndex["u"][-1] + 1, vIndex["u"][-1] + n + 1))
    vIndex["w"] = [vIndex["vp"][-1] + 1]
    vIndex["zp"] = list(range(vIndex["w"][-1] + 1, vIndex["w"][-1] + n + 1))
    vIndex["zw"] = [vIndex["zp"][-1] + 1]
    vIndex["epsp"] = list(range(vIndex["zw"][-1] + 1, vIndex["zw"][-1] + n + 1))
    vIndex["espw"] = [vIndex["epsp"][-1] + 1]
    vIndex["delp"] = list(range(vIndex["espw"][-1] + 1, vIndex["espw"][-1] + n + 1))
    vIndex["delw"] = [vIndex["delp"][-1] + 1]
    vIndex["x"] = list(range(vIndex["delw"][-1] + 1, vIndex["delw"][-1] + nUpdated + 2))
    nVars = vIndex["x"][-1] + 1

    for u in vIndex["u"]:
        problem.addVariable(
            index=u, prefix="u_", lb="infinity", ub="infinity", vtype="continuous"
        )

    for vp in vIndex["vp"]:
        problem.addVariable(
            index=vp, prefix="vp_", lb=0, ub="infinity", vtype="continuous"
        )

    for w in vIndex["w"]:
        problem.addVariable(
            index=w, prefix="w_", lb=0, ub="infinity", vtype="continuous"
        )

    for zp in vIndex["zp"]:
        problem.addVariable(index=zp, prefix="zp_", lb=0, ub=1, vtype="binary")

    for zw in vIndex["zw"]:
        problem.addVariable(index=zw, prefix="zw_", lb=0, ub=1, vtype="binary")

    for epsp in vIndex["epsp"]:
        problem.addVariable(
            index=epsp, prefix="epsp_", lb=0, ub="infinity", vtype="continuous"
        )

    for espw in vIndex["espw"]:
        problem.addVariable(index=espw, prefix="espw_", lb=0, ub=0, vtype="continuous")

    for delp in vIndex["delp"]:
        problem.addVariable(
            index=delp, prefix="delp_", lb=0, ub="infinity", vtype="continuous"
        )

    for delw in vIndex["delw"]:
        problem.addVariable(index=delw, prefix="delw_", lb=0, ub=0, vtype="continuous")

    for x in vIndex["x"][:-1]:
        problem.addVariable(
            index=x, prefix="x_", lb=0, ub="infinity", vtype="continuous"
        )
    problem.addVariable(
        index=vIndex["x"][-1], prefix="x_", lb=1000, ub="infinity", vtype="continuous"
    )

    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]
    vGroups["eps"] = vIndex["epsp"] + vIndex["espw"]
    vGroups["del"] = vIndex["delp"] + vIndex["delw"]

    cIndex = {}

    cIndex["Ndual"] = list(range(0, nUpdated))
    cIndex["forceBioCons"] = [cIndex["Ndual"][-1] + 1]
    cIndex["forceKnockOut"] = [cIndex["forceBioCons"][-1] + 1]

    numberOfLinearCombinations = matrixS.shape[0] + matrixK.shape[0] + 1 + 1
    cIndex["linearCombination"] = list(
        range(
            cIndex["forceKnockOut"][-1] + 1,
            cIndex["forceKnockOut"][-1] + numberOfLinearCombinations,
        )
    )

    cIndex["forceLength"] = [cIndex["linearCombination"][-1] + 1]
    nCons = cIndex["forceLength"][-1] + 1

    A = scipy.sparse.lil_matrix((nCons, nVars))
    A[np.ix_(cIndex["Ndual"], vIndex["u"])] = matrixS.transpose()
    A[np.ix_(cIndex["Ndual"], vIndex["vp"])] = matrixK.transpose()
    A[np.ix_(cIndex["Ndual"], vIndex["w"])] = -t
    A[np.ix_(cIndex["forceBioCons"], vIndex["w"])] = -targetB
    A[np.ix_(cIndex["forceKnockOut"], vIndex["vp"])] = rxnKnockOutIndex

    auxS = scipy.sparse.hstack(
        [matrixS, scipy.sparse.lil_matrix(np.zeros((mUpdated, 1)))]
    )
    auxK = scipy.sparse.hstack([matrixK, scipy.sparse.lil_matrix(np.zeros((n, 1)))])
    auxT = scipy.sparse.hstack([scipy.sparse.lil_matrix(-t), scipy.sparse.lil_matrix(targetB)])
    A[np.ix_(cIndex["linearCombination"], vIndex["x"])] = scipy.sparse.vstack(
        [auxS, auxK, auxT]
    )

    matrixEp1 = scipy.sparse.lil_matrix(
        (mUpdated, (len(vIndex["vp"]) + len(vIndex["w"])))
    )
    matrixEp2 = scipy.sparse.eye((len(vIndex["vp"]) + len(vIndex["w"])), format="lil")
    auxEp = scipy.sparse.vstack([matrixEp1, -matrixEp2])
    A[np.ix_(cIndex["linearCombination"], vGroups["eps"])] = auxEp

    A[np.ix_(cIndex["linearCombination"], vGroups["del"])] = -auxEp

    rhs = np.zeros(nCons)
    rhs[cIndex["Ndual"]] = None
    rhs[cIndex["forceBioCons"][0]] = -cConstant
    rhs[cIndex["forceKnockOut"]] = None
    rhs[
        cIndex["linearCombination"][mUpdated : len(cIndex["linearCombination"]) - 1]
    ] = rxnKnockOutIndex

    lhs = np.zeros(nCons)
    lhs[cIndex["Ndual"]] = 0
    lhs[cIndex["forceBioCons"][0]] = -1
    lhs[cIndex["forceKnockOut"]] = 1e-3
    lhs[
        cIndex["linearCombination"][mUpdated : len(cIndex["linearCombination"]) - 1]
    ] = rxnKnockOutIndex

    sense = np.array(["E"] * nCons)
    sense[cIndex["Ndual"]] = "G"
    sense[cIndex["forceBioCons"][0]] = "G"
    sense[cIndex["forceKnockOut"]] = "G"

    A = A.tocsr()
    problem.addConstraintsFromSparseMatrix(A, lhs, sense)

    varsForceBioCons = A[cIndex["forceBioCons"][0]].indices
    coefficients = A[cIndex["forceBioCons"][0]].data
    constraintName = f'c{cIndex["forceBioCons"][0]}'
    problem.modifyConstraint(
        constraintName=constraintName,
        variables=varsForceBioCons,
        coefficients=coefficients,
        sense="L",
        rhs=-cConstant,
    )
    # problem.addConstraint(name="forceBioCons_2",variables=varsForceBioCons, coefficients=coefficients, sense="L", rhs=-cConstant)

    problem.addObjective(vIndex["zp"], np.ones(len(vIndex["zp"])), sense="minimize")

    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]
    vGroups["eps"] = vIndex["epsp"] + vIndex["espw"]
    vGroups["del"] = vIndex["delp"] + vIndex["delw"]

    # Indicators
    # z = 1  -->  v >= alpha
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOn_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        sense="G",
        constant=ALPHA,
        activeWhen=1,
    )

    # z = 0  -->  v <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOff_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        sense="L",
        constant=0,
        activeWhen=0,
    )

    # z = 1  -->  epsilon <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOnEpsilon_",
        variables=vGroups["eps"],
        indicatorVariables=vGroups["z"],
        sense="L",
        constant=0,
        activeWhen=1,
    )

    # z = 0  -->  epsilon <= M
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOffEpsilon_",
        variables=vGroups["eps"],
        indicatorVariables=vGroups["z"],
        sense="L",
        constant=1e3,
        activeWhen=0,
    )

    return problem


@beartype
def buildDictionaryFBA(cobraModel: cobra.Model):
    """Build FBA problem dictionary
    The basic idea of FBA is to use constraints based on the mass balance and stoichiometry of reactions to determine the
    steady-state fluxes through the metabolic network. The steady-state assumption assumes that the concentrations of metabolites
    do not change over time, which is reasonable for many metabolic systems.
    """

    # build Stoichiometric matrix
    matrixS = cobra.util.array.create_stoichiometric_matrix(
        cobraModel, array_type="lil"
    )

    [numberMetabolites, numberReactions] = matrixS.shape
    ## m = number of metabolites
    ## n = number of reactions

    # reactions bounds
    lb, ub = zip(*[(r._lower_bound, r._upper_bound) for r in cobraModel.reactions])
    lb = np.array(lb)
    ub = np.array(ub)

    problem = OptimizationProblem()

    variableNames = {}
    vIndices = {}

    vIndices["u"] = list(range(0, numberReactions))
    nVars = vIndices["u"][-1] + 1

    prefix = "u_"
    for u in vIndices["u"]:
        problem.addVariable(
            index=u, prefix=prefix, lb=lb[u], ub=ub[u], vtype="continuous"
        )

    cIndices = {}

    cIndices["S"] = list(range(0, numberMetabolites))
    nCons = cIndices["S"][-1] + 1
    matrixS = matrixS.tocsr()

    for i in cIndices["S"]:
        constraintName = f"Srow_{i}"
        variablesIndices = matrixS[i].indices
        coefs = matrixS[i].data
        problem.addConstraint(
            name=constraintName,
            variables=variablesIndices,
            coefficients=coefs,
            sense="E",
            rhs=0,
        )

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = [objectiveVector.index(1)]

    problem.addObjective(indexC, np.ones(len(indexC)), sense="maximize")

    return problem


def buildDictionaryFBAWithDeletions(
    cobraModel: cobra.Model, gMCS: List[str] = None, isoformSeparator: str = None
):
    """Build FBA problem dictionary
    The basic idea of FBA is to use constraints based on the mass balance and stoichiometry of reactions to determine the
    steady-state fluxes through the metabolic network. The steady-state assumption assumes that the concentrations of metabolites
    do not change over time, which is reasonable for many metabolic systems.
    """

    # build Stoichiometric matrix
    matrixS = cobra.util.array.create_stoichiometric_matrix(
        cobraModel, array_type="lil"
    )

    [numberMetabolites, numberReactions] = matrixS.shape
    ## m = number of metabolites
    ## n = number of reactions

    if isoformSeparator is not None:
        # Get the GPR of each reaction
        GPRs = [
            [r, r.gpr]
            for r in cobraModel.reactions
            for gene in r.gpr.genes
            if gene.split(isoformSeparator)[0] in gMCS
        ]
        # modelReactions = cobraModel.reactions
        # modelGenes = [g.id for g in cobraModel.genes]
        reactionsToKnockOut = []
        for reaction, gpr in GPRs:
            if len(gpr.genes) == 1:
                reactionsToKnockOut.append(reaction.id)
            else:
                expr = gpr.to_string()
                for gene in gpr.genes:
                    # Substitute the gene name for True in the expression
                    geneSplit = gene.split(isoformSeparator)[0]
                    if geneSplit in gMCS:
                        expr = expr.replace(gene, "False")
                    else:
                        expr = expr.replace(gene, "True")
                # Evaluate the expression
                outcome = eval(expr)
                if not outcome:
                    reactionsToKnockOut.append(reaction.id)

    else:
        # Get the GPR of each reaction
        GPRs = [
            [r, r.gpr]
            for r in cobraModel.reactions
            for gene in gMCS
            if gene in r.gpr.genes
        ]
        # modelReactions = cobraModel.reactions
        # modelGenes = [g.id for g in cobraModel.genes]
        reactionsToKnockOut = []
        for reaction, gpr in GPRs:
            if len(gpr.genes) == 1:
                reactionsToKnockOut.append(reaction.id)
            else:
                expr = gpr.to_string()
                for gene in gpr.genes:
                    # Substitute the gene name for True in the expression
                    if gene in gMCS:
                        expr = expr.replace(gene, "False")
                    else:
                        expr = expr.replace(gene, "True")
                # Evaluate the expression
                outcome = eval(expr)
                if not outcome:
                    reactionsToKnockOut.append(reaction.id)

    """     if isoformSeparator is not None:
        for gene in cobraModel.genes:
            gene.id = gene.id.split(isoformSeparator)[0]
            if gene.id in gMCS:
                for reaction in gene.reactions:
                    reactionsToKnockOut.append(reaction.id)
    else:
        for gene in cobraModel.genes:
            if gene.id in gMCS:
                for reaction in gene.reactions:
                    reactionsToKnockOut.append(reaction.id) """

    reactionsToKnockOut = list(set(reactionsToKnockOut))

    # reactions bounds
    lb = []
    ub = []
    for reaction in cobraModel.reactions:
        if reaction.id in reactionsToKnockOut:
            lb.append(0)
            ub.append(0)
        else:
            lb.append(reaction.lower_bound)
            ub.append(reaction.upper_bound)

    lb = np.array(lb)
    ub = np.array(ub)

    problem = OptimizationProblem()

    variableNames = {}
    vIndices = {}

    vIndices["u"] = list(range(0, numberReactions))
    nVars = vIndices["u"][-1] + 1

    prefix = "u_"
    for u in vIndices["u"]:
        problem.addVariable(
            index=u, prefix=prefix, lb=lb[u], ub=ub[u], vtype="continuous"
        )

    cIndices = {}

    cIndices["S"] = list(range(0, numberMetabolites))
    nCons = cIndices["S"][-1] + 1
    matrixS = matrixS.tocsr()

    for i in cIndices["S"]:
        constraintName = f"Srow_{i}"
        variablesIndices = matrixS[i].indices
        coefs = matrixS[i].data
        problem.addConstraint(
            name=constraintName,
            variables=variablesIndices,
            coefficients=coefs,
            sense="E",
            rhs=0,
        )

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = [objectiveVector.index(1)]

    problem.addObjective(indexC, np.ones(len(indexC)), sense="maximize")

    return problem

def buildDictionaryReactionsFBAWithDeletions(
    cobraModel: cobra.Model, reactions: List[str] = None
):
    """Build FBA problem dictionary
    The basic idea of FBA is to use constraints based on the mass balance and stoichiometry of reactions to determine the
    steady-state fluxes through the metabolic network. The steady-state assumption assumes that the concentrations of metabolites
    do not change over time, which is reasonable for many metabolic systems.
    """

    # build Stoichiometric matrix
    matrixS = cobra.util.array.create_stoichiometric_matrix(
        cobraModel, array_type="lil"
    )

    [numberMetabolites, numberReactions] = matrixS.shape
    ## m = number of metabolites
    ## n = number of reactions

    reactionsToKnockOut = reactions

    # reactions bounds
    lb = []
    ub = []
    for reaction in cobraModel.reactions:
        if reaction.id in reactionsToKnockOut:
            lb.append(0)
            ub.append(0)
        else:
            lb.append(reaction.lower_bound)
            ub.append(reaction.upper_bound)

    lb = np.array(lb)
    ub = np.array(ub)

    problem = OptimizationProblem()

    variableNames = {}
    vIndices = {}

    vIndices["u"] = list(range(0, numberReactions))
    nVars = vIndices["u"][-1] + 1

    prefix = "u_"
    for u in vIndices["u"]:
        problem.addVariable(
            index=u, prefix=prefix, lb=lb[u], ub=ub[u], vtype="continuous"
        )

    cIndices = {}

    cIndices["S"] = list(range(0, numberMetabolites))
    nCons = cIndices["S"][-1] + 1
    matrixS = matrixS.tocsr()

    for i in cIndices["S"]:
        constraintName = f"Srow_{i}"
        variablesIndices = matrixS[i].indices
        coefs = matrixS[i].data
        problem.addConstraint(
            name=constraintName,
            variables=variablesIndices,
            coefficients=coefs,
            sense="E",
            rhs=0,
        )

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = [objectiveVector.index(1)]

    problem.addObjective(indexC, np.ones(len(indexC)), sense="maximize")

    return problem

def buildDictionaryMCSGeneProblem(
    cobraModel: cobra.Model,
    gDict: dict,
    matrixG: scipy.sparse.csc_matrix,
    relationships: dict,
    maxKOLength: int,
    numberNewGenesByKO: dict,
    minKOLength: int = 0,
    forceLength: bool = False,
    verbose: int = 0,
    targetB: float = 1e-3,
    cConstant: float = 1e-3,
    logPath: str = None,
):
    """
    Builds the dictionary that defines the MILP problem to be solved by the MCS algorithm.

    :param modelStruct: Cobra Model, the metabolic model to be operated on. Only Cobra models are supported.
    :param forceLength: boolean. If true, the length of the solution is forced to be off a given length each iteration, otherwise the length is not forced.
    :param targetB: float. Desired activity level of the metabolic task to be disrupted.
    :param cConstant: float.  Used to activate w variable
    :param rxnToForce: list of strings. Reactions to be forced to be part of the solution. Empty by default.

    If the bounds are left as None they are ment to be infinity (None = Infinity), so each solver can interpret it as it is defined.

    :return: problemDict. Dictionary that defines the MILP problem to be solved by any solver defined in gMCSpy.
    """

    problem = OptimizationProblem()
    problemMethods = MinimalCutSetProblem()

    # Alpha. Used to relate the lower bound of v variables with z variables
    ALPHA = 1

    numberOfPossibleKO = matrixG.shape[0]
    possibleKOs = list(gDict.keys())

    if relationships is not None:
        numberOfRelationships = sum(len(lst) for lst in relationships.values())
    else:
        numberOfRelationships = 0

    matrices = problemMethods.prepareConstraintMatrices(
        cobraModel=cobraModel, companionMatrix=matrixG
    )
    matrixS = matrices["S"]
    matrixG = matrices["Companion"]

    # S matrix changed size
    [metabolites, reactions] = matrixS.shape

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = objectiveVector.index(1)

    # Target vector
    t = np.zeros(reactions)
    t[indexC] = 1

    # Build variables
    vGroups = {}
    vIndex = {}

    vIndex["u"] = list(range(0, metabolites))
    vIndex["vp"] = list(
        range(vIndex["u"][-1] + 1, vIndex["u"][-1] + numberOfPossibleKO + 1)
    )
    vIndex["w"] = [vIndex["vp"][-1] + 1]
    vIndex["zp"] = list(
        range(vIndex["w"][-1] + 1, vIndex["w"][-1] + numberOfPossibleKO + 1)
    )
    vIndex["zw"] = [vIndex["zp"][-1] + 1]

    name = "u_"
    for u in range(0, metabolites):
        problem.addVariable(
            index=u, prefix=name, lb="infinity", ub="infinity", vtype="continuous"
        )

    name = "vp_"
    for i, vp in enumerate(range(u + 1, u + numberOfPossibleKO + 1)):
        problem.addVariable(
            index=vp,
            prefix=name,
            lb=0,
            ub="infinity",
            vtype="continuous",
            relatedTo=possibleKOs[i],
        )

    name = "w_"
    w = vp + 1
    problem.addVariable(index=w, prefix=name, lb=0, ub="infinity", vtype="continuous")

    name = "zp_"
    for i, zp in enumerate(range(w + 1, w + numberOfPossibleKO + 1)):
        problem.addVariable(
            index=zp, prefix=name, lb=0, ub=1, vtype="binary", relatedTo=possibleKOs[i]
        )

    name = "zw_"
    zw = zp + 1
    problem.addVariable(index=zw, prefix=name, lb=0, ub=1, vtype="binary")

    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]

    # Build constraints
    cIndex = {}
    cIndex["Ndual"] = list(range(0, reactions))
    cIndex["forceBioCons"] = [cIndex["Ndual"][-1] + 1]
    if numberOfRelationships > 0:
        cIndex["relationships"] = list(
            range(
                cIndex["forceBioCons"][-1] + 1,
                cIndex["forceBioCons"][-1] + numberOfRelationships + 1,
            )
        )
        cIndex["forceLength"] = [cIndex["relationships"][-1] + 1]
    else:
        cIndex["forceLength"] = [cIndex["forceBioCons"][-1] + 1]

    matrixS = matrixS.transpose().tocsr()
    matrixG = matrixG.transpose().tocsr()
    tIndexConstraint = np.where(t != 0)[0]
    for i in cIndex["Ndual"]:
        constraintName = f"Ndual_constraint{i}"
        variableIndicesU = matrixS[i].indices
        vars = [vIndex["u"][j] for j in variableIndicesU]
        variableIndicesVp = matrixG[i].indices
        vars.extend([vIndex["vp"][j] for j in variableIndicesVp])
        coefs_s = matrixS[i].data
        coefs_g = matrixG[i].data
        coefs = np.concatenate((coefs_s, coefs_g))

        problem.addConstraint(constraintName, vars, coefs, "G", 0)

    constraintName = f"Ndual_constraint{tIndexConstraint[0]}"
    constraint = problem.getConstraint(constraintName)
    constraint["variables"].extend([w])
    coefs = np.concatenate((constraint["coefficients"], -t[tIndexConstraint]))
    problem.modifyConstraint(
        constraintName=constraintName,
        variables=constraint["variables"],
        coefficients=coefs,
        sense=constraint["sense"],
        rhs=constraint["rhs"],
    )

    vpVariableNames = {}
    for index, name in zip(vIndex["zp"], gDict.keys()):
        vpVariableNames[name] = index

    if numberOfRelationships > 0:
        subIndex = 0
        for key, values in relationships.items():
            vpContainer = vpVariableNames[key]
            for value in values:
                constraintName = f"relationship_constraint{subIndex}"
                vpContained = vpVariableNames[value]
                vars = [vpContainer, vpContained]
                coefs = [-1, 1]
                problem.addConstraint(constraintName, vars, coefs, "G", 0)
                subIndex += 1

    # Build objective
    objectiveVars = [
        vIndex["zp"][i]
        for i in range(0, len(possibleKOs))
    ]
    problem.addObjective(objectiveVars, list(numberNewGenesByKO.values()), sense="minimize")

    if forceLength:
        constraintName = "forceLength"
        vars = objectiveVars
        coefs = list(numberNewGenesByKO.values())
        problem.addConstraint(constraintName, vars, coefs, "L", 1)

    constraintName = "forceBioCons_constraint"
    vars = vIndex["w"]
    coefs = [-targetB]
    problem.addConstraint(constraintName, vars, coefs, "L", -cConstant)

    # Indicators
    # z = 1  -->  v >= alpha
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOn_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=1,
        sense="G",
        constant=ALPHA,
    )

    # z = 0  -->  v <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOff_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=0,
        sense="L",
        constant=0,
    )

    return problem


def buildDictionaryMCSGeneTargetedProblem(
    cobraModel: cobra.Model,
    gDict: dict,
    matrixG: scipy.sparse.csc_matrix,
    relationships: dict,
    maxKOLength: int,
    numberNewGenesByKO: dict,
    geneSubset: list,
    targetKOs: list,
    minKOLength: int = 0,
    forceLength: bool = False,
    verbose: int = 0,
    targetB: float = 1e-3,
    cConstant: float = 100,
):
    """
    Builds the dictionary that defines the MILP problem to be solved by the MCS algorithm.

    :param modelStruct: Cobra Model, the metabolic model to be operated on. Only Cobra models are supported.
    :param forceLength: boolean. If true, the length of the solution is forced to be off a given length each iteration, otherwise the length is not forced.
    :param targetB: float. Desired activity level of the metabolic task to be disrupted.
    :param cConstant: float.  Used to activate w variable
    :param rxnToForce: list of strings. Reactions to be forced to be part of the solution. Empty by default.

    If the bounds are left as None they are meant to be infinity (None = Infinity), so each solver can interpret it as it is defined.

    :return: problemDict. Dictionary that defines the MILP problem to be solved by any solver defined in gMCSpy.
    """

    problem = OptimizationProblem()
    problemMethods = MinimalCutSetProblem()

    # Alpha. Used to relate the lower bound of v variables with z variables
    ALPHA = 1

    numberOfPossibleKO = matrixG.shape[0]
    possibleKOs = list(gDict.keys())

    if relationships is not None:
        numberOfRelationships = sum(len(lst) for lst in relationships.values())
    else:
        numberOfRelationships = 0

    matrices = problemMethods.prepareConstraintMatrices(
        cobraModel=cobraModel, companionMatrix=matrixG
    )
    matrixS = matrices["S"]
    matrixG = matrices["Companion"]

    # Find maximum value of matrixS
    maxS = np.max(matrixS.data)
    maxS = 10 ** math.ceil(math.log10(maxS))
    # targetB should 5 orders of manitude smaller than the maximum value of matrixS
    targetB = maxS * 1e-8

    # S matrix changed size
    [metabolites, reactions] = matrixS.shape

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = objectiveVector.index(1)

    # Target vector
    t = np.zeros(reactions)
    t[indexC] = 1

    if geneSubset:
        if verbose > 0:
            print("The solution space is made up by these genes: " + str(geneSubset))

        # Find the interventions that contain the genes in the geneSubset

        interventionsToKeep = []
        for gene in geneSubset:
            for intervention in possibleKOs:
                if gene.issubset(intervention):
                    interventionsToKeep.append(intervention)

        interventionsIndex = [
            possibleKOs.index(intervention) for intervention in interventionsToKeep
        ]
        interventionsIndex = sorted(list(set(interventionsIndex)))
        matrixG = matrixG[interventionsIndex, :]
        n = matrixG.shape[0]

    targetKOIndex = np.zeros(len(possibleKOs))
    for gene in targetKOs:
        for intervention in possibleKOs:
            if gene == intervention:
                targetKOIndex[possibleKOs.index(intervention)] = 1

    vIndex = {}
    vGroups = {}
    vIndex["u"] = list(range(0, metabolites))
    vIndex["vp"] = list(range(vIndex["u"][-1] + 1, vIndex["u"][-1] + n + 1))
    vIndex["w"] = [vIndex["vp"][-1] + 1]
    vIndex["zp"] = list(range(vIndex["w"][-1] + 1, vIndex["w"][-1] + n + 1))
    vIndex["zw"] = [vIndex["zp"][-1] + 1]
    vIndex["epsp"] = list(range(vIndex["zw"][-1] + 1, vIndex["zw"][-1] + n + 1))
    vIndex["espw"] = [vIndex["epsp"][-1] + 1]
    vIndex["delp"] = list(range(vIndex["espw"][-1] + 1, vIndex["espw"][-1] + n + 1))
    vIndex["delw"] = [vIndex["delp"][-1] + 1]
    vIndex["x"] = list(
        range(vIndex["delw"][-1] + 1, vIndex["delw"][-1] + reactions + 2)
    )
    nVars = vIndex["x"][-1] + 1

    for u in vIndex["u"]:
        problem.addVariable(
            index=u, prefix="u_", lb="infinity", ub="infinity", vtype="continuous"
        )

    for vp in vIndex["vp"]:
        problem.addVariable(
            index=vp, prefix="vp_", lb=0, ub="infinity", vtype="continuous"
        )

    for w in vIndex["w"]:
        problem.addVariable(
            index=w, prefix="w_", lb=0, ub="infinity", vtype="continuous"
        )

    for zp in vIndex["zp"]:
        problem.addVariable(index=zp, prefix="zp_", lb=0, ub=1, vtype="binary")

    for zw in vIndex["zw"]:
        problem.addVariable(index=zw, prefix="zw_", lb=0, ub=1, vtype="binary")

    for epsp in vIndex["epsp"]:
        problem.addVariable(
            index=epsp, prefix="epsp_", lb=0, ub="infinity", vtype="continuous"
        )

    for espw in vIndex["espw"]:
        problem.addVariable(index=espw, prefix="espw_", lb=0, ub=0, vtype="continuous")

    for delp in vIndex["delp"]:
        problem.addVariable(
            index=delp, prefix="delp_", lb=0, ub="infinity", vtype="continuous"
        )

    for delw in vIndex["delw"]:
        problem.addVariable(index=delw, prefix="delw_", lb=0, ub=0, vtype="continuous")

    for x in vIndex["x"][:-1]:
        problem.addVariable(
            index=x, prefix="x_", lb=0, ub="infinity", vtype="continuous"
        )
    problem.addVariable(
        index=vIndex["x"][-1], prefix="x_", lb=1000, ub="infinity", vtype="continuous"
    )

    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]
    vGroups["eps"] = vIndex["epsp"] + vIndex["espw"]
    vGroups["del"] = vIndex["delp"] + vIndex["delw"]

    cIndex = {}

    cIndex["Ndual"] = list(range(0, reactions))
    cIndex["forceBioCons"] = [cIndex["Ndual"][-1] + 1]
    cIndex["forceKnockOut"] = [cIndex["forceBioCons"][-1] + 1]

    numberOfLinearCombinations = matrixS.shape[0] + matrixG.shape[0] + 1 + 1
    cIndex["linearCombination"] = list(
        range(
            cIndex["forceKnockOut"][-1] + 1,
            cIndex["forceKnockOut"][-1] + numberOfLinearCombinations,
        )
    )

    cIndex["forceLength"] = [cIndex["linearCombination"][-1] + 1]
    nCons = cIndex["forceLength"][-1] + 1

    A = scipy.sparse.lil_matrix((nCons, nVars))
    A[np.ix_(cIndex["Ndual"], vIndex["u"])] = matrixS.transpose()
    A[np.ix_(cIndex["Ndual"], vIndex["vp"])] = matrixG.transpose()
    A[np.ix_(cIndex["Ndual"], vIndex["w"])] = -t
    A[np.ix_(cIndex["forceBioCons"], vIndex["w"])] = -targetB
    A[np.ix_(cIndex["forceKnockOut"], vIndex["vp"])] = targetKOIndex

    auxS = scipy.sparse.hstack(
        [matrixS, scipy.sparse.lil_matrix(np.zeros((metabolites, 1)))]
    )
    auxK = scipy.sparse.hstack([matrixG, scipy.sparse.lil_matrix(np.zeros((n, 1)))])
    auxT = scipy.sparse.hstack([scipy.sparse.lil_matrix(-t), targetB])
    A[np.ix_(cIndex["linearCombination"], vIndex["x"])] = scipy.sparse.vstack(
        [auxS, auxK, auxT]
    )

    matrixEp1 = scipy.sparse.lil_matrix(
        (metabolites, (len(vIndex["vp"]) + len(vIndex["w"])))
    )
    matrixEp2 = scipy.sparse.eye((len(vIndex["vp"]) + len(vIndex["w"])), format="lil")
    auxEp = scipy.sparse.vstack([matrixEp1, -matrixEp2])
    A[np.ix_(cIndex["linearCombination"], vGroups["eps"])] = auxEp

    A[np.ix_(cIndex["linearCombination"], vGroups["del"])] = -auxEp

    rhs = np.zeros(nCons)
    rhs[cIndex["Ndual"]] = None
    rhs[cIndex["forceBioCons"][0]] = -cConstant
    rhs[cIndex["forceKnockOut"]] = None
    rhs[
        cIndex["linearCombination"][metabolites : len(cIndex["linearCombination"]) - 1]
    ] = targetKOIndex

    lhs = np.zeros(nCons)
    lhs[cIndex["Ndual"]] = 0
    lhs[cIndex["forceBioCons"][0]] = -1
    lhs[cIndex["forceKnockOut"]] = 10
    lhs[
        cIndex["linearCombination"][metabolites : len(cIndex["linearCombination"]) - 1]
    ] = targetKOIndex

    sense = np.array(["E"] * nCons)
    sense[cIndex["Ndual"]] = "G"
    sense[cIndex["forceBioCons"][0]] = "G"
    sense[cIndex["forceKnockOut"]] = "G"

    A = A.tocsr()
    problem.addConstraintsFromSparseMatrix(A, lhs, sense)

    varsForceBioCons = A[cIndex["forceBioCons"][0]].indices
    coefficients = A[cIndex["forceBioCons"][0]].data
    constraintName = f'c{cIndex["forceBioCons"][0]}'
    problem.modifyConstraint(
        constraintName=constraintName,
        variables=varsForceBioCons,
        coefficients=coefficients,
        sense="L",
        rhs=-cConstant,
    )
    # problem.addConstraint(name="forceBioCons_2",variables=varsForceBioCons, coefficients=coefficients, sense="L", rhs=-cConstant)

    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]
    vGroups["eps"] = vIndex["epsp"] + vIndex["espw"]
    vGroups["del"] = vIndex["delp"] + vIndex["delw"]

    # Indicators
    # z = 1  -->  v >= alpha
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOn_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        sense="G",
        constant=ALPHA,
        activeWhen=1,
    )

    # z = 0  -->  v <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOff_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        sense="L",
        constant=0,
        activeWhen=0,
    )

    # z = 1  -->  epsilon <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOnEpsilon_",
        variables=vGroups["eps"],
        indicatorVariables=vGroups["z"],
        sense="L",
        constant=0,
        activeWhen=1,
    )

    # z = 0  -->  epsilon <= M
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOffEpsilon_",
        variables=vGroups["eps"],
        indicatorVariables=vGroups["z"],
        sense="L",
        constant=10,
        activeWhen=0,
    )

    vpVariableNames = {}
    for index, name in zip(vIndex["zp"], gDict.keys()):
        vpVariableNames[name] = index

    if numberOfRelationships > 0:
        subIndex = 0
        for key, values in relationships.items():
            vpContainer = vpVariableNames[key]
            for value in values:
                constraintName = f"relationship_constraint{subIndex}"
                vpContained = vpVariableNames[value]
                vars = [vpContainer, vpContained]
                coefs = [-1, 1]
                problem.addConstraint(constraintName, vars, coefs, "G", 0)
                subIndex += 1

    # Build objective
    objectiveVars = [
        vIndex["zp"][i]
        for i in range(0, len(possibleKOs))
    ]
    problem.addObjective(objectiveVars, list(numberNewGenesByKO.values()), sense="minimize")

    if forceLength:
        constraintName = "forceLength"
        vars = objectiveVars
        coefs = np.ones(len(objectiveVars))
        problem.addConstraint(constraintName, vars, coefs, "L", 1)

    return problem

def buildDictionaryRegNetwork(
    cobraModel: cobra.Model,
    forceLength: bool,
    compressedNodes: List[str],
    rxnToForce: List[str],
    verbose: int = 0,
    targetB: float = 1e-3,
    cConstant: float = 1e-3,
):
    """
    Builds the dictionary that defines the MILP problem to be solved by the MCS algorithm in a Regulatory network.

    :param modelStruct: Cobra Model, the metabolic model to be operated on. Only Cobra models generated by expanded network are supported.
    :param forceLength: boolean. If true, the length of the solution is forced to be off a given length each iteration, otherwise the length is not forced.
    :param targetB: float. Desired activity level of the metabolic task to be disrupted.
    :param cConstant: float.  Used to activate w variable
    :param rxnToForce: list of strings. Reactions to be forced to be part of the solution. Empty by default.

    If the bounds are left as None they are ment to be infinity (None = Infinity), so each solver can interpret it as it is defined.

    :return: problemDict. Dictionary that defines the MILP problem to be solved by any solver defined in gMCSpy.
    """

    problemMethods = MinimalCutSetProblem()
    problem = OptimizationProblem()

    # Alpha. Used to relate the lower bound of v variables with z variables
    ALPHA = 1

    matrices = problemMethods.prepareConstraintMatrices(cobraModel, "K")

    matrixS = matrices["S"]
    matrixK = matrices["Companion"]

    # S matrix changed size
    [metabolites, reactions] = matrixS.shape

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = objectiveVector.index(1)

    # Target vector
    t = np.zeros(reactions)
    t[indexC] = 1

    if rxnToForce:
        if verbose > 0:
            print("Forcing reactions to be part of the solution: " + str(rxnToForce))
        # Find the index of the reactions in the model
        rxnIndex = [
            cobraModel.reactions.index(cobraModel.reactions.get_by_id(rxn))
            for rxn in rxnToForce
        ]
        matrixK = matrixK[rxnIndex, :]
        reactionsK = matrixK.shape[0]

    # Add variables to the problem

    vIndex = {}
    vIndex["u"] = list(range(0, metabolites))
    vIndex["vp"] = list(range(vIndex["u"][-1] + 1, vIndex["u"][-1] + reactionsK + 1))
    vIndex["w"] = [vIndex["vp"][-1] + 1]
    vIndex["zp"] = list(range(vIndex["w"][-1] + 1, vIndex["w"][-1] + reactionsK + 1))
    vIndex["zw"] = [vIndex["zp"][-1] + 1]
    name = "u_"
    for u in vIndex["u"]:
        problem.addVariable(
            index=u, prefix=name, lb="infinity", ub="infinity", vtype="continuous"
        )

    name = "vp_"
    for vp in vIndex["vp"]:
        problem.addVariable(
            index=vp, prefix=name, lb=0, ub="infinity", vtype="continuous"
        )

    name = "w_"
    w = vIndex["w"][0]
    problem.addVariable(index=w, prefix=name, lb=0, ub="infinity", vtype="continuous")

    count = 0
    for zp in vIndex["zp"]:
        name = f"zp_{cobraModel.reactions[rxnIndex[count]].id}"
        problem.addVariable(
            index=zp, prefix=name, lb=0, ub=1, vtype="binary", relatedTo=count
        )
        count += 1

    name = "zw_"
    zw = vIndex["zw"][0]
    problem.addVariable(index=zw, prefix=name, lb=0, ub=1, vtype="binary")

    vGroups = {}
    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]

    cIndex = {}

    cIndex["Ndual"] = list(range(0, reactions))
    cIndex["forceBioCons"] = [cIndex["Ndual"][-1] + 1]
    cIndex["forceLength"] = [cIndex["forceBioCons"][-1] + 1]

    matrixS = matrixS.transpose().tocsr()
    matrixK = matrixK.transpose().tocsr()
    tIndexConstraint = np.where(t != 0)[0]
    for i in cIndex["Ndual"]:
        constraintName = f"Ndual_constraint{i}"
        variableIndicesU = matrixS[i].indices
        vars = []
        vars = [vIndex["u"][j] for j in variableIndicesU]
        variableIndicesVp = matrixK[i].indices
        vars.extend([vIndex["vp"][j] for j in variableIndicesVp])
        coefs = np.concatenate((matrixS[i].data, matrixK[i].data))
        problem.addConstraint(
            name=constraintName, variables=vars, coefficients=coefs, sense="G", rhs=0
        )

    constraintName = f"Ndual_constraint{tIndexConstraint[0]}"
    constraint = problem.getConstraint(constraintName)
    constraint["variables"].extend([w])
    coefs = np.concatenate((constraint["coefficients"], -t[tIndexConstraint]))
    problem.modifyConstraint(
        constraintName=constraintName,
        variables=constraint["variables"],
        coefficients=coefs,
        sense=constraint["sense"],
        rhs=constraint["rhs"],
    )

    problem.addConstraint(
        name="forceBioCons",
        variables=[w],
        coefficients=[-targetB],
        sense="L",
        rhs=-cConstant,
    )

    # obatin de indices of the zp_on 
    rxn_on = [rxnToForce.index(i) for i in rxnToForce if 'on' in i] 
    zp_on = [vIndex["zp"][i] for i in rxn_on]
    
    # Force at least one reaction
    problem.addConstraint(
        name="forceAtLeastOne",
        variables=zp_on,
        coefficients=np.ones(len(zp_on)),
        sense="G",
        rhs=1,
    )


    if forceLength:
        variables = zp_on
        problem.addConstraint(
            name="forceLength",
            variables=variables,
            coefficients=np.ones(len(variables)),
            sense="E",
            rhs=1,
        )
    
    # Relationships contraints 
    for node in compressedNodes:
        # Obtain the index of the reactions that are related to the node
        rxn_index = [rxnToForce.index(i) for i in rxnToForce if node in i]
        zp_index = [vIndex["zp"][i] for i in rxn_index]
        problem.addConstraint(
            name=f"relationship_{node}",
            variables=zp_index,
            coefficients=np.ones(len(zp_index)),
            sense="E",
            rhs=1,
        )
        
        
    # Build objective
    problem.addObjective(zp_on, np.ones(len(zp_on)), sense="minimize")

    # Indicators
    # z = 1  -->  v >= alpha
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOn_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=1,
        sense="G",
        constant=ALPHA,
    )

    # z = 0  -->  v <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOff_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=0,
        sense="L",
        constant=0,
    )

    return problem

def buildDictionaryDossageNetwork(
    cobraModel: cobra.Model,
    forceLength: bool,
    compressedNodes: List[str],
    rxnToForce: List[str],
    verbose: int = 0,
    targetB: float = 1e-3,
    cConstant: float = 1e-3,
):
    """
    Builds the dictionary that defines the MILP problem to be solved by the MCS algorithm in a Regulatory network.

    :param modelStruct: Cobra Model, the metabolic model to be operated on. Only Cobra models generated by expanded network are supported.
    :param forceLength: boolean. If true, the length of the solution is forced to be off a given length each iteration, otherwise the length is not forced.
    :param targetB: float. Desired activity level of the metabolic task to be disrupted.
    :param cConstant: float.  Used to activate w variable
    :param rxnToForce: list of strings. Reactions to be forced to be part of the solution. Empty by default.

    If the bounds are left as None they are ment to be infinity (None = Infinity), so each solver can interpret it as it is defined.

    :return: problemDict. Dictionary that defines the MILP problem to be solved by any solver defined in gMCSpy.
    """

    problemMethods = MinimalCutSetProblem()
    problem = OptimizationProblem()

    # Alpha. Used to relate the lower bound of v variables with z variables
    ALPHA = 1

    matrices = problemMethods.prepareConstraintMatrices(cobraModel, "K")

    matrixS = matrices["S"]
    matrixK = matrices["Companion"]

    # S matrix changed size
    [metabolites, reactions] = matrixS.shape

    # Objective vector of model
    objectiveVector = [r.objective_coefficient for r in cobraModel.reactions]
    indexC = objectiveVector.index(1)

    # Target vector
    t = np.zeros(reactions)
    t[indexC] = 1

    if rxnToForce:
        if verbose > 0:
            print("Forcing reactions to be part of the solution: " + str(rxnToForce))
        # Find the index of the reactions in the model
        rxnIndex = [
            cobraModel.reactions.index(cobraModel.reactions.get_by_id(rxn))
            for rxn in rxnToForce
        ]
        matrixK = matrixK[rxnIndex, :]
        reactionsK = matrixK.shape[0]

    # Add variables to the problem

    vIndex = {}
    vIndex["u"] = list(range(0, metabolites))
    vIndex["vp"] = list(range(vIndex["u"][-1] + 1, vIndex["u"][-1] + reactionsK + 1))
    vIndex["w"] = [vIndex["vp"][-1] + 1]
    vIndex["zp"] = list(range(vIndex["w"][-1] + 1, vIndex["w"][-1] + reactionsK + 1))
    vIndex["zw"] = [vIndex["zp"][-1] + 1]
    name = "u_"
    for u in vIndex["u"]:
        problem.addVariable(
            index=u, prefix=name, lb="infinity", ub="infinity", vtype="continuous"
        )

    name = "vp_"
    for vp in vIndex["vp"]:
        problem.addVariable(
            index=vp, prefix=name, lb=0, ub="infinity", vtype="continuous"
        )

    name = "w_"
    w = vIndex["w"][0]
    problem.addVariable(index=w, prefix=name, lb=0, ub="infinity", vtype="continuous")

    count = 0
    for zp in vIndex["zp"]:
        name = f"zp_{cobraModel.reactions[rxnIndex[count]].id}"
        problem.addVariable(
            index=zp, prefix=name, lb=0, ub=1, vtype="binary", relatedTo=count
        )
        count += 1

    name = "zw_"
    zw = vIndex["zw"][0]
    problem.addVariable(index=zw, prefix=name, lb=0, ub=1, vtype="binary")

    vGroups = {}
    vGroups["v"] = vIndex["vp"] + vIndex["w"]
    vGroups["z"] = vIndex["zp"] + vIndex["zw"]

    cIndex = {}

    cIndex["Ndual"] = list(range(0, reactions))
    cIndex["forceBioCons"] = [cIndex["Ndual"][-1] + 1]
    cIndex["forceLength"] = [cIndex["forceBioCons"][-1] + 1]

    matrixS = matrixS.transpose().tocsr()
    matrixK = matrixK.transpose().tocsr()
    tIndexConstraint = np.where(t != 0)[0]
    for i in cIndex["Ndual"]:
        constraintName = f"Ndual_constraint{i}"
        variableIndicesU = matrixS[i].indices
        vars = []
        vars = [vIndex["u"][j] for j in variableIndicesU]
        variableIndicesVp = matrixK[i].indices
        vars.extend([vIndex["vp"][j] for j in variableIndicesVp])
        coefs = np.concatenate((matrixS[i].data, matrixK[i].data))
        problem.addConstraint(
            name=constraintName, variables=vars, coefficients=coefs, sense="G", rhs=0
        )

    constraintName = f"Ndual_constraint{tIndexConstraint[0]}"
    constraint = problem.getConstraint(constraintName)
    constraint["variables"].extend([w])
    coefs = np.concatenate((constraint["coefficients"], -t[tIndexConstraint]))
    problem.modifyConstraint(
        constraintName=constraintName,
        variables=constraint["variables"],
        coefficients=coefs,
        sense=constraint["sense"],
        rhs=constraint["rhs"],
    )

    problem.addConstraint(
        name="forceBioCons",
        variables=[w],
        coefficients=[-targetB],
        sense="L",
        rhs=-cConstant,
    )


    # Relationships contraints (11) (12) (13) 
    subFixes_KO = ['KO_on_redflux', 'KO_off_redflux']
    subFixes_KI = ['KI_on_redflux', 'KI_off_redflux']
    subFixes_KO_KI = ['KO_on_redflux', 'KI_off_redflux']
    objectiveNodes = []
    for node in compressedNodes:
            
        # Obtain the index of the reactions that are related to the node and are KO (11)
        rxn_index = [rxnToForce.index(f'{node}_{i}') for i in subFixes_KO]
        zp_index = [vIndex["zp"][i] for i in rxn_index]
        problem.addConstraint(
            name=f"relationship_{node}_KO",
            variables=zp_index,
            coefficients=np.ones(len(zp_index)),
            sense="L",
            rhs=1,
        )
        # Obtain the index of the reactions that are related to the node and are KI (12)
        rxn_index = [rxnToForce.index(f'{node}_{i}') for i in subFixes_KI]
        zp_index = [vIndex["zp"][i] for i in rxn_index]
        problem.addConstraint(
            name=f"relationship_{node}_KI",
            variables=zp_index,
            coefficients=np.ones(len(zp_index)),
            sense="L",
            rhs=1,
        )
        # Obtain the index of the reactions that are related to the node and are KI_off and KO_on (13)
        rxn_index = [rxnToForce.index(f'{node}_{i}') for i in subFixes_KO_KI]
        zp_index = [vIndex["zp"][i] for i in rxn_index]
        problem.addConstraint(
            name=f"relationship_{node}_KI_KO",
            variables=zp_index,
            coefficients=np.ones(len(zp_index)),
            sense="L",
            rhs=1,
        )
        objectiveNodes.extend(zp_index)
        
    # At least one active node constraint (14)
    problem.addConstraint(
        name='SumNodes',
        variables=objectiveNodes,
        coefficients=np.ones(len(objectiveNodes)),
        sense='G',
        rhs=1
    )       

    if forceLength:
        problem.addConstraint(
            name="forceLength",
            variables=objectiveNodes,
            coefficients=np.ones(len(objectiveNodes)),
            sense="E",
            rhs=1,
        )
     
        
    # Build objective
    problem.addObjective(objectiveNodes, np.ones(len(objectiveNodes)), sense="minimize")

    # Indicators
    # z = 1  -->  v >= alpha
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOn_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=1,
        sense="G",
        constant=ALPHA,
    )

    # z = 0  -->  v <= 0
    problemMethods.buildIndicatorsIntoProblem(
        problem=problem,
        prefix="indicatorZOff_",
        variables=vGroups["v"],
        indicatorVariables=vGroups["z"],
        activeWhen=0,
        sense="L",
        constant=0,
    )

    return problem

def readjustmentProblem(rules, gMIs):
    opt = OptimizationProblem()
    variables = []
    constraints = {}
    for id, rule in rules.items():
        x, eq = rule['rule'].split(' = ')
        if x not in variables:
            variables.append(x)
            variableIndex = len(variables) - 1
        else:
            variableIndex = variables.index(x)
        
        if '&' in eq:
            count = 0
            eq = eq.split(' & ')
            rightTerms = []
            
            coefficients = []
            negatedTerms = 0
            
            for i in eq:
                
                negated = False
                if i.startswith('!'):
                    i = i[1:]
                    negatedTerms += 1
                    negated = True
                    
                if i not in variables:
                    variables.append(i)
                    rightTerms.append(len(variables) - 1)
                else:
                    rightTerms.append(variables.index(i))
                idString = f'{id}_{count}'
                
                if negated:
                    coef = 1
                    rhs = 1
                else:
                    coef = -1
                    rhs = 0
                
                coefficients.append(coef)
                
                constraints[idString] = {
                    'variables': [variableIndex] + [rightTerms[count]],
                    'coefficients': [1] + [coef],
                    'sense': 'L',
                    'rhs': rhs
                }
                count += 1
            constraints[id] = {
                'variables': [variableIndex] + rightTerms,
                'coefficients': [1] + coefficients,
                'sense': 'G',
                'rhs': -len(rightTerms) + 1 + negatedTerms
            }
        elif '|' in eq:
            eq = eq.split(' | ')
            rightTerms = []
            
            coefficients = []
            negatedTerms = 0
            
            count = 0
            for i in eq:
                negated = False
                if i.startswith('!'):
                    i = i[1:]
                    negatedTerms += 1
                    negated = True
                    
                if i not in variables:
                    variables.append(i)
                    rightTerms.append(len(variables) - 1)
                else:
                    rightTerms.append(variables.index(i))
                idString = f'{id}_{count}'
                
                if negated:
                    coef = 1
                    rhs = 1
                else:
                    coef = -1
                    rhs = 0

                coefficients.append(coef)
                
                constraints[idString] = {
                    'variables': [variableIndex] + [rightTerms[count]],
                    'coefficients': [1] + [coef],
                    'sense': 'G',
                    'rhs': rhs
                }
                count += 1
                    
            constraints[id] = {
                'variables': [variableIndex] + rightTerms,
                'coefficients': [1] + coefficients,
                'sense': 'L',
                'rhs': negatedTerms
            }
        elif eq.startswith('!'):
            eq = eq[1:]
            if eq not in variables:
                variables.append(eq)
                rightTerm = len(variables) - 1
            else:
                rightTerm = variables.index(eq)
            constraints[id] = {
                'variables': [variableIndex] + [rightTerm],
                'coefficients': [1, 1],
                'sense': 'E',
                'rhs': 1
            }
    
    for i, variable in enumerate(variables):
        opt.addVariable(index=i, prefix=f'{variable}--', lb=0, ub=1, vtype='binary')
        
    for id, constraint in constraints.items():
        opt.addConstraint(name=id, **constraint)
        
    objIndex = []
    objCoefficients = []
    for intervention in gMIs:
        id, coef = intervention.split('_K')
        if id in variables:
            vIndex = variables.index(id)
            objIndex.append(vIndex)
            if coef == 'I':
                objCoefficients.append(-1)
            elif coef == 'O':
                objCoefficients.append(1)
    
    opt.addObjective(objIndex, objCoefficients, sense='minimize')
    
    return opt
        
        
    
    