from .Utilities import setSolver
from .Utilities import pd
from .ProblemInterpreters import cplexProblemInterpreter
from .ProblemInterpreters import gurobiProblemInterpreter
from .ProblemDefinitions import buildDictionaryMCSProblem
from .ProblemDefinitions import buildDictionaryMCSTargetedProblem
from .ProblemDefinitions import cobra
from .ProblemInterpreters import np
from .Utilities import beartype
from .ProblemDefinitions import List
from tqdm import tqdm

from .MinimalCutSetClass import MinimalCutSetProblem


def calculateMCS(
    cobraModel: cobra.Model,
    MAX_MCS: int,
    MAX_LENGTH_MCS: int,
    solver_params: dict = None,
    rxnSubset: List[str] = None,
    rxnKnockOut: List[str] = None,
    targetB=1e-3,
    timelimit=60,
    forceLength=True,
    numWorkers=0,
    verbose=0,
    useIndicators=True,
    solver: str = "gurobi",
    **kwargs,
):
    """Calculates the Minimal Cut Sets (MCS) of a metabolic network.
    Calculates the MCS of a metabolic network using a mixed integer linear programming (MILP) formulation.
    Among all reactions or with a given subset of reactions, the function finds the minimal set of reactions
    to stop biomass production.

    Required arguments:
        - **cobraModel:** A metabolic model in cobrapy format. See `cobra models <https://cobrapy.readthedocs.io/en/latest/io.html#Reading- and- Writing- Models>`_.
        - **MAX_MCS:** Maximum number of MCS to be found.
        - **MAX_LENGTH_MCS:** Maximum length of MCS to be found. (i.e. number of reactions in the MCS)

    Optional arguments:
        - **rxnToForce:** List of reactions to be considered in the MCS problem. The list is expected to be of reaction names. Default value is an empty list, which means that all reactions will be considered. e.g. rxnToForce = ["R1", "R2", "R3"] will only consider reactions R1, R2 and R3 as possible intervention points. Only ["R1"], ["R2"], ["R3"], ["R1", "R2"], ["R1", "R3"], ["R2", "R3"] and ["R1", "R2", "R3"] could be solutions.

        - **rxnToKnockOut:** List of reactions that you want to be part of the MCS solutions. The list is expected to be of reaction names. Default value is an empty list, which means that no reactions will be forced to be part of the MCS solutions. e.g. rxnToKnockOut = ["R1", "R2", "R3"] will force the MCS solutions to contain at leat one reaction from the list.
        
        - **targetB:** Desired activity level of the metabolic task to be disrupted. Default value is 1e- 3.

        - **timelimit:** Maximum time allowed for the solver to find a solution. Default value is 1e75.

        - **forceLength:** If True, the solver will find MCS of length 1, 2, 3, ..., MAX_LENGTH_MCS. If False, the solver will find MCS of length MAX_LENGTH_MCS. Default value is True.

        - **numWorkers:** Number of workers to be used for parallel processing. Default value is 0. The solvers cplex and gurobi already use parallel processing.

        - **verbose:** If 0, the solver will not print any information. Verbose  = 1, log will be printed into the console. Verbose = 2, log will be printed into the console and a log file will be created. Default value is 0.

        - **useIndicators:** If True, the solver will use binary variables to indicate if a reaction is part of the MCS.

        - **solver:** Solver to be used. Current compatibility includes cplex, gurobi, pulp and optlang_gurobi, optlang_cplex. Default value is cplex.

    """
    # Unpack the kwargs
    if kwargs:
        for key, value in kwargs.items():
            if key == "solver":
                solver = value
            elif key == "useIndicators":
                useIndicators = value
            elif key == "verbose":
                verbose = value
            elif key == "numWorkers":
                numWorkers = value
            elif key == "forceLength":
                forceLength = value
            elif key == "timelimit":
                timelimit = value
            elif key == "targetB":
                targetB = value
            elif key == "MAX_LENGTH_MCS":
                MAX_LENGTH_MCS = value
            elif key == "MAX_MCS":
                MAX_MCS = value
            else:
                raise ValueError("Unknown argument: " + key)

    gurobi_default_params = {
        "MIPFocus": 3,
        "TimeLimit": max(1, timelimit),
        "Threads": numWorkers,
        "PoolSolutions": 100,
        "PoolGap": 0.001,
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

    if solver_params is not None:
        for key, value in solver_params.items():
            if key in params:
                params[key] = value

    # Set the solver based on the user input, current compatibility includes cplex, gurobi, pulp and optlang_gurobi, optlang_cplex
    [interface, useIndicators] = setSolver(solver, useIndicators)

    # Constants
    integralityTolerance = 1e-5
    cConstant = 1e-3  # Used to activate w variable
    infinity = 1e30

    # Initialize the reactions to be knocked out, e.g. reactions that you want to be part of the MCS
    if rxnKnockOut is None:
        rxnKnockOut = []

    # Initialize the reactions to be forced, e.g. reactions that you want to force the flux through.
    if rxnSubset is None:
        rxnSubset = []

    # If there are reactions to be knocked out, then the problem is a targeted problem
    if rxnKnockOut:
        # Check that rxnToKnockOut is in the model and if rxnToForce is not empty check that rxnToKnockOut is in rxnToForce
        if not set(rxnSubset).issubset(set(cobraModel.reactions.list_attr("id"))):
            raise ValueError("One or more reactions in rxnToForce are not in the model")

        # merge rxnToForce and rxnKnockOut, or rxnToKnockOut with all reactions in the model
        if not rxnSubset:
            rxnSubset = cobraModel.reactions.list_attr("id")
        else:
            if not set(rxnKnockOut).issubset(set(rxnSubset)):
                raise ValueError(
                    "One or more reactions in rxnKnockOut are not in rxnToForce"
                )
        # Build the problem object
        problem = buildDictionaryMCSTargetedProblem(
            cobraModel=cobraModel,
            forceLength=forceLength,
            rxnSubset=rxnSubset,
            rxnKnockOut=rxnKnockOut,
        )
        # Set the solver
        problem.setSolver(solver)
        # Interpret the problem, depending on the solver each problem has to be interpreted differently into the solver's interface
        optimization = problem.interpretProblem(params)

        # Import the common problem methods
        problemMethods = MinimalCutSetProblem()

        # Get the objective variables, they are the variables that we are interested in
        obj = problem.getObjective()

        # Set the bounds for the solver, e.g. maximum number of solutions, maximum numbers in a solution, or if we want to force the step by step calculation etc.
        bounds = {
            "MAX_SOLUTION_LENGTH": MAX_LENGTH_MCS,
            "MAX_SOLUTIONS": MAX_MCS,
            "forceLength": False,
        }
        reactions = rxnSubset
        obj = obj["variables"]
        solutionVariables = {"indices": obj, "names": reactions}

        # Solve the problem iteratively
        solutionDict = problemMethods.iterateOverSolver(
            problem=problem,
            problemInterpreted=optimization,
            solver=solver,
            solutionVariables=solutionVariables,
            bounds=bounds,
        )      

        return solutionDict

    else:
        # If there are no specific reactions to be forced through all reactions in the model are included.
        if not rxnSubset:
            rxnSubset = cobraModel.reactions.list_attr("id")

        problem = buildDictionaryMCSProblem(
            cobraModel=cobraModel, forceLength=forceLength, rxnToForce=rxnSubset
        )

        problem.setSolver(solver)
        # Add the solver specific parameters

        optimization = problem.interpretProblem(params)
        problemMethods = MinimalCutSetProblem()
        obj = problem.getObjective()["variables"]
        reactions = rxnSubset
        solutionVariables = {"indices": obj, "names": reactions}
        bounds = {
            "MAX_SOLUTION_LENGTH": MAX_LENGTH_MCS,
            "MAX_SOLUTIONS": MAX_MCS,
            "forceLength": forceLength,
        }
        solutionDict = problemMethods.iterateOverSolver(
            problem=problem,
            problemInterpreted=optimization,
            solver=solver,
            solutionVariables=solutionVariables,
            bounds=bounds,
        )

        return solutionDict
