import numpy as np
from .Utilities import createLog
from .Utilities import beartype


def cplexProblemInterpreter(
    problemDict: dict, optModel, interface, timelimit=1e5, verbose=0
):
    """
    The purpose of this function is to interpret the problem dictionary into cplex, so we can solve it.
    Any problem with the correct format can be solved with this function.

    :param problemDict: Dictionary with the problem information
    :param interface: Interface to use to solve the problem, this is automatically set by the function setSolver.
        See the :func:`gMCSpy.Utilities.setSolver` function for more information.

    :returns:
        - **problemInterpreted:** An cplex model that contains the solver and the problem to solve.
    """

    # Setting cplex parameters
    logName = f"{problemDict.logPath}/solverLog/cplex_{problemDict.getModelName()}_{problemDict.getCharacteristicParameter()}.log"
    if verbose == 2:
        filename = createLog(logName)
        cplexlog = open(filename, "w")
        optModel.set_log_stream(cplexlog)
        optModel.set_error_stream(cplexlog)
        optModel.set_warning_stream(cplexlog)
        optModel.set_results_stream(cplexlog)
    elif verbose == 0:
        optModel.set_log_stream(None)
        optModel.set_error_stream(None)
        optModel.set_warning_stream(None)
        optModel.set_results_stream(None)

    # Transform the problem bounds to cplex readables
    infinity = interface.infinity

    variables = problemDict.getVariables()

    for key, var in variables.items():
        optModel.variables.add(names=[var["name"]])
        if var["lb"] == "infinity":
            var["lb"] = -infinity
        if var["ub"] == "infinity":
            var["ub"] = infinity
        optModel.variables.set_lower_bounds(var["name"], var["lb"])
        optModel.variables.set_upper_bounds(var["name"], var["ub"])

        if var["vtype"] == "binary":
            optModel.variables.set_types(var["name"], optModel.variables.type.binary)

        elif var["vtype"] == "continuous":
            optModel.variables.set_types(
                var["name"], optModel.variables.type.continuous
            )

    if verbose > 0:
        print("Variables added")

    for key, data in problemDict.getConstraints().items():
        variableNamesCon = [variables[j]["name"] for j in data["variables"]]
        coefficients = data["coefficients"]
        sense = data["sense"]
        optModel.linear_constraints.add(
            lin_expr=[interface.SparsePair(ind=variableNamesCon, val=coefficients)],
            rhs=[data["rhs"]],
            senses=[sense],
            names=[key],
        )

    if verbose > 0:
        print("Constraints added")

    if problemDict.getIndicators():
        for key, data in problemDict.getIndicators().items():
            indicatorVariable = variables[data["indicatorVariable"]]["name"]
            variable = [variables[data["variable"]]["name"]]

            if data["activeWhen"] == 1:
                complemented = 0
                senseCon = data["sense"]

            else:
                complemented = 1
                senseCon = data["sense"]

            expr = interface.SparsePair(ind=variable, val=np.ones(len(variable)))
            optModel.indicator_constraints.add(
                indvar=indicatorVariable,
                complemented=complemented,
                lin_expr=expr,
                rhs=data["b"],
                sense=senseCon,
                name=key,
            )
        if verbose > 0:
            print("Indicators added")

    # Define objective function of the model
    if problemDict.getSense() == "maximize":
        optModel.objective.set_sense(optModel.objective.sense.maximize)
    elif problemDict.getSense() == "minimize":
        optModel.objective.set_sense(optModel.objective.sense.minimize)
    else:
        raise ValueError("Objective sense not recognized")

    objective = problemDict.getObjective()

    for i, coeff in enumerate(objective["coefficients"]):
        optModel.objective.set_linear(objective["variables"][i], int(coeff))

    if verbose > 0:
        print("Objective Function added")

    return optModel


def gurobiProblemInterpreter(
    problemDict: dict, env, interface, timelimit=1e5, verbose=0
):
    """
    The purpose of this function is to interpret the problem dictionary into gurobi, so we can solve it.
    Any problem with the correct format can be solved with this function.

    :param problemDict: Dictionary with the problem information
    :param interface: Interface to use to solve the problem, this is automatically set by the function setSolver.
        See the :func:`gMCSpy.Utilities.setSolver` function for more information.

    :returns:
        - **problemInterpreted:** An gurobi model that contains the solver and the problem to solve.
    """


    logName = f"{problemDict.logPath}/solverLog/gurobi_{problemDict.getModelName()}_{problemDict.getCharacteristicParameter()}.log"

    # verbose control
    if verbose == 2:
        env.setParam("OutputFlag", 1)
        env.setParam("LogToConsole", 0)
        filename = createLog(logName)
        env.setParam("LogFile", filename)
    elif verbose == 0:
        env.setParam("OutputFlag", 0)
    else:
        env.setParam("OutputFlag", 1)

    env.start()

    optModel = interface.Model(name="MCS", env=env)

    # Transform the problem bounds to gurobi readables
    infinity = interface.GRB.INFINITY

    # change "infinity" to value of infinity

    for key, var in problemDict.getVariables().items():
        if var["vtype"] == "binary":
            optModel.addVar(
                name=var["name"],
                vtype=interface.GRB.BINARY,
                lb=var["lb"],
                ub=var["ub"],
            )
        elif var["vtype"] == "continuous":
            if var["lb"] == "infinity":
                var["lb"] = -infinity
            if var["ub"] == "infinity":
                var["ub"] = infinity
            optModel.addVar(
                name=var["name"],
                vtype=interface.GRB.CONTINUOUS,
                lb=var["lb"],
                ub=var["ub"],
            )

    optModel.update()
    if verbose > 0:
        print("Variables added")

    for key, data in problemDict.getConstraints().items():
        variableNames = [optModel.getVars()[j] for j in data["variables"]]
        expr = interface.LinExpr(data["coefficients"], variableNames)
        if data["sense"] == "E":
            senseCon = interface.GRB.EQUAL
        elif data["sense"] == "G":
            senseCon = interface.GRB.GREATER_EQUAL
        elif data["sense"] == "L":
            senseCon = interface.GRB.LESS_EQUAL

        optModel.addLConstr(expr, senseCon, data["rhs"], name=key)

    optModel.update()
    if verbose > 0:
        print("Constraints added")

    if problemDict.getIndicators():
        for key, data in problemDict.getIndicators().items():
            indicatorVariable = optModel.getVars()[data["indicatorVariable"]]
            variable = optModel.getVars()[data["variable"]]
            exprCon = interface.LinExpr(variable)

            if data["sense"] == "G":
                senseCon = interface.GRB.GREATER_EQUAL
            else:
                senseCon = interface.GRB.LESS_EQUAL

            optModel.addGenConstrIndicator(
                indicatorVariable, data["activeWhen"], exprCon, senseCon, data["b"], key
            )

        if verbose > 0:
            print("Indicators added")

    optModel.update()

    # Define objective function of the model
    if problemDict.getSense() == "maximize":
        optModel.ModelSense = interface.GRB.MAXIMIZE
    elif problemDict.getSense() == "minimize":
        optModel.ModelSense = interface.GRB.MINIMIZE
    else:
        raise ValueError("Objective sense not recognized")

    objective = problemDict.getObjective()
    objectiveVars = [optModel.getVars()[i] for i in objective["variables"]]
    objtiveExpr = interface.LinExpr(objective["coefficients"], objectiveVars)

    optModel.setObjective(objtiveExpr)
    optModel.update()
    if verbose > 0:
        print("Objective Function added")
   
    return optModel


def scipProblemInterpreter(
    problemDict: dict,
    interface,
    timelimit=1e5,
    verbose=0,
):
    model = interface.Model("MCS")
    modelVars = {}
    for key, var in problemDict.getVariables().items():
        if var["vtype"] == "binary":
            modelVars[key] = model.addVar(
                name=var["name"],
                vtype="B",
                lb=var["lb"],
                ub=var["ub"],
            )
        elif var["vtype"] == "continuous":
            if var["lb"] == "infinity" and var["ub"] == "infinity":
                lb = None
                ub = None
            elif var["ub"] == "infinity":
                ub = None
                lb = var["lb"]
            elif var["lb"] == "infinity":
                lb = None
                ub = var["ub"]
            else:
                ub = var["ub"]
                lb = var["lb"]
            modelVars[key] = model.addVar(
                name=var["name"],
                vtype="C",
                lb=lb,
                ub=ub,
            )

    if verbose > 0:
            print("Variables added")
    
    for key, constraint in problemDict.getConstraints().items():
        if len(constraint['variables']) == 0:
            Warning(f"Constraint {key} has no variables")
            continue
        
        if len(constraint['coefficients']) == 0:
            Warning(f"Constraint {key} has no coefficients")
            continue
        
        if constraint["sense"] == "E":
            sense = "=="
        elif constraint["sense"] == "G":
            sense = ">="
        elif constraint["sense"] == "L":
            sense = "<="
        constraint_string = f"{' + '.join([f'{coeff}*modelVars[{var}]' for var, coeff in zip(constraint['variables'], constraint['coefficients'])])} {sense} {constraint['rhs']}"
        code = f"model.addCons({constraint_string}, modifiable=True, name='{key}')"
        exec(code)
    
    if verbose > 0:                  
        print("Constraints added")
        
    if len(problemDict.getIndicators()) > 0:
        for key, indicator in problemDict.getIndicators().items():
            indicatorVariable = indicator["indicatorVariable"]
            if indicator["sense"] == "G":
                sense = ">="
            elif indicator["sense"] == "L":
                sense = "<="
            elif indicator["sense"] == "E":
                sense = "=="
            if indicator["activeWhen"] == 1:
                active = True
            else:
                active = False
            constraint_string = f"{' + '.join([f'modelVars[{var}]' for var in [indicator['variable']]])} {sense} {indicator['b']}"
            
            code = f"model.addConsIndicator({constraint_string}, binvar=modelVars[{indicatorVariable}], activeone={active}, name='{key}')"
            exec(code)
    
    if verbose > 0:        
        print("Indicators added")
    objective = problemDict.getObjective()
    
    constraint_string = f"{' + '.join([f'{coeff}*modelVars[{var}]' for var, coeff in zip(objective['variables'], objective['coefficients'])])}"
    code = f"model.setObjective({constraint_string},'{problemDict.getSense()}')"
    exec(code)
    
    if verbose > 0:
        print("Objective added")
    
    model.setHeuristics(3)
    model.hideOutput()
    #model.setBoolParam("constraints/countsols/collect", True)
    #model.setPresolve(1)
    return model