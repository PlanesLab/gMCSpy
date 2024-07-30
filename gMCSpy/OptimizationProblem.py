from .Utilities import setSolver
from .ProblemInterpreters import cplexProblemInterpreter
from .ProblemInterpreters import gurobiProblemInterpreter
from .ProblemInterpreters import scipProblemInterpreter


"""
The `OptimizationProblem` class represents an optimization problem. It allows users to define variables, constraints, indicators, and an objective function for the problem. It also provides methods to set the solver, interpret the problem, and modify the problem. The class encapsulates the problem data and provides convenient access to its components.

Example Usage:
    # Create an instance of the OptimizationProblem class
    problem = OptimizationProblem()

    # Add variables to the problem
    problem.addVariable(0, 0, 'infinity', 'continuous')
    problem.addVariable(1, 0, 'infinity', 'continuous')
    problem.addVariable(2, 0, 'infinity', 'continuous')

    # Add constraints to the problem
    problem.addConstraint('c0',variables=[0,1,2], coefficients=[1,1,1], sense='L', rhs=100)
    problem.addConstraint('c1',variables=[0,1,2], coefficients=[10,4,5], sense='L', rhs=600)
    problem.addConstraint('c2',variables=[0,1,2], coefficients=[2,2,6], sense='L', rhs=300)

    # Add an objective function to the problem
    problem.addObjective(variables=[0,1,2], coefficients=[10,6,4], sense='maximize')

    # Set the solver for the problem
    problem.setSolver("cplex")

    # Interpret the problem and get the solver interface
    [problemInterpreted, _] = problem.interpretProblem()

    # Solve the problem using the solver interface
    problemInterpreted.solve()

    # Get the optimal solution
    optimalSolution = problemInterpreted.solution.get_objective_value()
    # Print the optimal solution
    print(optimalSolution)
    
Main functionalities:
- Define variables, constraints, indicators, and an objective function for an optimization problem
- Set the solver for the problem
- Interpret the problem and get the solver interface
- Modify the problem by adding, modifying, or slicing variables and constraints
- Get information about the problem, such as variables, constraints, objective function, and solver

Methods:
- addVariable(index, lb, ub, vtype="continuous", prefix="x_", relatedTo=None): Adds a variable to the problem
- addConstraint(name, variables, coefficients, sense, rhs): Adds a constraint to the problem
- addIndicator(name, indicatorVariable, activeWhen, variable, coefficients, sense, b): Adds an indicator constraint to the problem
- addObjective(variables, coefficients, sense): Sets the objective function of the problem
- setSolver(solver): Sets the solver for the problem
- addConstraintsFromSparseMatrix(matrix, rhs, sense): Adds constraints to the problem from a sparse matrix representation
- sliceVariables(variables): Returns a subset of variables from the problem
- modifyConstraint(constraintName, variables, coefficients, sense, rhs): Modifies an existing constraint in the problem
- getProblem(): Returns a dictionary representation of the problem
- getVariables(): Returns a dictionary of variables in the problem
- getConstraints(): Returns a dictionary of constraints in the problem
- getIndicators(): Returns a dictionary of indicators in the problem
- getObjective(): Returns the objective function of the problem
- getSense(): Returns the sense of the objective function
- getSolver(): Returns the solver used for the problem
- getVariable(index): Returns the variable with the specified index
- getConstraint(name): Returns the constraint with the specified name
- getIndicator(name): Returns the indicator constraint with the specified name
- getVariableNames(variablesIndex=None): Returns a dictionary of variable names and indices in the problem
- interpretProblem(parameters={}, verbose=0): Interprets the problem and returns the solver interface and the solver interface object
- setModelName(name): Sets the name of the model
- getModelName(): Returns the name of the model
- setCharacteristicParameter(parameter): Sets the characteristic parameter of the model
- getCharacteristicParameter(): Returns the characteristic parameter of the model
- getRelatedVariable(index): Returns the related variable of the variable with the specified index
"""
class OptimizationProblem:
    def __init__(self):
        self.variables = {}
        self.constraints = {}
        self.indicators = {}
        self.sense = None
        self.objective = None
        self.solver = None
        self.interface = None 
        self.modelName = None
        self.characteristicParameter = None
        self.logPath = None
        self.parameters = None

    def addVariable(self, index, lb, ub, vtype="continuous", prefix="x_", relatedTo=None):
        prefix = f"{prefix}{index}"
        self.variables[index] = {"name": prefix, "lb": lb, "ub": ub, "vtype": vtype, "relatedTo": relatedTo}

    def addConstraint(self, name, variables, coefficients, sense, rhs):
        self.constraints[name] = {
            "variables": variables,
            "coefficients": coefficients,
            "sense": sense,
            "rhs": rhs,
        }

    def addIndicator(
        self, name, indicatorVariable, activeWhen, variable, coefficients, sense, b
    ):
        self.indicators[name] = {
            "indicatorVariable": indicatorVariable,
            "activeWhen": activeWhen,
            "variable": variable,
            "coefficients": coefficients,
            "sense": sense,
            "b": b,
        }

    def addObjective(self, variables, coefficients, sense):
        self.objective = {"variables": variables, "coefficients": coefficients}
        self.sense = sense

    def setSolver(self, solver):
        self.solver = solver

    def addConstraintsFromSparseMatrix(self, matrix, rhs, sense):
        if matrix.shape[0] == len(rhs) == len(sense):
            for i in range(matrix.shape[0]):
                self.addConstraint(
                    "c" + str(i), matrix[i].indices, matrix[i].data, sense[i], rhs[i]
                )
        else:
            raise Exception(
                "Matrix, rhs and sense must have the same length, rhs and sense must be 1D arrays and matrix must be a sparse matrix"
            )

    def sliceVariables(self, variables):
        return {k: self.variables[k] for k in variables}
    
    def modifyConstraint(self, constraintName, variables, coefficients, sense, rhs):
        self.constraints[constraintName]["variables"] = variables
        self.constraints[constraintName]["coefficients"] = coefficients
        self.constraints[constraintName]["sense"] = sense
        self.constraints[constraintName]["rhs"] = rhs
        
    def modifyVariable(self, index, lb=None, ub=None, vtype=None, relatedTo=None):
        if lb is not None:
            self.variables[index]["lb"] = lb
        if ub is not None:
            self.variables[index]["ub"] = ub
        if vtype is not None:
            self.variables[index]["vtype"] = vtype
        if relatedTo is not None:
            self.variables[index]["relatedTo"] = relatedTo
            
    def modifyObjective(self, variables=None, coefficients=None, sense=None):
        if variables is not None:
            self.objective["variables"] = variables
        if coefficients is not None:
            self.objective["coefficients"] = coefficients
        if sense is not None:
            self.sense = sense

    def getProblem(self):
        return {
            "variables": self.variables,
            "constraints": self.constraints,
            "indicators": self.indicators,
            "sense": self.sense,
            "objective": self.objective,
            "solver": self.solver,
        }
        
    def getVariables(self):
        return self.variables

    def getConstraints(self):
        return self.constraints
    
    def getIndicators(self):
        return self.indicators
    
    def getObjective(self):
        return self.objective
    
    def getSense(self):
        return self.sense
    
    def getSolver(self):
        return self.solver
    
    def getVariable(self, index):
        return self.variables[index]
    
    def getConstraint(self, name):
        return self.constraints[name]
    
    def getIndicator(self, name):
        return self.indicators[name]
    
    def getVariableNames(self, variablesIndex=None):
        variableNames = {}
        if variablesIndex is None:
            for index, variable in self.variables.items():
                variableNames[variable["name"]] = index
            return variableNames
        else:
            for index in variablesIndex:
                variableNames[self.variables[index]["name"]] = index
            return variableNames
    
    def interpretProblem(self, parameters=None, verbose=0):  
        if self.solver is None:
            raise Exception("No solver provided, assign a solver to the problem")
        else:
            [interface, _] = setSolver(self.solver)
            self.interface = interface 
        
        self.parameters = parameters
                 
        if self.solver == "gurobi":
            env = self.interface.Env(empty=True)
            if parameters:
                for name, value in parameters.items():
                    try:
                        env.setParam(name, value)
                    except:
                        Warning(f"Parameter {name} could not be set")

            return [gurobiProblemInterpreter(interface=self.interface, env=env, problemDict=self, verbose=verbose), interface]
        elif self.solver == "cplex":
            optModel = interface.Cplex()
            if parameters:
                for name, value in parameters.items():
                    try:
                        code_string = f"optModel.parameters.{name}.set({value})"
                        exec(code_string)
                    except:
                        Warning(f"Parameter {name} could not be set")
            return [cplexProblemInterpreter(interface=interface, optModel=optModel, problemDict=self, verbose=verbose), interface]
        elif self.solver == "scip":
            return [scipProblemInterpreter(interface=interface, problemDict=self, verbose=verbose), interface]
    
    def setModelName(self, name):
        self.modelName = name
        
    def getModelName(self):
        return self.modelName
    
    def setCharacteristicParameter(self, parameter):
        self.characteristicParameter = parameter
    
    def getCharacteristicParameter(self):
        return self.characteristicParameter
    
    def getRelatedVariable(self, index):
        return self.variables[index]["relatedTo"]
    
    def getRelatedVariables(self):
        relatedVariables = {}
        for index, variable in self.variables.items():
            if variable["relatedTo"] is not None:
                variableInfo = {}
                variableInfo['index'] = index
                variableInfo['relatedTo'] = variable['relatedTo']
                relatedVariables[variable['name']] = variableInfo
        return relatedVariables
    
    def getContainedVariables(self, key):
        variablesOfInterest = []
        for index, variable in self.variables.items():
            if variable["relatedTo"] is not None:
                if key.issubset(variable["relatedTo"]): 
                    variablesOfInterest.append(variable["name"])
        return variablesOfInterest
        
    def getConstraintByName(self, name):
        return self.constraints[name]
            
    def setLogPath(self, logPath=str):
        self.logPath = logPath
        
    def getLogPath(self):
        return self.logPath

    def getTimeLimitParameter(self):
        if self.parameters is not None:
            if self.solver == "gurobi":
                return self.parameters["TimeLimit"]
            elif self.solver == "cplex":
                return self.parameters["timelimit"]
        return None
            
    

        
    