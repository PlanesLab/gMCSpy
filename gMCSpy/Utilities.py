from cobra import Model, Reaction, Metabolite
from beartype import beartype
import os
from datetime import datetime
import pandas as pd
import re

def createToyModel():
    """Create a toy model for testing purposes.
    Creates a model with 9 reactions and 5 metabolites, the biomass reaction is **r8** based on the model in
    this `paper <https://academic.oup.com/bioinformatics/article/32/13/2001/1742871?login=true#394386436>`_.

    .. image:: ./figs/toymodel.jpeg
        :height: 350px
        :width: 400 px
        :align: center




    """
    model = Model("ToyModel")
    metaboliteA = Metabolite("m1", name="metaboliteA", compartment="c")
    metaboliteB = Metabolite("m2", name="metaboliteB", compartment="c")
    metaboliteC = Metabolite("m3", name="metaboliteC", compartment="c")
    metaboliteD = Metabolite("m4", name="metaboliteD", compartment="c")
    metaboliteE = Metabolite("m5", name="metaboliteE", compartment="c")

    reaction1 = Reaction("r1")
    reaction1.name = "reaction1"
    reaction1.subsystem = "subsystem1"
    reaction1.lower_bound = 0
    reaction1.upper_bound = 1000.0
    reaction1.add_metabolites({metaboliteA: 1.0})
    reaction1.gene_reaction_rule = "( gene1 and gene2 and gene3)"

    reaction2 = Reaction("r2")
    reaction2.name = "reaction2"
    reaction2.subsystem = "subsystem1"
    reaction2.lower_bound = -1000.0
    reaction2.upper_bound = 1000.0
    reaction2.add_metabolites({metaboliteA: -1.0, metaboliteD: 1.0})
    reaction2.gene_reaction_rule = "( gene1 or gene2 or gene4 or gene5 or gene6)"

    reaction3 = Reaction("r3")
    reaction3.name = "reaction3"
    reaction3.subsystem = "subsystem1"
    reaction3.lower_bound = 0.0
    reaction3.upper_bound = 1000.0
    reaction3.add_metabolites({metaboliteC: -1.0, metaboliteD: 1.0})
    reaction3.gene_reaction_rule = "( gene3 or (gene5 and gene6))"

    reaction4 = Reaction("r4")
    reaction4.name = "reaction4"
    reaction4.subsystem = "subsystem1"
    reaction4.lower_bound = 0.0
    reaction4.upper_bound = 1000.0
    reaction4.add_metabolites({metaboliteA: -1.0, metaboliteC: 1.0})
    reaction4.gene_reaction_rule = '((gene1 and gene2 and gene3 and gene4 and gene5) or (gene2 and gene3 and gene4 and gene6 and gene5) or (gene1 and gene2 and gene3 and gene7 and gene5) or (gene8 and gene2 and gene3 and gene4 and gene5) or (gene8 and gene2 and gene3 and gene7 and gene5) or (gene9 and gene2 and gene3 and gene7 and gene5) or (gene2 and gene3 and gene7 and gene6 and gene5) or (gene9 and gene2 and gene3 and gene4 and gene5)) and ((gene10 and gene11 and gene12 and gene13 and gene14 and gene15 and gene16 and gene17) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene14 and gene16 and gene17) or (gene10 and gene11 and gene12 and gene13 and gene14 and gene19 and gene15 and gene17) or (gene10 and gene11 and gene12 and gene13 and gene20 and gene14 and gene15 and gene16) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene20 and gene14 and gene16) or (gene10 and gene11 and gene12 and gene13 and gene15 and gene16 and gene21 and gene17) or (gene10 and gene11 and gene12 and gene13 and gene19 and gene15 and gene21 and gene17) or (gene10 and gene11 and gene12 and gene13 and gene20 and gene15 and gene16 and gene21) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene20 and gene22 and gene16) or (gene10 and gene11 and gene12 and gene13 and gene19 and gene15 and gene22 and gene17) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene20 and gene19 and gene21) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene19 and gene21 and gene17) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene20 and gene16 and gene21) or (gene10 and gene11 and gene12 and gene13 and gene15 and gene22 and gene16 and gene17) or (gene10 and gene11 and gene12 and gene13 and gene20 and gene15 and gene22 and gene16) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene16 and gene21 and gene17) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene22 and gene16 and gene17) or (gene10 and gene11 and gene12 and gene13 and gene20 and gene19 and gene15 and gene22) or (gene10 and gene11 and gene12 and gene13 and gene20 and gene19 and gene15 and gene21) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene20 and gene19 and gene22) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene19 and gene22 and gene17) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene14 and gene19 and gene17) or (gene10 and gene11 and gene12 and gene18 and gene13 and gene20 and gene14 and gene19) or (gene10 and gene11 and gene12 and gene13 and gene20 and gene14 and gene19 and gene15))'

    reaction5 = Reaction("r5")
    reaction5.name = "reaction5"
    reaction5.subsystem = "subsystem1"
    reaction5.lower_bound = 0.0
    reaction5.upper_bound = 1000.0
    reaction5.add_metabolites({metaboliteA: -1.0, metaboliteB: 1.0})
    reaction5.gene_reaction_rule = "( gene5 and gene2)"

    reaction6 = Reaction("r6")
    reaction6.name = "reaction6"
    reaction6.subsystem = "subsystem1"
    reaction6.lower_bound = 0.0
    reaction6.upper_bound = 1000.0
    reaction6.add_metabolites({metaboliteB: -1.0, metaboliteC: 1.0})
    reaction6.gene_reaction_rule = "( gene3 )"

    reaction7 = Reaction("r7")
    reaction7.name = "reaction7"
    reaction7.subsystem = "subsystem1"
    reaction7.lower_bound = 0.0
    reaction7.upper_bound = 1000.0
    reaction7.add_metabolites({metaboliteD: -1.0, metaboliteE: 1.0})
    reaction7.gene_reaction_rule = "( gene6 )"

    reaction8 = Reaction("r8")
    reaction8.name = "reaction8"
    reaction8.subsystem = "subsystem1"
    reaction8.lower_bound = 0.0
    reaction8.upper_bound = 1000.0
    reaction8.add_metabolites({metaboliteE: -1.0})
    reaction8.gene_reaction_rule = "( gene8 and gene9 and gene10)"

    reaction9 = Reaction("r9")
    reaction9.name = "reaction9"
    reaction9.subsystem = "subsystem1"
    reaction9.lower_bound = 0.0
    reaction9.upper_bound = 1000.0
    reaction9.add_metabolites({metaboliteB: -1.0})
    reaction9.gene_reaction_rule = "((gene2 and (gene5 or gene6)) or gene7)"

    model.add_reactions(
        [
            reaction1,
            reaction2,
            reaction3,
            reaction4,
            reaction5,
            reaction6,
            reaction7,
            reaction8,
            reaction9,
        ]
    )

    model.objective = "r8"

    return model

@beartype
def saveSolutions(
    problemDict: dict, solutionDict: dict, modelStruct: Model, filename: str, path: str
):
    """
    Saves the solutions to a csv file.

    .. warning::
        This function has been design to save the solution of our optimization problems, it may not work for other problems.

    :param problemDict: Dictionary containing the problem information.
    :param solutionDict: Dictionary containing the solution information.
    :param modelStruct: The cobra Model used to build the problem dictionary.
    :param filename: Name of the file to save the solutions.
    :param path: Path to save the file.


    """
    names = []
    orders = []
    reactionSets = []
    for key, data in solutionDict.items():
        vars = []
        for i in data:
            solPos = problemDict["variables"]["variableNames"].index(str(i))
            corresPondingVar = modelStruct.reactions[
                problemDict["variables"]["variableGroups"]["z"].index(solPos)
            ].id
            vars.append(corresPondingVar)
        orderNumber = key.split("order")[1].split("_")[0]
        names.append(key)
        orders.append(orderNumber)
        reactionSets.append(vars)

    df = pd.DataFrame({"solName": names, "order": orders, "Reactions": reactionSets})
    df.to_csv(path + "/" + filename + ".csv", index=False)

@beartype
def saveSolutionDict(dict: dict, filename: str, path: str):
    dataframe = pd.DataFrame.from_dict(dict, orient="columns")
    dataframe.to_csv(path + "/" + filename + ".csv", index=False)

def setSolver(defaultSolver, useIdicators: bool = True):
    """Set the solver to be used for the optimization problem.
    This function is used automatically set the correct interface based on the solver chosen.
    Current options are 'cplex', 'gurobi' and 'scip'.

    :param defaultSolver: The solver to be used for the optimization problem. Must be one of the options above.
    :param indicator: If True, the solver will be set to use indicators. If False, the solver will be set to not use indicators.
        The function also checks if the solver has the capability to use indicators. If the solver does not have the capability
        it returns turns the indicators off and the Big M method is implemented.

    :return:
        - The interface to be used for the optimization problem.
        - A boolean indicating if the solver will use indicators or not.
    """
    if defaultSolver == "cplex":
        try:
            import cplex

            return [cplex, True]
        except ImportError :
            print(
                "cplex could not be loaded, check if it is installed and the path is correct"
            )
    elif defaultSolver == "gurobi":
        try:
            import gurobipy

            return [gurobipy, True]
        except ImportError :
            print(
                "gurobi could not be loaded, check if it is installed and the path is correct"
            )

    elif defaultSolver == "scip":
        try: 
            import pyscipopt

            return [pyscipopt, True]
        except ImportError :
            print(
                "pyscipopt could not be loaded, check if it is installed and the path is correct"
            )
    else:
        print(
            "The default solver is not supported, please select (cplex, gurobi, scip)"
        )

def createLog(folderLogPath: str = "logs/log"):
    """Creates a log file with the current date and time, and returns the filename.
    This functions also creates the folder to save the log file if it does not exist.
    It is used to save the log of the optimization process depending on the verbosity level."""
    filename = folderLogPath
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    print(filename)
    return filename

def loadSolutions(path: str = None):
    """Load the solutions generated from the calculation of gMCSs.

    :param path: Path to the folder with the solutions, if no path is given the most recent calculation in the logs folder is used.

    :return: A dictionary with the solutions.
    """
    if path is None:
        path = 'logs'
        # Get a list of all directories in the given path
        directories = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
        
        # Sort the directories by their modification time in descending order
        directories.sort(key=lambda x: os.path.getmtime(os.path.join(path, x)), reverse=True)
        
        if directories:
            # The most recent folder is the first one in the sorted list
            most_recent = directories[0]
            complePath = path + "/" + most_recent + "/solutions/solutions.log"
            return read_sol_log(complePath)
            
        else:
            raise Warning("No solutions found in the given path.")
    else:
        return read_sol_log(path)  

def read_sol_log(logFile):
    # Read the log file
    with open(logFile, "r") as file:
        # First line are headers
        lines = file.readlines()[1:]
        solDict = {}
        for line in lines:
            line = line.strip()
            # Split the line by the separator
            lineList = line.split(":")
            # Get the genes
            geneList = add_quotes_inside_brackets(lineList[0])
            geneList = eval(geneList)
            geneList = frozenset(geneList)
            # Get the reactions in the solution
            reactionList = eval(add_quotes_inside_brackets(lineList[1]))
            solDict[geneList] = reactionList
        
        return solDict
        



def add_quotes_inside_brackets(input_string):
    # Use regular expressions to find text inside square brackets
    pattern = r'\[([^]]*)\]'
    
    # Define a function to add quotes around the found text
    def add_quotes(match):
        # Split the text inside brackets by ',' and add quotes to each element
        elements = match.group(1).split(',')
        quoted_elements = [f'"{element.strip()}"' for element in elements]
        return '[' + ','.join(quoted_elements) + ']'
    
    # Use re.sub() to apply the function to the input string
    result = re.sub(pattern, add_quotes, input_string)
    
    return result