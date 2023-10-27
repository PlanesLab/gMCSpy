Examples
===========================

Here we present a set of example of how to use gMCSpy with genome-scale metabolic models (GEMS). 

.. toctree::
   :maxdepth: 2

|

Calculating the Genetic Minimal Cut Sets (gMCS) of E. coli core GEM
---------------------------------------------------------------------------

To replicate this example, make sure you have completed the `installation guide <installation.rst>`_, the E. coli core model can be downloaded from the `E. coli core Model <http://bigg.ucsd.edu/models/e_coli_core>`_. The following code-block shows how to calculate the gMCS of the E. coli core model.

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateGeneMCS

    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load the E. coli core model
    mini_mat_path = (Path(".") / "./data/e_coli_core.mat")
    model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
    # Setting the id of the model to "e_coli_core" is not necessary, but makes it easier to identify the results
    model_ecoli.id = "e_coli_core"

    # Calculate the gMCS of the E. coli core model
    # By default, the gMCS are calculated using the gurobi solver, 
    # with maxKOLength=3 (i.e. max number of genes in a gMCS is 3)) and with a maximum of 1e6 solutions
    calculateGeneMCS(model_ecoli)

    # The results are saved in the working directory in a folder named "logs"

    # The results can be retrieved by loading the results file
    from gMCSpy.Utilities import loadSolutions

    # If no path is specified, the most recent calculation in the logs folder is used.
    solutions = loadSolutions()
    # Or you could specify the path to the results file you desire
    solutions = loadSolutions("logs/gMCSpy_ToyModel_gurobi_20231024-130644/solutions/solutions.log")


|
| 

Calculating the Minimal Cut Sets (MCS) of E. coli core GEM
---------------------------------------------------------------------------

To replicate this example, make sure you have completed the `installation guide <installation.rst>`_, the E. coli core model can be downloaded from the `E. coli core Model <http://bigg.ucsd.edu/models/e_coli_core>`_. The following code-block shows how to calculate the MCS of the E. coli core model.

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateMCS

    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load the E. coli core model
    mini_mat_path = (Path(".") / "./data/e_coli_core.mat")
    model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
    model_ecoli.id = "e_coli_core"

    # Calculate the MCS of the E. coli core model
    # Gurobi is the default solver, but you can also use CPLEX or GLPK
    mcs =  calculateMCS(model_ecoli, MAX_LENGTH_MCS=3, MAX_MCS=1000)
    print(mcs)

|
| 

Calculating the Nutrient Genetic Minimal Cut Sets (ngMCS) of E. coli core GEM
------------------------------------------------------------------------------------

The concept of nutrient-genetic Minimal Cut Sets is a novel analysis tool in the study of GEMs. These sets are crucial for determining the nutrients and genes necessary to sustain the production of biomass within a system. When the parameter "isNutrient" is set to true, it includes all exchange reactions as genes. Exchange reactions represent the uptake or secretion of metabolites, we only include uptakes in this analysis, mirroring the intake of nutrients and the release of products by an organism. By identifying the minimal cut sets, researchers can pinpoint the combinations of nutrients and genetic elements that are necessary for the synthesis of biomass, shedding light on the interplay between genetics and metabolism. This insight is valuable for understanding the requirements of an organism and has applications in biotechnology and medicine.
The E. coli core model can be downloaded from the `E. coli core Model <http://bigg.ucsd.edu/models/e_coli_core>`_.

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateGeneMCS
    from gMCSpy.Utilities import loadSolutions
    
    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load the E. coli core model
    mini_mat_path = (Path(".") / "./data/e_coli_core.mat")
    model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
    # Setting the id of the model to "e_coli_core" is not necessary, but makes it easier to identify the results
    model_ecoli.id = "e_coli_core"

    # When using isNutrient=True, all the exchange reactions are included in the model as artificial genes. 
    # Therefore, we can calculate which nutrients are required in conjunction with which genes to maintain
    # the production of biomass. 
    calculateGeneMCS(model_ecoli, isNutrient=True)


    # The results are saved in the working directory in a folder named "logs"
    # The results can be retrieved by loading the results file
    solutions = loadSolutions()


|
| 

Calculating the Genetic Minimal Cut Sets (gMCS) of E. coli core GEM targeting a specific set of genes
-------------------------------------------------------------------------------------------------------

Targeting a specific set of genes is a common use case for gMCSpy. This can be done by specifying the target genes in the parameter "targetKOs". The E. coli core model can be downloaded from the `E. coli core Model <http://bigg.ucsd.edu/models/e_coli_core>`_.

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateGeneMCS
    from gMCSpy.Utilities import loadSolutions
    
    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load the E. coli core model
    mini_mat_path = (Path(".") / "./data/e_coli_core.mat")
    model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
    # Setting the id of the model to "e_coli_core" is not necessary, but makes it easier to identify the results
    model_ecoli.id = "e_coli_core"

    # We can target a specific set of genes by specifying the target genes in the parameter "targetKOs"
    #gene: b2284, name: nuoF, function: NADH16/electron transport chain in bacteria
    targets = ['b2284']
    calculateGeneMCS(model_ecoli, targetKOs=targets)


    # The results are saved in the working directory in a folder named "logs"
    # The results can be retrieved by loading the results file
    solutions = loadSolutions()

The resulting gene KOs are presented in the following table, we can see that the gMCS are including targeted gene and filling the set with genes that in conjunction are required for proliferation.

.. list-table:: 
   :widths: 50 50
   :header-rows: 1
   :align: center
   
   * - Gene
     - Reactions
   * - b2284, b3919
     - NADH16, TPI
   * - b2284, b4025
     - PGI, NADH16
   * - b0902, b2284, b3952
     - PFL, NADH16
   * - b0902, b2284, b3951
     - PFL, NADH16
   * - b1723, b3916, b2284
     - NADH16, PFK
   * - b2492, b2284, b0904
     - FORt2, NADH16, FORt
   * - b2284, b3734, b0116
     - AKGDH, NADH16, ATPS4r, PDH
   * - b2284, s0001, b3734
     - H2Ot, ATPS4r, NH4t, NADH16, O2t, CO2t, ACALDt
   * - b3733, b2284, s0001
     - H2Ot, ATPS4r, NH4t, NADH16, O2t, CO2t, ACALDt
   * - b2284, s0001, b3732
     - H2Ot, ATPS4r, NH4t, NADH16, O2t, CO2t, ACALDt
   * - b2284, s0001, b3736
     - H2Ot, ATPS4r, NH4t, NADH16, O2t, CO2t, ACALDt
   * - b2284, s0001, b3737
     - H2Ot, ATPS4r, NH4t, NADH16, O2t, CO2t, ACALDt


|
| 

Validating the Genetic Minimal Cut Sets (gMCS) of E. coli core GEM targeting a specific set of genes
-------------------------------------------------------------------------------------------------------

To validate that the actual solutions are minimal cut sets we have include a function that blocks the reactions associated with those genes and runs an FBA. If the growth rate is zero, the solution is a valid gMCS. 

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateGeneMCS
    from gMCSpy import checkGMCSParallel
    from gMCSpy.Utilities import loadSolutions

    
    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load the E. coli core model
    mini_mat_path = (Path(".") / "./data/e_coli_core.mat")
    model_ecoli = load_matlab_model(str(mini_mat_path.resolve()))
    # Setting the id of the model to "e_coli_core" is not necessary, but makes it easier to identify the results
    model_ecoli.id = "e_coli_core"

    calculateGeneMCS(model_ecoli, solver="gurobi", maxKOLength=2)

    # The results are saved in the working directory in a folder named "logs"
    # The results can be retrieved by loading the results file
    solutions = loadSolutions()

    # Check the solutions generates its own and saves a log file each time is called
    checkGMCSParallel(solutions, model_ecoli, workers=16)

    # Another option is setting the parameter "checkSolutions" to True in the calculateGeneMCS function
    # The checks would be included in the log file of the calculation
    calculateGeneMCS(model_ecoli, solver="gurobi", maxKOLength=2, checkSolutions=True)







|
| 

Advance Example, Calculating gMCS of Human-GEM with custom parameters
----------------------------------------------------------------------

In this example we will calculate all gMCSs up to length 4 of the `Human-GEM <https://github.com/SysBioChalmers/Human-GEM/releases>`_, using cplex and increasing the time limit to 1e5.
All parameters can be found in the `calculateGeneMCS <gMCSpy.calculateGMCS.html#gMCSpy.calculateGMCS.calculateGeneMCS>`_ function.

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateGeneMCS
    from gMCSpy.Utilities import loadSolutions
    
    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load the E. coli core model
    mini_mat_path = (Path(".") / "./data/humanGEM.mat")
    model_human = load_matlab_model(str(mini_mat_path.resolve()))
    # Setting the id of the model to "e_coli_core" is not necessary, but makes it easier to identify the results
    model_human.id = "humanGEM"

    # Additional parameters can be specified in the calculateGeneMCS function
    calculateGeneMCS(model_human, 
                     solver="cplex", 
                     maxKOLength=4, 
                     timeLimit=1e5)


    # The results are saved in the working directory in a folder named "logs"
    # The results can be retrieved by loading the results file
    solutions = loadSolutions()