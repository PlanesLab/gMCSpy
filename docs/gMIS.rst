Genetic Minimal Intervention Sets (gMIS)
===========================

Here we present a set of examples of how to use gMCSpy with genome-scale metabolic models (GEMS) and Regulatory networks. 

.. toctree::
   :maxdepth: 2

|

Calculating the gMIS of Human-GEM v.14.0
---------------------------------------------------------------------------

To replicate this example, make sure you have completed the `installation guide <installation.rst>`_, The Human-GEM model can be downloaded from the `Human-GEM <https://github.com/SysBioChalmers/Human-GEM/releases>`_. We also need a regulatory network, this will allow us to integrate activation and blockage intereactions. The following code-block shows how to calculate the gMIS of the Human-GEM model.

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateGMIS

    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load model
    mini_mat_path = (Path(".") / "./data/human-gem-14.mat")
    model_human = load_matlab_model(str(mini_mat_path.resolve()))
    # Setting the id of the model to "humanGEM-14" is not necessary, but makes it easier to identify the results
    model_human.id = "humanGEM-14"

    # Load regulatory network
    regulatory_network = pd.read_csv(Path('signaling_data', 'Omnipath.txt'), sep=' ')

    # Calculate the gMIS of the Human-GEM model
    # By default, the gMIS are calculated using the gurobi solver, 
    # with maxKOLength=3 (i.e. max number of interventions in a set is 3)) and with a maximum of 1e6 solutions
    # Another important parameter is the number of layers, this parameter is used to control how many layers of the regulatory network are used in the calculation
    # 1 layer means that only the direct interactions are considered, 2 layers means that the direct interactions and the interactions of the direct interactors # are considered, and so on.
    calculateGMIS(model_human, regulatory_dataframe=regulatory_network, solver='gurobi', num_layers=2, maxKOLength=3)
    # The results are saved in the working directory in a folder named "logs"

    # The results can be retrieved by loading the results file
    from gMCSpy.Utilities import loadSolutions

    # If no path is specified, the most recent calculation in the logs folder is used.
    solutions = loadSolutions()
    # Or you could specify the path to the results file you desire
    solutions = loadSolutions("logs/gMCSpy_ToyModel_gurobi_20231024-130644/solutions/solutions.log")


|
| 

Calculating the readjustment of Human-GEM v.14.0 interventions
---------------------------------------------------------------------------

When calculating the gMIS, we can also calculate the readjustment of the interventions. Sometimes an intervention of the set cannot be knocked down, or the knocked down of another intervention in the set causes the activation of the intervention. In simple terms, not all interventions are independent and we need to check if the interventions are readaptable.
The following code-block shows how to calculate the readjustment of the interventions of the Human-GEM model.

.. code-block:: python

    from cobra.io import load_matlab_model
    from pathlib import Path
    from gMCSpy import calculateReadjustment

    # We recommend setting a working directory to easily find the results
    import os
    os.chdir("path/to/working/directory")
    
    # Load model
    mini_mat_path = (Path(".") / "./data/human-gem-14.mat")
    model_human = load_matlab_model(str(mini_mat_path.resolve()))
    # Setting the id of the model to "humanGEM-14" is not necessary, but makes it easier to identify the results
    model_human.id = "humanGEM-14"

    # Load regulatory network
    regulatory_network = pd.read_csv(Path('signaling_data', 'Omnipath.txt'), sep=' ')

    # Calculate the readjustment of the interventions of the Human-GEM model

    # We need to load al the gMIS calculated previously, or a custome list of gMIS
    gMIS = [['ENSG00000001630_KO', 'ENSG00000132196_KO'], # not readptation
        ['ENSG00000002330_KI', 'ENSG00000145284_KO'], # not readptation
        ['ENSG00000004487_KO', 'ENSG00000145284_KO'], # not readptation
        ['ENSG00000027075_KI', 'ENSG00000082701_KO'], 
        ['ENSG00000067606_KI', 'ENSG00000082701_KO'],
        ['ENSG00000072062_KI', 'ENSG00000082701_KO'],
        ['ENSG00000072062_KI', 'ENSG00000141736_KO'],
        ['ENSG00000078061_KI', 'ENSG00000082701_KO'],
        ['ENSG00000080824_KO', 'ENSG00000145284_KO'],]

    calculateReadjustment(gMISList=gMIS, model=model_human, regulatory_dataframe=regulatory_network, num_layers=2, solver='gurobi')
|
| 

**Note:** Only the interventions with readjustments are shown in the results. The results are saved in the working directory in a folder named "logs"

E.g. 
{'ENSG00000027075_KI': {'obj': 1.0, 'readjustment': True}, 'ENSG00000082701_KO': {'obj': None}, 'readjustment': True}
{'ENSG00000067606_KI': {'obj': 1.0, 'readjustment': True}, 'ENSG00000082701_KO': {'obj': None}, 'readjustment': True}
{'ENSG00000072062_KI': {'obj': 1.0, 'readjustment': True}, 'ENSG00000082701_KO': {'obj': None}, 'readjustment': True}
{'ENSG00000072062_KI': {'obj': 1.0, 'readjustment': True}, 'ENSG00000141736_KO': {'obj': None}, 'readjustment': True}
{'ENSG00000078061_KI': {'obj': 1.0, 'readjustment': True}, 'ENSG00000082701_KO': {'obj': None}, 'readjustment': True}
