Installation Guide 
====================

gMCSpy is a Python package that provides a framework for finding genetic minimal cut sets (gMCS), it requires optimization solvers; this installation guide will walk you through the process of installing CPLEX, Gurobi and SCIP (with PySCIPOpt)

Creating a conda environment:
-----------------------------
.. code-block:: console

    conda create -n gmcspy python=3.9

    pip install gmcspy


Installing CPLEX:
-----------------

1. IBM CPLEX provides a free license for academic users. To obtain a license, visit the IBM Academic Initiative website: `CPLEX <https://community.ibm.com/community/user/ai-datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students>`_ 

2. Download the appropriate version of CPLEX for your operating system.

3. Follow the installation instructions provided for your specific platform.

4. Locate the folder where CPLEX was installed. Example: `C:\Program Files\IBM\ILOG\CPLEX_Studio2211\python`

5. Using the python in your environment run the following command:

.. code-block:: console

    python C:/Program Files/IBM/ILOG/CPLEX_Studio2211/python/setup.py install

Installing Gurobi:
------------------

1. Visit the Gurobi download page: https://www.gurobi.com/downloads/

2. Download the Gurobi Optimizer for your operating system.

3. Follow the installation instructions provided for your specific platform, `Gurobi Installation Guide <https://www.gurobi.com/documentation/quickstart.html>`_

4. Make sure to activate your Gurobi license.

Installing SCIP (with conda for PySCIPOpt):
--------------------------------------------

1. In your environment, run the following command:

.. code-block:: console

    conda install PySCIPOpt=4.3.0
   
