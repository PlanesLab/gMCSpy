.. raw :: html
   <head>
   <meta name="google-site-verification" content="XAVEwjTjhotOHuQDSJKqwWTR1AGqwel2lCMHWnabMKI" />
   </head>

=============================================
Genetic Minimal Cut Sets in python (gMCSpy)
=============================================

gMCSpy is a Python package that provides tools for computing genetic minimal cut sets (gMCSs) in large metabolic networks.
gMCSs are minimal sets of gene knockouts that disable a desired metabolic function, such as biomass production or target product secretion. 
gMCSs can be used to identify metabolic vulnerabilities in cancer cells or to design optimal gene knockout strategies for metabolic engineering.

gMCS is based on the concept of minimal cut sets (MCSs), which are minimal sets of reaction knockouts that disable a desired metabolic function. 
However, gMCS takes into account the gene-protein-reaction (GPR) rules that link reactions to genes in genome-scale metabolic models.
Therefore, gMCS computes MCSs at the gene level, which are more biologically relevant and feasible than reaction-level MCSs.

gMCSpy is compatible with the COBRApy package, which is a popular Python library for constraint-based modeling of metabolic networks.
gMCSpy can read and write metabolic models in all formats supported by COBRApy, such as SBML, JSON, and MAT. 


**Authors: Carlos Rodríguez, Naroa Barrena, Danel Olaverri-Mendizabal, Idoia Ochoa, Luis Valcárcel, and Francisco J. Planes** 



.. important::

    This sample documentation was generated on |today|.


Quickstart
==========

1. Install gMCSpy with pip:

   .. code-block:: console

       $ pip install gMCSpy

2. Check our tutorial to learn how to use gMCSpy:

   **gMCS of E. coli core** :doc:`examples` 

.. toctree::
   :maxdepth: 1
   :caption: Documentation:

   modules
   examples
   installation





   