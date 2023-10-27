# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'gMCSpy'
copyright = '2023, Carlos Rodriguez'
author = 'Carlos Rodríguez, Naroa Barrena, Danel Olaverri-Mendizabal, Luis Valcárcel, and Francisco J. Planes'

# The package is one up directory
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 
              'sphinx.ext.napoleon',
              'sphinx.ext.coverage',
              'sphinx.ext.autosummary',
               'sphinx.ext.autosectionlabel',
               'sphinx.ext.viewcode']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

add_module_names = False

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['my_theme.css']
