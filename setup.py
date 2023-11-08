from distutils.core import setup
setup(
  name = 'gMCSpy',         
  packages = ['gMCSpy'],   
  version = '0.21',      
  license='MIT',        
  description = 'gMCSpy is a python package for the calculation of Genetic Minimal Cut sets (GMCS). In simple terms the idea is to take a metabolic model and calculate the genetic vulnerabilities that will render the biomass production impossible. This is done through a Mixed-Integer Linear problem (MILP) formultion and the use of a linear solver.',   
  author = 'Carlos Javier Rodriguez',                   
  author_email = 'cjrodriguezf@unav.es',      
  url = 'https://github.com/PlanesLab/gMCSpy',   
  download_url = 'https://github.com/PlanesLab/gMCSpy/archive/refs/tags/v.0.2.tar.gz',    
  keywords = ['gMCS', 'Genetic Minimal Cut Sets', 'Gurobi', 'CPLEX'],   
  install_requires=[            
          'cobra',
          'beartype',
          'scipy',
          'tqdm',
          'joblib',
          'gurobipy',
          'cplex'
          
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',      # Define that your audience are developers
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3.8',      
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',      
    'Programming Language :: Python :: 3.11',
  ],
)