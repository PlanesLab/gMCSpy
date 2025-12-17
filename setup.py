from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="gMCSpy",
    version="1.1.5",
    packages=find_packages(),  # include your package
    license="MIT",
    description=(
        "gMCSpy is a python package for the calculation of Genetic Minimal Cut sets (GMCS). "
        "In simple terms, the idea is to take a metabolic model and calculate the genetic "
        "vulnerabilities that will render the biomass production impossible. "
        "This is done through a Mixed-Integer Linear problem (MILP) formulation and the use of a linear solver."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Carlos Javier Rodriguez",
    author_email="cjrodriguezf@unav.es",
    url="https://github.com/PlanesLab/gMCSpy",
    download_url="https://github.com/PlanesLab/gMCSpy/archive/refs/tags/v1.1.5.tar.gz",
    keywords=["gMCS", "Genetic Minimal Cut Sets", "Gurobi", "CPLEX"],
    install_requires=[
        "cobra",
        "beartype",
        "scipy",
        "tqdm",
        "joblib",
        "gurobipy",
        "cplex",
        "bidict"
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
)
