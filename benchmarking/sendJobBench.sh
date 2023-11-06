#!/bin/bash
###################################################################################
### sendgmcspy
###################################################################################
#
#

#SBATCH --job-name=GMCSpy
#SBATCH --partition=xlong
#SBATCH --mail-type=END
#SBATCH --mail-user=cjrodriguezf@unav.es
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=atlas-316
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=144:00:00
#SBATCH -o /scratch/a905383/logs/gmcspy/bench_%j.out

echo ===================================
echo ===     Load the Packages       ===
echo ===================================
echo `date`

module load Miniconda3/py37_4.8.2

source activate /scratch/a905383/gmcspy/conda

python -c 'import sys; print(sys.version_info[:])'

echo $PATH

export PATH="/scratch/a905383/gmcspy/conda/bin:$PATH"

python -c 'import sys; print(sys.version_info[:])'

export LD_LIBRARY_PATH="/scratch/a905383/gmcspy/conda/lib:$LD_LIBRARY_PATH"

echo $LD_LIBRARY_PATH

export GRB_LICENSE_FILE=/scratch/a905383/gmcspy/gurobi.lic

export JAVA_HOME=/scratch/a905383/Software/java/jdk1.8.0_371

echo ===================================
echo ===   Run Python script   ===
echo ===================================

python /scratch/a905383/gmcspy/Gemcuts/benchGMCSpy.py


echo ===================================
echo ===   Run Matlab   ===
echo ===================================

module load MATLAB/R2020b
export ILOG_CPLEX_PATH="/scratch/a905383/Software/cplex1210/cplex/matlab/x86-64_linux"
cd /scratch/a905383/gmcspy/matlab_codes/
matlab -nojvm -nodisplay -nosplash -r "run('benchGMCS.m'); exit"

echo ===================================
echo ===   Run Python script   ===
echo ===================================

python /scratch/a905383/gmcspy/Gemcuts/benchSD.py