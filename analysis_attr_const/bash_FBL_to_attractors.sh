#!/bin/bash

#SBATCH --nodes=1 # request node
#SBATCH --cpus-per-task=50  # ask for cpu
#SBATCH --mem=31G # Maximum amount of memory this job will be given, try to estimate this to the best of your ability.
#SBATCH --time=03:00:00 # ask that the job be allowed to run for time.
#SBATCH --array=0-99%50 #specify how many times you want a job to run
#--partition=biocrunch

# everything below this line is optional, but are nice to have quality of life things
cd /work/LAS/ckadelka-lab/BooleanGRNs/analysis_FBL_attractors

# under this line, we can load any modules if necessary
module purge
module load gcc/7.3.0-xegsmw4
module load python/3.7.8-vzxxni7

# set virtual environment
python -mvenv env3

# activate python virtual environment
. env3/bin/activate

#install dependencies
pip install numpy
pip install scipy
pip install matplotlib
pip install networkx
pip install sympy
pip install pandas

#below this line is where we can place our commands
python Black_Computational_Analysis_v01.py $SLURM_ARRAY_TASK_ID

