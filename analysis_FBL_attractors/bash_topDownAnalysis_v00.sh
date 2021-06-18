#!/bin/bash

#SBATCH --nodes=1 # request one node
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=1G # Maximum amount of memory this job will be given, try to estimate this to the best of your ability. This asks for 1 GB of ram.
#SBATCH --time=0-00:59:00 # ask that the job be allowed to run for 30 minutes.
#SBATCH --array=0-99%50 #specify how many times you want a job to run, we have a total of 7 array spaces

# everything below this line is optional, but are nice to have quality of life things
cd /work/LAS/ckadelka-lab/BooleanGRNs

# under this line, we can load any modules if necessary
module load gcc/7.3.0-xegsmw4
module load python/3.7.8-vzxxni7

# activate python virtual environment
. env3/bin/activate

#below this line is where we can place our commands, in this case it will just simply output the task ID of the array
python FBL_to_attractors.py $SLURM_ARRAY_TASK_ID 3

