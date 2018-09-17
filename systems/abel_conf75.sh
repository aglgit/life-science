#!/bin/bash

#SBATCH --job-name=CONF75
#SBATCH --account=nn4654k
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1 --ntasks-per-node=1

## Set up job environment
source /cluster/bin/jobsetup

# Run job
module load gromacs
./run_conf.sh conf75

