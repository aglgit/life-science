#!/bin/bash

#SBATCH --job-name=CONF
#SBATCH --account=nn4654k
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1 --ntasks-per-node=4

## Set up job environment
source /cluster/bin/jobsetup

# Run job
module load gromacs
./run_conf.sh CONF

exit 1

