#!/bin/bash

#SBATCH --job-name=CONF
#SBATCH --account=nn4654k
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1 --ntasks-per-node=16
#SBATCH --output=log/slurm-CONF.out

## Set up job environment
source /cluster/bin/jobsetup

# Run job
module load gromacs
./run_conf.sh CONF

exit 1

