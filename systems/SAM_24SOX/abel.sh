#!/bin/bash
set -e

CONF=$1

#SBATCH --job-name=$CONF
#SBATCH --account=nn4654k
#SBATCH --time=05:30:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1 --ntasks-per-node=1

## Set up job environment
source /cluster/bin/jobsetup

# Run job
module load gromacs
./run_conf.sh $CONF

exit 1

