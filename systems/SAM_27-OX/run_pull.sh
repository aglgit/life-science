#!/bin/bash
set -e

mkdir -p confs/
cd run/

gmx grompp -c pull_eql.gro -p topol.top -f ../mdp/pull_prd.mdp -o pull_prd.tpr -n index.ndx
gmx mdrun -deffnm pull_prd

echo 0 | gmx trjconv -s pull_prd.tpr -f pull_prd.xtc -o ../confs/conf.gro -sep

exit 1

