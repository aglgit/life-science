#!/bin/bash
set -e

CONF=$1

mkdir -p umbrella/
mkdir -p run/
cd run/
cp ../topol.top topol.top
cp ../index.ndx index.ndx

gmx grompp -c ../confs/${CONF}.gro -p topol.top -f ../mdp/energy_eql.mdp -o ${CONF}_eql.tpr -n index.ndx
gmx mdrun -deffnm ${CONF}_eql

gmx grompp -c ${CONF}_eql.gro -p topol.top -f ../mdp/energy_prd.mdp -o ${CONF}_prd.tpr -n index.ndx
gmx mdrun -deffnm ${CONF}_prd -pf ../umbrella/${CONF}f.xvg -px ../umbrella/${CONF}x.xvg

cp ${CONF}_prd.tpr ../umbrella/${CONF}_prd.tpr

exit 1

