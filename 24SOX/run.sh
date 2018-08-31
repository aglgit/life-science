set -e

L=2.3
NMOL=10

cd run/
gmx insert-molecules -ci ../24SOX.pdb -o box.gro -box $L $L $L -nmol $NMOL
gmx grompp -f ../mdp/min.mdp -c box.gro -p ../topol.top -o min.tpr
# gmx genion -s min.tpr -o conf.gro -p min.top -neutral
gmx mdrun -deffnm min
cd ../
