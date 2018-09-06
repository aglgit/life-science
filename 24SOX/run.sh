set -e

L=5.0
NMOL=1

mkdir -p run/
cd run/
cp ../topol.top topol.top

gmx insert-molecules -ci ../24SOX.pdb -o box.gro -box $L $L $L -nmol $NMOL
gmx solvate -cp box.gro -o box_water.gro -cs tip4p -p topol.top

gmx grompp -f ../mdp/min.mdp -p topol.top -c box_water.gro
gmx genion -s topol.tpr -o box_water.gro -p topol.top -neutral

gmx grompp -c box_water.gro -p topol.top -f ../mdp/min.mdp -o min.tpr
gmx mdrun -deffnm min

gmx grompp -c box_water.gro -p topol.top -f ../mdp/eql.mdp -o eql.tpr
gmx mdrun -deffnm eql

gmx grompp -c box_water.gro -p topol.top -f ../mdp/prd.mdp -o prd.tpr
gmx mdrun -deffnm prd

gmx trjconv -f prd.xtc -s prd.tpr -pbc mol -o prd-mol.xtc
