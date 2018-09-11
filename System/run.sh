set -e

mkdir -p run/
cd run/
cp ../topol.top topol.top
cp ../index.ndx index.ndx

gmx editconf -f ../system.pdb -o box.gro -box 4.33 4.50 6.5

gmx grompp -c box.gro -p topol.top -f ../mdp/min.mdp -o min.tpr -n index.ndx
gmx mdrun -deffnm min

gmx grompp -c min.gro -p topol.top -f ../mdp/eql.mdp -o eql.tpr -n index.ndx
gmx mdrun -deffnm eql

gmx grompp -c eql.gro -p topol.top -f ../mdp/prd.mdp -o prd.tpr -n index.ndx
gmx mdrun -deffnm prd

gmx trjconv -f prd.xtc -s prd.tpr -pbc mol -o prd-mol.xtc

