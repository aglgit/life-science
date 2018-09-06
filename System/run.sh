set -e

L=10.0
NMOL=10

mkdir -p run/
cd run/
cp ../topol.top topol.top

gmx editconf -f ../ch3_thr_nopep.pdb -o box.gro -c -d 1.0 -bt cubic
gmx insert-molecules -ci ../../24SOX/24SOX.pdb -o box.gro -f box.gro -nmol $NMOL

gmx grompp -f ../mdp/min.mdp -p topol.top -c box.gro
gmx genion -s topol.tpr -o box.gro -p topol.top -neutral

gmx grompp -c box.gro -p topol.top -f ../mdp/min.mdp -o min.tpr
gmx mdrun -deffnm min 

gmx grompp -c box.gro -p topol.top -f ../mdp/eql.mdp -o eql.tpr
gmx mdrun -deffnm eql 

gmx grompp -c box.gro -p topol.top -f ../mdp/eql2.mdp -o eql.tpr
gmx mdrun -deffnm eql 

gmx grompp -c box.gro -p topol.top -f ../mdp/prd.mdp -o prd.tpr
gmx mdrun -deffnm prd 
