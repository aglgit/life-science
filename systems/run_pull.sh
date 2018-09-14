set -e

mkdir -p run/
cd run/
cp ../topol.top topol.top
cp ../index.ndx index.ndx

gmx editconf -f ../../pdb/SAM_24SOX.pdb -o box.gro -box 4.33 4.50 6.5

gmx grompp -c box.gro -p topol.top -f ../mdp/pull_min.mdp -o pull_min.tpr -n index.ndx
gmx mdrun -deffnm pull_min

gmx grompp -c pull_min.gro -p topol.top -f ../mdp/pull_eql.mdp -o pull_eql.tpr -n index.ndx
gmx mdrun -deffnm pull_eql

gmx grompp -c pull_eql.gro -p topol.top -f ../mdp/pull_prd.mdp -o pull_prd.tpr -n index.ndx
gmx mdrun -deffnm pull_prd

exit 1

