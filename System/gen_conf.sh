set -e

mkdir -p run/
cd run/
cp ../topol.top topol.top
cp ../index.ndx index.ndx

gmx editconf -f ../24SOX_C12.pdb -o box.gro -box 4.33 4.50 7.0

gmx grompp -c box.gro -p topol.top -f ../mdp/conf_min.mdp -o conf_min.tpr -n index.ndx
gmx mdrun -deffnm conf_min

gmx grompp -c conf_min.gro -p topol.top -f ../mdp/conf_eql.mdp -o conf_eql.tpr -n index.ndx
gmx mdrun -deffnm conf_eql

gmx grompp -c conf_eql.gro -p topol.top -f ../mdp/conf_pull.mdp -o conf_pull.tpr -n index.ndx
gmx mdrun -deffnm conf_pull

../get_distances.sh

