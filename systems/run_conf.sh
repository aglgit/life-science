set -e

CONF=$1

mkdir -p force/
cd run/
cp ../topol.top topol.top
cp ../index.ndx index.ndx

gmx grompp -c ../confs/${CONF}.gro -p topol.top -f ../mdp/energy_eql.mdp -o ${CONF}_eql.tpr -n index.ndx
gmx mdrun -deffnm ${CONF}_eql

gmx grompp -c ${CONF}_eql.gro -p topol.top -f ../mdp/energy_prd.mdp -o ${CONF}_prd.tpr -n index.ndx
gmx mdrun -deffnm ${CONF}_prd -pf ../force/${CONF}f.xvg -px ${CONF}x.xvg

exit 1

