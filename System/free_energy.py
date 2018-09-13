import pandas as pd
import numpy as np
import os

names = ['conf', 'dist']
dtype={'conf': np.int32, 'dist': np.float64}
df = pd.read_csv('run/summary_distances.dat', names=names, delim_whitespace=True,
                 dtype=dtype)
print(df)

dist = df['dist'].values
distances = np.abs(dist[1:] - dist[:-1])
print(distances)
cumsum_distances = np.cumsum(distances)
print(np.insert(cumsum_distances, 0, 0))

confs = [0, 4, 8, 11, 14, 16, 18, 19]
conf_names = ['conf%d' % conf for conf in confs]

os.chdir('run/')
for name in conf_names:
    os.system('gmx grompp -c %s.gro -p topol.top -f ../mdp/energy_eql.mdp -o %s.tpr -n index.ndx' % (name, name))
    os.system('gmx mdrun -deffnm %s -pf %sf.xvg -px %sx.xvg' % (name, name, name))

    os.system('gmx grompp -c %s.gro -p topol.top -f ../mdp/energy_prd.mdp -o %s.tpr -n index.ndx' % (name, name))
    os.system('gmx mdrun -deffnm %s -pf %sf.xvg -px %sx.xvg' % (name, name, name))

