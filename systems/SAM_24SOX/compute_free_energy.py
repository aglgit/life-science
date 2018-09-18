import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from scipy.integrate import simps

files = ['force/%s' % f for f in os.listdir('force') if not f.startswith('.')]
ordered_files = sorted(files, key=lambda x: (int(re.sub('\D','',x)),x))

mean_force = np.zeros(len(ordered_files))
std_force = np.zeros(len(ordered_files))

for i, file in enumerate(ordered_files):
    with open(file, 'r') as infile:
        for line in infile:
            if '@TYPE xy' in line:
                break
        names = ['time', 'force']
        df = pd.read_csv(infile, delim_whitespace=True, names=names)
        mean_force[i] = df['force'].mean()
        std_force[i] = df['force'].std()

files = ['pos/%s' % f for f in os.listdir('pos') if not f.startswith('.')]
ordered_files = sorted(files, key=lambda x: (int(re.sub('\D','',x)),x))

mean_pos = np.zeros(len(ordered_files))
std_pos = np.zeros(len(ordered_files))

for i, file in enumerate(ordered_files):
    with open(file, 'r') as infile:
        for line in infile:
            if '@TYPE xy' in line:
                break
        names = ['time', 'pos']
        df = pd.read_csv(infile, delim_whitespace=True, names=names)
        mean_pos[i] = df['pos'].mean()
        std_pos[i] = df['pos'].std()


plt.errorbar(mean_pos, mean_force, xerr=None, yerr=std_force)
plt.title(r'$\langle \frac{dA}{dx} \rangle$')
plt.xlabel('distance from surface [nm]')
plt.ylabel('restraint force [kJ/mol/nm]')
plt.show()

free_energy = simps(mean_force, mean_pos)
print(free_energy)
