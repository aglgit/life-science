import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
from scipy.integrate import simps

files = ['force/%s' % f for f in os.listdir('force')]
ordered_files = sorted(files, key=lambda x: (int(re.sub('\D','',x)),x))

mean_force = np.zeros(len(ordered_files))

for i, file in enumerate(ordered_files):
    with open(file, 'r') as infile:
        for line in infile:
            if '@TYPE xy' in line:
                break
        names = ['time', 'force']
        df = pd.read_csv(infile, delim_whitespace=True, names=names)
        mean_force[i] = df['force'].mean()

files = ['pos/%s' % f for f in os.listdir('pos')]
ordered_files = sorted(files, key=lambda x: (int(re.sub('\D','',x)),x))

mean_pos = np.zeros(len(ordered_files))

for i, file in enumerate(ordered_files):
    with open(file, 'r') as infile:
        for line in infile:
            if '@TYPE xy' in line:
                break
        names = ['time', 'pos']
        df = pd.read_csv(infile, delim_whitespace=True, names=names)
        mean_pos[i] = df['pos'].mean()


plt.plot(mean_pos, mean_force)
plt.show()

free_energy = simps(mean_force, mean_pos)
print(free_energy)
