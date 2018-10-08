import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

systems = ['24SOX', '25-OX', '27-OX']
system_paths = ['../systems/SAM_%s/umbrella/hist.xvg' % s for s in systems]
names = ['z [nm]', 'PMF [$kJ mol^{-1} nm^{-1}$]']

def read_xvg(filename, names):
    with open(filename, 'r') as infile:
        for line in infile:
            if '@TYPE xy' in line:
                break
        df = pd.read_csv(infile, delim_whitespace=True, names=names)

        return df

zeros = [-47.0, -5.0, -4.5]
for i, path in enumerate(system_paths):
    df = read_xvg(path, names)
    df[names[1]] -= zeros[i]
    plt.plot(df[names[0]], df[names[1]], label=systems[i])

plt.xlabel(names[0])
plt.ylabel(names[1])
plt.title('Potential of mean force')
plt.legend()
plt.savefig('images/free_energy.png')
plt.show()
