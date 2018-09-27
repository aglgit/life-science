import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

systems = ['24SOX', '25-OX', '27-OX']
system_paths = ['../systems/SAM_%s/umbrella/clust-size.xvg' % s for s in systems]
names = ['Cluster #', '# Structures']

def read_xvg(filename, names):
    with open(filename, 'r') as infile:
        for line in infile:
            if '@g0 type bar' in line:
                break
        df = pd.read_csv(infile, delim_whitespace=True, names=names)

        return df

for i, path in enumerate(system_paths):
    df = read_xvg(path, names)
    total = df[names[1]].sum()
    p = df[names[1]].values / total
    entropy = np.sum(-p*np.log(p))

    plt.plot(df[names[0]], df[names[1]], label='%s, S=%.3f' % (systems[i], entropy))
plt.xlabel(names[0])
plt.ylabel(names[1])
plt.title("Cluster analysis of molecules at surface, $d_{SH} = 1.8$ nm")
plt.legend()
plt.show()
