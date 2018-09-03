import pandas as pd

names = 'nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB'
names = names.split()
df = pd.read_csv('../24SOX.itp', skiprows=0, skipfooter=-1, delim_whitespace=True, names=names, engine='python')
print(df)

sum_charge = df['charge'].sum()
print("Sum charge: " + str(sum_charge))
