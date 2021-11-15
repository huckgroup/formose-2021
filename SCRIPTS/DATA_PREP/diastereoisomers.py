import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
repository_dir = Path(__file__).parents[2]

from helpers import chem_info 

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
exp_info = pd.read_csv(exp_info_dir, index_col = 0)

df = pd.read_csv(derived_parameters_dir/'AverageData.csv')
# remove empty columns
df = df.dropna(axis = 1)

from helpers.chem_info import compound_numbering

comp_rev = {compound_numbering[n]:n for n in compound_numbering}
exp_code = df.iloc[:,0]
compounds = df.columns[1:]
df = df.drop(df.columns[0], axis = 1)

data = df.to_numpy()
comp_by_num = {}
comp_by_name = {}
for f in range(0,len(data[0])):
    name = compounds[f].split('/')[0]
    if name in compound_numbering:
        comp_by_num[compound_numbering[name]] = data[:,f]
        comp_by_name[name] = data[:,f]

selections = ['14', '15']
markers = ['.', 'D', '^', 'x']
print(exp_info.index)
x_axis = [exp_info['[CaCl2]/ M'][e] for e in exp_code]
factor = 1000

fig, ax = plt.subplots()

for c,s in enumerate(selections):
    name = comp_rev[s]
    print(name)
    colour = chem_info.colour_assignments[name]
    ax.scatter(x_axis, factor*comp_by_num[s], c = colour, label = s, marker =
            markers[c]) 
ax.set_xlabel('[CaCl$_2$]/ M')
ax.set_ylabel('Concentration/ mM')
ax.legend()
plt.show()
