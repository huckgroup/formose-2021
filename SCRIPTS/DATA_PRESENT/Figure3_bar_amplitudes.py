'''
Create a plot of the shortest path in 
a reaction network vs. the amplitude of
the compound.
'''
import sys
import pickle
import pandas as pd
import networkx as nx
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes
from helpers import chem_info as info_params
from helpers.network_load_helper import convert_to_networkx

data_folder = repository_dir/'DATA'
determined_params_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

# Define experiment code
exp_code = 'FRN089C'
factor = 1000

data_file = determined_params_dir/'AmplitudeData.csv'
# Load in data as a guide for detected species
amplitude_data = pd.read_csv(data_file, index_col = 0)
amplitude_data = amplitude_data.dropna(axis = 1)
compounds = [x.split('/')[0] for x in amplitude_data.columns]
data_entry = amplitude_data.loc[exp_code].to_numpy()
detected_indices = np.argwhere(data_entry > 0.0)
detected_compounds = [compounds[i[0]] for i in detected_indices]
amplitude_series = data_entry[detected_indices]


# load compound numbering scheme
compound_numbering = {}
with open(repository_dir/'COMPOUND_INFO/compound_numbering.txt', 'r') as f:
    for line in f:
        ins = line.strip('\n').split(',')
        compound_numbering[ins[0]] = int(ins[1])

numbers_to_compounds = {compound_numbering[c]:c for c in compound_numbering}

sel_numbers = [3,9,12,13,18]
sel_compounds = [numbers_to_compounds[n] for n in sel_numbers]
amps = [factor*amplitude_data.loc[exp_code,c+'/ M'] for c in sel_compounds]
clrs = [info_params.colour_assignments[x] for x in sel_compounds]
compound_numbers = [compound_numbering[x] for x in sel_compounds]

width = 5/2.54
height = 5/2.54
print(amps)
fig, ax = plt.subplots(figsize = (width, height))

ax.bar(sel_compounds, amps, color = clrs)
ax.set_ylabel('Concentration/ mM')
ax.set_xticklabels(sel_numbers, fontweight = 'bold')
plt.show()
