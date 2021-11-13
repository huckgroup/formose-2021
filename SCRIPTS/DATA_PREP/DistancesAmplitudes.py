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

# load network from which reactions can be used
formose_file = repository_dir/'FORMOSE_REACTION'/'FormoseReactionNetwork.pickle'
with open(formose_file, 'rb') as f:
	FormoseNetwork = pickle.load(f)

# Load in reaction list and convert to
# a networkx graph
reaction_file = reaction_list_directory/f'FRN088B_reaction_list.txt'

with open(reaction_file, 'r') as f:
    lines = f.readlines()

reactions = [x.strip('\n') for x in lines]

rxn_objs = [FormoseNetwork.NetworkReactions[r] for r in reactions]
network = Classes.Network(rxn_objs, exp_code, '')
graph = convert_to_networkx(network)

# edit network to remove secondary reactants
node_removals = ['C=O', 'O', '[OH-]']
[graph.remove_node(node)
    for node in node_removals
        if node in graph.nodes]

# find shortest paths from input to detected species
search_list = [x for x in detected_compounds
                    if x not in node_removals]

# container for results
shortest_paths = {c:[] for c in detected_compounds}

root_node = 'O=C(CO)CO'
for d in search_list:
    if d == root_node:
        continue
    if d in graph.nodes:
        path = nx.shortest_path(graph, source = root_node, target = d)
        shortest_paths[d] = path
    else:
        print(d)

path_lengths = {d:len(shortest_paths[d]) for d in shortest_paths}

x_ax = [path_lengths[d] for d in path_lengths]
y_ax = [amplitude_data.loc[exp_code,c+'/ M'] for c in detected_compounds]
clrs = [info_params.colour_assignments[x] for x in detected_compounds]
compound_numbers = [compound_numbering[x] for x in detected_compounds]

width = 5/2.54
height = 5/2.54

fig, ax = plt.subplots(figsize = (width, height))

ax.bar(compound_numbers, y_ax, color = clrs)

plt.show()
