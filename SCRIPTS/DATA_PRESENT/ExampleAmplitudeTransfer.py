'''
Script for creating elements 
for a figure panel explaining how 
input modulations are transferred 
to reaction products.
'''

import sys
import pickle
import pandas as pd
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
from helpers.layout import graphviz_layout
from NorthNet.network_visualisation import coordinates as c_ops
from matplotlib.patches import FancyArrowPatch

data_folder = repository_dir/'DATA'
determined_params_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

# load in experiment information.
exp_info = pd.read_csv(exp_info_dir, index_col = 0)
# sequences of data set keys
series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)
# average data for sizing nodes
average_data = pd.read_csv(determined_params_dir/'AverageData.csv', index_col = 0)
average_data = average_data.dropna(axis = 1)
# amplitude data for sizing nodes
amplitude_data = pd.read_csv(determined_params_dir/'AmplitudeData.csv', index_col = 0)
amplitude_data = amplitude_data.dropna(axis = 1)

# loading in the formose reaction as a NorthNet Network Object
formose_file = repository_dir/'FORMOSE_REACTION/FormoseReactionNetwork.pickle'
with open(formose_file, 'rb') as f:
	FormoseNetwork = pickle.load(f)

data_set = 'FRN093C'

# Get the compound names and amplitudes for the
# data set selected
compounds = amplitude_data.columns
amps = amplitude_data.loc[data_set].to_numpy()
averages = average_data.loc[data_set].to_numpy()

# Get the reaction list determined for the data set
directed_fname = f'{data_set}_reaction_list.txt'
with open(reaction_list_directory/directed_fname, 'r') as f:
    directed_lines = f.readlines()

# Convert the reaction list to a Network object
# before converting it to a networkx DiGraph
rxns = []
for l in directed_lines:
    rxns.append(FormoseNetwork.NetworkReactions[l.strip('\n')])

n_net = Classes.Network(rxns, data_set, '')

exp_network = convert_to_networkx(n_net)

# Get the reaction list determined 
# from an undirected search
# of the data set
undirected_fname = f'{data_set}_undirected_search_reaction_list.txt'
with open(reaction_list_directory/undirected_fname, 'r') as f:
    undirected_lines = f.readlines()

# Convert the reaction list to a Network object
# before converting it to a networkx DiGraph
rxns = []
for l in undirected_lines:
    rxns.append(FormoseNetwork.NetworkReactions[l.strip('\n')])

undirected_n_net = Classes.Network(rxns, data_set+'undirected', '')

undirected_network = convert_to_networkx(undirected_n_net)

# remove secondary reactant nodes
node_removals = ['C=O', 'O', '[OH-]']
[exp_network.remove_node(node) for node in node_removals
                                    if node in exp_network.nodes]

[undirected_network.remove_node(node) for node in node_removals
                                    if node in undirected_network.nodes]

# Generate a layout for the graph
# and replace coordinates
F = nx.compose(exp_network, undirected_network)
pos = graphviz_layout(F, render_engine = 'neato')
c_ops.set_network_coords(F,pos)
c_ops.normalise_network_coordinates(F)
pos = {n:F.nodes[n]['pos'] for n in F.nodes}
c_ops.set_network_coords(exp_network,pos)
c_ops.set_network_coords(undirected_network,pos)

# Create scatter plots with circles whose diameters are
# proportional to the amplitude of compounds

scale_factor = 1000000

fig, ax = plt.subplots(figsize=(5,5))

compound_scatter_x = []
compound_scatter_y = []
for n in undirected_network.nodes:
    if '>>' in n:
        pass
    else:
        compound_scatter_x.append(pos[n][0])
        compound_scatter_y.append(pos[n][1])

ax.scatter(compound_scatter_x, compound_scatter_y, s = 4, linewidth = 0.4, 
                    edgecolor = '#6a6b6c', facecolor = 'none')

for x in range(0,len(amps)):
    compound = compounds[x].split('/')[0]
    if amps[x] > 0.0 and compound in exp_network.nodes:
        radius = amps[x]*scale_factor
        colour = info_params.colour_assignments[compound]
        x = exp_network.nodes[compound]['pos'][0]
        y = exp_network.nodes[compound]['pos'][1]
        plt.scatter(x,y, c = colour, s = radius)

reaction_scatter_x = []
reaction_scatter_y = []
for n in undirected_network.nodes:
    if '>>' in n:
        reaction_scatter_x.append(pos[n][0])
        reaction_scatter_y.append(pos[n][1])

ax.scatter(reaction_scatter_x, reaction_scatter_y, marker = 'D', c = 'k', s = 2)

for e in exp_network.edges:
    arrow = FancyArrowPatch(exp_network.nodes[e[0]]['pos'],
                            exp_network.nodes[e[1]]['pos'],
                            arrowstyle='-|>',
                            path = None,
                            connectionstyle='Arc',
                            facecolor = 'k',
                            edgecolor = 'k',
                            linewidth = 1,
                            mutation_scale = 5,
                            shrinkA = 2,
                            shrinkB = 1,
                            alpha = 1,
                            zorder = 1)

    ax.add_patch(arrow)

for e in undirected_network.edges:
    arrow = FancyArrowPatch(undirected_network.nodes[e[0]]['pos'],
                            undirected_network.nodes[e[1]]['pos'],
                            arrowstyle='-|>',
                            path = None,
                            connectionstyle='Arc',
                            facecolor = '#b3b8bd',
                            edgecolor = '#b3b8bd',
                            linewidth = 1,
                            mutation_scale = 5,
                            shrinkA = 2,
                            shrinkB = 1,
                            alpha = 1,
                            zorder = 0)

    ax.add_patch(arrow)

ax.set_axis_off()
plt.savefig(plot_folder/f'{data_set}_amplitudes.svg')

