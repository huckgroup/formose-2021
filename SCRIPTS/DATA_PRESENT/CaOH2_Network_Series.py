'''
Example reaction networks determined from amplitude data. Specifically
illustrating the variation of reaction pathways as the formaldehyde
concentration is variated. Figure 2B.
'''
import sys
import pickle
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes
from helpers.layout import graphviz_layout
from helpers import chem_info as info_params
from helpers.network_load_helper import convert_to_networkx
from NorthNet.network_visualisation import coordinates as c_ops

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

series_sel = 'Ca_OH_grid_series'
file_name = 'CaOHGrid'

# get the experiment codes for the series
data_keys = series_seqs.loc[series_sel]
data_set_selections = list(data_keys.dropna())

# load in the reaction lists determined for
# modulated data sets.
# use dictionary insertion ordering to
# add network reactions into a
networks = {}
for e in exp_info.index:
    for d in data_set_selections:
        fname = '{}_reaction_list.txt'.format(d)
        with open(reaction_list_directory/fname, 'r') as f:
            for line in f:
                lines = f.readlines()
        rxns = []
        for l in lines:
            rxns.append(FormoseNetwork.NetworkReactions[l.strip('\n')])

        n_net = Classes.Network(rxns, e, '')

        networks[d] = convert_to_networkx(n_net)

# remove secondary reactant nodes
node_removals = ['C=O', 'O', '[OH-]']
for n in networks:
    [networks[n].remove_node(node) for node in node_removals
                                        if node in networks[n]]
# create a network merging all of the networks
F = nx.DiGraph()
for n in networks:
    F = nx.compose(F,networks[n])

# create a layout for F (this will be the layout for each plotted network)
pos = graphviz_layout(F, render_engine = 'neato')

# use F to process coordinate system
c_ops.set_network_coords(F,pos)
c_ops.normalise_network_coordinates(F)

# get line plto for merged network
base_net_plot = c_ops.get_network_lineplot(F)

# create new position container from F
pos_norm = {n:F.nodes[n]['pos'] for n in F}

for n in networks:
    c_ops.set_network_coords(networks[n],pos_norm)

'''Add colouring information into the networks'''
for n in networks:
    for node in networks[n].nodes:
        if '>>' in node:
            networks[n].nodes[node]['color'] = "#000000"
        else:
            networks[n].nodes[node]['color'] = info_params.colour_assignments[node]

for n in networks:
    for edge in networks[n].edges:
        for e in edge:
            if '>>' in e:
                col = info_params.reaction_colours[e]
                networks[n].edges[edge]['color'] = col

'''Add sizing information into networks'''
min_node_size = 5
node_size_factor = 1e4
for n in networks:
    for node in networks[n].nodes:
        if '>>' in node:
            networks[n].nodes[node]['size'] = min_node_size
        elif node + '/ M' in average_data.columns:
            conc = average_data.loc[n,node +'/ M']
            amp = amplitude_data.loc[n,node +'/ M']
            networks[n].nodes[node]['size'] = conc*node_size_factor
        else:
            networks[n].nodes[node]['size'] = min_node_size

'''Plotting series in four panels'''
fig_width = 20/2.54 # cm conversion to inches for plt
fig_height = 10/2.54 # cm conversion to inches for plt

base_linew = 0.5

fig,ax = plt.subplots(nrows = 2, ncols = 4,
                figsize = (fig_width, fig_height))
axes = ax.flatten()

for c,n in enumerate(networks):
    axes[c].set_title(n)
    axes[c].plot(base_net_plot[0],base_net_plot[1],
                c = '#acb5ad',
                linewidth = base_linew,
                zorder = 0)

    for e in networks[n].edges:
        arrow = FancyArrowPatch(networks[n].nodes[e[0]]['pos'],
                                networks[n].nodes[e[1]]['pos'],
                                arrowstyle='-|>',
                                path = None,
                                connectionstyle='Arc',
                                facecolor = networks[n].edges[e]['color'],
                                edgecolor = networks[n].edges[e]['color'],
                                linewidth = 1,
                                mutation_scale = 5,
                                shrinkA = 2,
                                shrinkB = 1,
                                alpha = 1,
                                zorder = 1)

        axes[c].add_patch(arrow)

    # build node scatter
    compound_nodes_x = []
    compound_nodes_y = []
    compound_node_colours = []
    compound_node_sizes = []

    reaction_nodes_x = []
    reaction_nodes_y = []

    for node in networks[n].nodes:
        if '>>' in node:
            reaction_nodes_x.append(networks[n].nodes[node]['pos'][0])
            reaction_nodes_y.append(networks[n].nodes[node]['pos'][1])
        else:
            compound_nodes_x.append(networks[n].nodes[node]['pos'][0])
            compound_nodes_y.append(networks[n].nodes[node]['pos'][1])
            compound_node_colours.append(networks[n].nodes[node]['color'])
            compound_node_sizes.append(networks[n].nodes[node]['size'])

    # plot solid scatter for compound concentrations
    axes[c].scatter(compound_nodes_x, compound_nodes_y,
                facecolors = compound_node_colours,
                s = compound_node_sizes,
                zorder = 2,
                edgecolors = 'None',
                alpha = 1)

    axes[c].scatter(reaction_nodes_x, reaction_nodes_y,
                    c = '#000000',
                    s = min_node_size,
                    marker = 'D',
                    edgecolors = 'None',
                    zorder = 2,
                    alpha = 1)

    axes[c].set_axis_off()

    # optional annotations
    # for node in networks[n].nodes:
    #   if node in info_params.compound_numbering:
    #       number = info_params.compound_numbering[node]
    #       axes[c].annotate(number, xy = networks[n].nodes[node]['pos'],
    #                       ha = 'center', va = 'center',
    #                       fontsize = 6)

fig.tight_layout()
plt.show()
quit()
plt.savefig(plot_folder/'{}.png'.format(file_name), dpi = 600)
plt.savefig(plot_folder/'{}.svg'.format(file_name))
