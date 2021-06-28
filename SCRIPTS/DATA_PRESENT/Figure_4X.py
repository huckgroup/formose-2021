import sys
import pickle
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from networkx.drawing.nx_agraph import graphviz_layout

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes
from helpers import chem_info as info_params
from helpers.network_load_helper import convert_to_networkx
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')
network_file = repository_dir/'FORMOSE_REACTION/FullFormoseReaction.txt'

# load in experiment information.
exp_info = pd.read_csv(exp_info_dir, index_col = 0)
# sequences of data set keys
series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)

# loading in the formose reaction as a NorthNet Network Object
formose_file = repository_dir/'FORMOSE_REACTION/FormoseReactionNetwork.pickle'
with open(formose_file, 'rb') as f:
	FormoseNetwork = pickle.load(f)

series_sel = 'Formaldehyde_paper_series'
condition_sel_x = '[C=O]/ M'
factor = 1000
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

# create a network merging all of the networks
F = nx.DiGraph()
for n in networks:
	F = nx.compose(F,networks[n])

# create a layout for F (this will be the layout for each plotted network)
pos = graphviz_layout(F, prog = 'neato')

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
			networks[n].nodes[node]['color'] = "#00000"
		else:
			networks[n].nodes[node]['color'] = info_params.colour_assignments[node]
	for edge in networks[n].edges:
		for e in edge:
			if '>>' in e:
				r_class = FormoseNetwork.get_reaction_name(e)
				col = info_params.reaction_class_colours[r_class]
				networks[n].edges[edge]['color'] = col

'''Plotting series in four panels'''
fig_width = 7.91/2.54 # cm conversion to inches for plt
fig_height = 8.42/2.54 # cm conversion to inches for plt

base_linew = 0.5

fig,ax = plt.subplots(nrows = 2, ncols = 2,
				figsize = (fig_width, fig_height))
axes = ax.flatten()

for c,n in enumerate(networks):
	axes[c].plot(base_net_plot[0],base_net_plot[1],
				c = '#acb5ad',
				linewidth = base_linew,
				zorder =0)

	for e in networks[n].edges:
	    arrow = FancyArrowPatch(networks[n].nodes[e[0]]['pos'],
	                            networks[n].nodes[e[1]]['pos'],
	                            arrowstyle='-|>',
	                            path = None,
	                            connectionstyle='Arc',
	                            zorder = 1,
	                            facecolor = networks[n].edges[e]['color'],
	                            edgecolor = networks[n].edges[e]['color'],
	                            linewidth = 1,
	                            mutation_scale = 10,
	                            shrinkA = 0,
	                            shrinkB = 0,
	                            alpha = 1)
	    axes[c].add_patch(arrow)

	axes[c].set_axis_off()

plt.show()
