import os
import sys
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes

import helpers.chem_info as info_params
from helpers.network_plotting import plot_network
from helpers.network_load_helper import load_reaction_list
from helpers.network_load_helper import convert_to_networkx

# name for output files
figname = 'Figure3X'
# set paths to files
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

# load in experiment information.
exp_info = pd.read_csv(exp_info_dir, index_col = 0)
series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)

# Import clusters for ordering the data sets
clusters = {}
with open(repository_dir/'RESOURCES/clusters.txt', 'r') as f:
	for c,line in enumerate(f):
		ins = line.strip('\n').split(',')
		clusters[c] = ins[1:]

# loading in the formose reaction as a NorthNet Network Object
with open(repository_dir/'FORMOSE_REACTION/FormoseReactionNetwork.pickle', 'rb') as f:
	FormoseNetwork = pickle.load(f)

# get a list of all of the reactions classes
# present in the search framework
reaction_classes = {}
for r in FormoseNetwork.NetworkReactions:
	class_name = FormoseNetwork.get_reaction_name(r)
	if class_name in reaction_classes:
		reaction_classes[class_name] += 1
	else:
		reaction_classes[class_name] = 1

class_names = [*reaction_classes]
class_names.sort()

# load in the reaction lists determined for
# modulated data sets.
# use dictionary insertion ordering to
# add network reactions into a
networks = []
for e in exp_info.index:
	if exp_info.loc[e,'Modulated_component'] != 'None':
		fname = '{}_reaction_list.txt'.format(e)
		with open(reaction_list_directory/fname, 'r') as f:
			for line in f:
				lines = f.readlines()
		rxns = []
		for l in lines:
			rxns.append(FormoseNetwork.NetworkReactions[l.strip('\n')])

		networks.append(Classes.Network(rxns, e, ''))

merged_network = Classes.Network([],'merged','')
for n in networks:
	for r in n.NetworkReactions:
		merged_network.add_reactions([n.get_reaction(r)])

observed_reaction_classes = {cls:0 for cls in class_names} # container for all observed reaction classes
for r in merged_network.NetworkReactions:
	class_name = merged_network.get_reaction_name(r)

	observed_reaction_classes[class_name] += 1

# get number of each reaction class in the networks
# store the results in an array
reaction_numbers = np.zeros((len(networks),len(reaction_classes)))
for c,n in enumerate(networks):
	for r in n.NetworkReactions:
		cls = n.get_reaction_name(r)
		idx = class_names.index(cls)
		reaction_numbers[c,idx] += 1

# normalise the reaction counts to the total
# number of the reactions of that class present
# in the search framework.

# this measure brings out much more detail
# in the data that other measures.
# it can be viewed as a measure on how the
# environmental conditions have 'sculpted'
# the observed pathways from the set
# for c,r in enumerate(reaction_classes):
# 	reaction_numbers[:,c] /= reaction_classes[r]


# normalise to scores to the total number of
# reaction classes observed in the networks
# how the reaction classes chang relative
# to their union
for c,r in enumerate(observed_reaction_classes):
	reaction_numbers[:,c] /= observed_reaction_classes[r]

# remove column containing zeroes or nan
reaction_numbers = np.nan_to_num(reaction_numbers)
zero_idx = np.argwhere(np.all(reaction_numbers[..., :] == 0, axis=0))
reaction_numbers = np.delete(reaction_numbers,zero_idx, axis = 1)
# update the class names
class_names = [c for i,c in enumerate(class_names) if i not in zero_idx]
# get experiment labels
exp_labels = [n.Name for n in networks]

# plot the results in a heatmap
# cividis is a good colourmap for this purpose
# as it is easier for those with
# colour vision deficiency to view compared to other
# colour maps, and it cycles between low values of blue
# and higher values of yellow
# it looks nice, to me, too.
# see Nuñez, PLoS, 2018
from mpl_toolkits.axes_grid1 import make_axes_locatable
fig, ax = plt.subplots(figsize = (35/2.54,15/2.54))

ax.tick_params(axis = 'both', which = 'both', length = 0)
# ax.set_position([0.25,0.25,0.7,0.7])
# mappable = ax.scatter(x_vals,y_vals, s=15, c = colours, zorder = 3, marker = 's', cmap = cm.seismic)
im = ax.imshow(reaction_numbers.T, cmap = 'cividis')
ax.set_xticks(np.arange(0,len(exp_labels),1))
ax.set_xticklabels(exp_labels, rotation = 90, fontsize = 6)

ax.set_yticks(np.arange(0,len(class_names),1))
ax.set_yticklabels(class_names, fontsize = 6)

# ax.set_position([0.1,0.1,0.8,0.8])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('Fractional expression', rotation=270, labelpad = 10)
cbar.ax.tick_params(labelsize= 6)
# ax.set_aspect('equal')
fig.tight_layout()
plt.savefig(repository_dir/'PLOTS/{}.png'.format(figname), dpi = 600)
