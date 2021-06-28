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
from helpers.network_load_helper import load_reaction_list
from helpers.network_load_helper import convert_to_networkx

# name for output files
figname = 'Figure3A'
# set paths to files
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

# load in experiment information.
exp_info = pd.read_csv(exp_info_dir, index_col = 0)

# loading in the formose reaction as a NorthNet Network Object
formose_file = repository_dir/'FORMOSE_REACTION/FormoseReactionNetwork.pickle'
with open(formose_file, 'rb') as f:
	FormoseNetwork = pickle.load(f)

# get a list of all of the reactions classes
# present in the search framework
reaction_classes = {}
for r in FormoseNetwork.NetworkReactions:
	rxn_name = FormoseNetwork.get_reaction_name(r)
	class_name = info_params.reaction_class_names[rxn_name]
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

# container for all observed reaction classes
observed_reaction_classes = {cls:0 for cls in class_names}
for r in merged_network.NetworkReactions:
	rxn_name = merged_network.get_reaction_name(r)
	class_name = info_params.reaction_class_names[rxn_name]
	observed_reaction_classes[class_name] += 1

# get number of each reaction class in the networks
# store the results in an array
reaction_numbers = np.zeros((len(networks),len(reaction_classes)))
for c,n in enumerate(networks):
	for r in n.NetworkReactions:
		rxn_name = n.get_reaction_name(r)
		cls = info_params.reaction_class_names[rxn_name]
		idx = class_names.index(cls)
		reaction_numbers[c,idx] += 1

# remove column containing zeroes or nan
reaction_numbers = np.nan_to_num(reaction_numbers)
zero_idx = np.argwhere(np.all(reaction_numbers[..., :] == 0, axis=0))
reaction_numbers = np.delete(reaction_numbers,zero_idx, axis = 1)
# update the class names
class_names = [c for i,c in enumerate(class_names) if i not in zero_idx]

# normalise the scores to the highest reaction class count for the experiment
reaction_numbers_normalised = np.zeros(reaction_numbers.shape)
for c,r in enumerate(class_names):
	reaction_numbers_normalised[:,c] = reaction_numbers[:,c]/reaction_numbers[:,c].max()

# get experiment labels
exp_names = [n.Name for n in networks]
exp_labels = [exp_info.loc[n.Name,'Experiment_entry'] for n in networks]

# organise zones of the heatmap by reaction expression
# load clusters upon which reaction ordering will be based.
with open('RESOURCES/clusters.txt', 'r') as f:
	lines = f.readlines()

# create a new experiment ordering based on clusters
cluster_order_exp_labels = []
for l in lines:
	line= l.strip('\n').split(',')
	cluster_order_exp_labels.extend([x for x in line[1:] if x in exp_names])

with open(repository_dir/'RESOURCES/leaf_list.txt', 'r') as f:
	lines = f.readlines()

cluster_order_exp_labels = lines[0].split(',')
cluster_order_exp_labels = [x for x in cluster_order_exp_labels
															if x in exp_names]

# get the indices which will sort the data according to
# the new experiment order.
idx = [exp_names.index(x) for x in cluster_order_exp_labels]

exp_labels_r_order = [exp_labels[i] for i in idx]
exp_names_r_order = [exp_names[i] for i in idx]

# use idx to re-order reaction_numbers
reaction_numbers_r_order = reaction_numbers[idx]
reaction_numbers_normalised_r_order = reaction_numbers_normalised[idx]

with open(repository_dir/'RESOURCES/reaction_expression.csv', 'w') as f:
	f.write(',')
	[f.write('{},'.format(rn)) for rn in class_names]
	f.write('\n')
	for c,v in enumerate(exp_names_r_order):
		f.write('{},'.format(v))
		[f.write('{},'.format(reaction_numbers_r_order[c,x])) for x in range(0,len(reaction_numbers_r_order[c]))]
		f.write('\n')

# plot the results in a heatmap
# cividis is a good colourmap for this purpose
# as it is easier for those with
# colour vision deficiency to view compared to other
# colour maps, and it cycles between low values of blue
# and higher values of yellow.
# It looks nice, to me, too.
# see Nuñez, Anderton, Renslow, PLoS One, 2018, 13, 1–14.
# - 'Color map poem' W. E. Robinson, 2021
from mpl_toolkits.axes_grid1 import make_axes_locatable
fig, ax = plt.subplots(figsize = (36/2.54,15/2.54))

ax.tick_params(axis = 'both', which = 'both', length = 0)
im = ax.imshow(reaction_numbers_normalised_r_order.T, cmap = 'cividis')
# 0.5 offset to centre the ticklabels
ax.set_xticks(np.arange(0.5,len(exp_labels)+0.5,1))
ax.set_xticklabels([])

ax.set_yticks(np.arange(0,len(class_names),1))
ax.set_yticklabels(class_names, fontsize = 8)

ax.set_xlabel('Experiment')
ax.set_ylabel('Reaction class')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('Fractional expression',labelpad = 10)
cbar.ax.tick_params(labelsize= 6)

fig.tight_layout()
plt.savefig(repository_dir/'PLOTS/{}.png'.format(figname), dpi = 600)
plt.close()
