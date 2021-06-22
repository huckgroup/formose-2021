import os
import pickle
from pathlib import Path

from NorthNet import Classes
from NorthNet import info_params

import __init__
from helpers.network_load_helper import load_reaction_list
from helpers.network_load_helper import convert_to_networkx
from helpers.network_plotting import plot_network

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from helpers import load_series

base_directory = Path('/Users/williamrobinson/Documents/Nijmegen')
# base_directory = Path(r'C:\Users\willi\Documents')
fname = base_directory/'safestore_DEP/Series_info.csv'
series_dict = load_series.load_series_sequences(fname)

'''Import clusters for ordering the data sets'''
clusters = {}
with open('information_sources/clusters.txt', 'r') as f:
	for c,line in enumerate(f):
		ins = line.strip('\n').split(',')
		clusters[c] = ins[1:]

with open('information_sources/FormoseReactionNetwork.pickle', 'rb') as f:
	FormoseNetwork = pickle.load(f)

directory = Path('information_sources/network_lists')
loaded_reactions = {}
for file in os.listdir(directory):
	loaded_reactions[file.split('_')[0]] = []
	with open(directory/file, 'r') as f:
		for line in f:
			ins = line.strip('\n')
			loaded_reactions[file.split('_')[0]].extend(ins.split(','))

reactions = {}
for c in clusters:
	for e in clusters[c]:
		if e in loaded_reactions:
			if len(loaded_reactions[e]) > 0:
				reactions[e] = loaded_reactions[e]

all_reactions = []
for r in reactions:
	addition = list(set(reactions[r]))
	all_reactions.extend(addition)

all_reactions = list(set(all_reactions))
all_reactions.sort()

full_set_reaction_classes = []
for r in FormoseNetwork.NetworkReactions:
	full_set_reaction_classes.append(FormoseNetwork.get_reaction_name(r))

all_reaction_classes = list(set(full_set_reaction_classes))
all_reaction_classes.sort()

base_reaction_class_numbers = [full_set_reaction_classes.count(x)
								for x in all_reaction_classes]

reaction_expression_stack = np.zeros((len(reactions),len(all_reaction_classes)))
for c,r in enumerate(reactions):
	reaction_classes = [FormoseNetwork.get_reaction_name(x) for x in reactions[r]]
	for c2,rc in enumerate(all_reaction_classes):
		reaction_expression_stack[c,c2] = reaction_classes.count(rc)

remove_idx = np.argwhere(np.all(reaction_expression_stack[:,...] == 0, axis=0))

reaction_expression_stack = np.delete(reaction_expression_stack, remove_idx, axis = 1)
all_reaction_classes = [x for i,x in enumerate(all_reaction_classes) if i not in remove_idx]

colours = np.array([])
for x in range(0,len(reaction_expression_stack)):
	for y in range(0,len(reaction_expression_stack[0])):
		colours = np.hstack((colours, reaction_expression_stack[x,y]/reaction_expression_stack[x].max()))

reaction_class_names = [info_params.reaction_class_names[r] for r in all_reaction_classes]

from mpl_toolkits.axes_grid1 import make_axes_locatable
fig, ax = plt.subplots(figsize = (35/2.54,15/2.54))

ax.tick_params(axis = 'both', which = 'both', length = 0)
# ax.set_position([0.25,0.25,0.7,0.7])
# mappable = ax.scatter(x_vals,y_vals, s=15, c = colours, zorder = 3, marker = 's', cmap = cm.seismic)
im = ax.imshow(reaction_expression_stack.T, cmap = 'cividis')
ax.set_xticks(np.arange(0,len(reactions),1))
ax.set_xticklabels([*reactions], rotation = 90, fontsize = 6)

ax.set_yticks(np.arange(0,len(reaction_class_names),1))
ax.set_yticklabels(reaction_class_names, fontsize = 6)

# ax.set_position([0.1,0.1,0.8,0.8])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('Fractional expression', rotation=270, labelpad = 10)
cbar.ax.tick_params(labelsize= 6)
# ax.set_aspect('equal')
fig.tight_layout()
plt.savefig('plots_for_paper/reaction_expression_heatmap.png', dpi = 600)
