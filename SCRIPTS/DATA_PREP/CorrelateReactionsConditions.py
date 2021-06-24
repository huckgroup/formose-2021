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

exp_names = [n.Name for n in networks]
# get number of each reaction class in the networks
# store the results in an array
reaction_numbers = np.zeros((len(networks),len(reaction_classes)))
for c,n in enumerate(networks):
	for r in n.NetworkReactions:
		cls = n.get_reaction_name(r)
		idx = class_names.index(cls)
		reaction_numbers[c,idx] += 1

# remove column containing zeroes or nan
zero_idx = np.argwhere(np.all(reaction_numbers[..., :] == 0, axis=0))
reaction_numbers = np.delete(reaction_numbers,zero_idx, axis = 1)
reaction_names = [r for i,r in enumerate(class_names) if i not in zero_idx]
# correlate the reaction counts against conditions
# only take the conditions for experiments
# for which reaction pathways have been determined
selected_exps = exp_info.loc[exp_names].iloc[:,2:-1]
# create numpy array from exp_info
conditions = selected_exps.to_numpy(dtype = np.float64)
# remove column containing zeroes or nan
zero_idx = np.argwhere(np.all(conditions[..., :] == 0, axis=0))
conditions = np.delete(conditions,zero_idx, axis = 1)
condition_names = [c for i,c in enumerate(selected_exps.columns)
                                                        if i not in zero_idx]
# transpose so conditions are row-wise
conditions = conditions.T

# transpose reaction_numbers so that reactions are row-wise
reaction_numbers = reaction_numbers.T

correlation_matrix = np.zeros((len(conditions), len(reaction_numbers)))

for x in range(0,len(conditions)):
    for y in range(0,len(reaction_numbers)):

        basis = np.vstack((conditions[x],reaction_numbers[y]))

        ans = np.corrcoef(basis)
        correlation_matrix[x,y] = ans[0,1]

from mpl_toolkits.axes_grid1 import make_axes_locatable
fig, ax = plt.subplots()

im = ax.imshow(correlation_matrix, cmap = 'seismic')
ax.set_xlabel('reactions')
ax.set_ylabel('conditions')
ax.set_xticks(np.arange(0,len(reaction_names),1))
ax.set_yticks(np.arange(0,len(condition_names),1))
ax.set_xticklabels(reaction_names, rotation = 90, fontsize = 6)
ax.set_yticklabels(condition_names, fontsize = 6)
ax.tick_params(axis = 'both', which = 'both', length = 0)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('Correlation coefficient',labelpad = 10)
cbar.ax.tick_params(labelsize= 6)
fig.tight_layout()
plt.show()
