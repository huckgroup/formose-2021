'''
Processing to derive a measure of the relative reaction expression of reaction
classes found in data-derived reaction networks.
'''
import os
import sys
import pickle
import numpy as np
import pandas as pd
from pathlib import Path

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

def getReactionClassNames(network):
	reactionClassNames = []
	for r in network.NetworkReactions:
		rxn_name = network.get_reaction_name(r)
		reactionClassNames.append(rxn_name)
	return reactionClassNames

def loadNetwork(filename):
	with open(filename, 'r') as f:
		for line in f:
			lines = f.readlines()
	reactions = []
	for l in lines:
		reactions.append(FormoseNetwork.NetworkReactions[l.strip('\n')])

	network = Classes.Network(reactions, e, '')

	return network

# set paths to files
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
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
formose_reaction_classes = getReactionClassNames(FormoseNetwork)

class_names = list(set(formose_reaction_classes))
class_names.sort()

reaction_classes = {c:formose_reaction_classes.count(c) for c in class_names}

# load in the reaction lists determined for
# modulated data sets.
# use dictionary insertion ordering to
# add network reactions into a
networks = []
for e in exp_info.index:
	if exp_info.loc[e,'Modulated_component'] != 'None':
		fname = '{}_reaction_list.txt'.format(e)

		network = loadNetwork(reaction_list_directory/fname)

		networks.append(network)

merged_network = Classes.Network([],'merged','')
for n in networks:
	for r in n.NetworkReactions:
		merged_network.add_reactions([n.get_reaction(r)])

# container for all observed reaction classes
experiment_r_classes = getReactionClassNames(merged_network)

exp_unique_r_classes = list(set(experiment_r_classes))
exp_unique_r_classes.sort()

observed_reaction_classes = {cls:experiment_r_classes.count(cls) for cls in exp_unique_r_classes}

# get number of each reaction class in the networks
# store the results in an array
reaction_numbers = np.zeros((len(networks),len(exp_unique_r_classes)))
for c,n in enumerate(networks):
	for r in n.NetworkReactions:
		rxn_name = n.get_reaction_name(r)
		idx = exp_unique_r_classes.index(rxn_name)
		reaction_numbers[c,idx] += 1

# normalise the the reaction class scores to the total number of reaction
# classes of the same type observed for all of the reaction systems analysed
reaction_numbers_normalised = np.zeros(reaction_numbers.shape)
for c,r in enumerate(class_names):
	reaction_numbers_normalised[:,c] = reaction_numbers[:,c]/observed_reaction_classes[r]
