import sys
import pickle
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes
from helpers.network_load_helper import convert_to_networkx

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
	pass
