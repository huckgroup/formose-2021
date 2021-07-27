'''
Searches for reaction pathways within an overarching formose reaction network
based on the transfer of modulated input concentrations to reaction products.
'''
import sys
import numpy as np
import pandas as pd
import networkx as nx
from rdkit import Chem
from pathlib import Path

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from helpers import chem_info as info_params
from helpers import pathway_helpers as path_hlp
from helpers.loading_helper import get_carbon_inputs
from helpers.network_load_helper import load_network_from_reaction_list
from helpers.network_load_helper import load_reaction_list,convert_to_networkx

header = [x+'/ M' for x in info_params.smiles_to_names]

# create Path objects for various information sources
data_dir = repository_dir/'DATA'
report_directory = data_dir/'DATA_REPORTS'
determined_params_dir = data_dir/'DERIVED_PARAMETERS'
exp_info_dir = repository_dir/'EXPERIMENT_INFO'
network_file = repository_dir/'FORMOSE_REACTION/FullFormoseReaction.txt'

exp_info = pd.read_csv(exp_info_dir/'Experiment_parameters.csv', index_col = 0)

experiment_names = list(exp_info.index)

average_data = pd.read_csv(determined_params_dir/'AverageData.csv', index_col = 0)
# remove empty columns
average_data = average_data.dropna(axis = 1)

amplitude_data = pd.read_csv(determined_params_dir/'AmplitudeData.csv', index_col = 0)
# remove empty columns
amplitude_data = amplitude_data.dropna(axis = 1)

experiment_amplitudes = {i:amplitude_data.loc[i].to_numpy()
												for i in amplitude_data.index}
experiment_averages = {i:average_data.loc[i].to_numpy()
												for i in average_data.index}

carbon_inputs = get_carbon_inputs(exp_info, average_data.columns)
for c in carbon_inputs:
	carbon_inputs[c] = [x.strip('/ M') for x in carbon_inputs[c]]

# Compounds below amplitude_filter will not be included in the search
amplitude_filter = 0.0 # / M
# define factor to scale nodes by for plotting
scale_factor = 2e5

# enolates are the only species with C=C bonds
# considered within the scope of this study.
enol_patt = Chem.MolFromSmarts('C=C')

# read the formose reaction network to be used as a
# basis for searches.
formose_reactions = load_reaction_list(network_file)
formose_network = load_network_from_reaction_list(formose_reactions)

# conversion of formose reaction network to networkx
# DiGraph for further analysis
network = convert_to_networkx(formose_network)

# Define 'reagent' type nodes from the network
# remove them for the network so they do not
# cause 'short circuiting'
remove_nodes = ['C=O', 'O', '[OH-]']
[network.remove_node(n) for n in remove_nodes]

# define compounds which were not detected experimentally
# (with consistency)
undetected = ['O', '[OH-]','O=CCO','O=C[C@H](O)CO',
			  'O=C[C@@H](O)CO','CC1COC(C)CO1']

# find and append enolates to the list of undetected compounds
# enolates are the only species with C=C bonds
# considered within the scope of this study.
enol_patt = Chem.MolFromSmarts('C=C')
for c in formose_network.NetworkCompounds:
	compound = formose_network.NetworkCompounds[c]
	if compound.Mol.HasSubstructMatch(enol_patt):
		undetected.append(c)

# create a container for the pathway search results
network_results = {}

# loop over amplitude data sets
for e in experiment_amplitudes:
	# if all amplitudes are equal to zero,
	# ingnore go to the next data set.
	if np.sum(experiment_amplitudes[e]) == 0.0:
		print('no amps', e)
		continue

	# remove C7 compounds: including C7 compounds
	# in the search framework creates a huge jump
	# in the size of the network.
	for c,p in enumerate(header):
		if p.count('C') == 7:
			experiment_amplitudes[e][c] = 0.0
			experiment_averages[e][c] = 0.0

	# sort amplitudes in decreasing magnitude
	idx = np.argsort(experiment_amplitudes[e])
	idx = np.flip(idx)
	sorted_amplitudes = experiment_amplitudes[e][idx]

	# create a similarly sorted list of compounds
	loc_header = [header[i] for i in idx]

	# filter amplitudes
	nonzero_idx_1 = np.argwhere(sorted_amplitudes > amplitude_filter).T[0]
	# loc_header becomes the ordered search list
	loc_header = [loc_header[i].split('/')[0] for i in nonzero_idx_1]

	# place the input reactants at the top of the search list
	for ci in carbon_inputs[e]:
		if ci in loc_header: loc_header.remove(ci)
		else: pass

	# Get a list of all detected compounds from the data vector
	nonzero_idx_2 = np.argwhere(experiment_averages[e] > 0.0).T[0]
	detected_compounds = [header[i].split('/')[0] for i in nonzero_idx_2]

	# include compounds which were not detected experimentally in the
	# 'detected' list
	# (if they could not be detected, their presence cannot be ruled out)
	detected_compounds.extend(undetected)

	# create a copy of the network to edit for pathway
	# searching.
	F = network.copy()

	# editing the network
	remove_reactions = []
	for r in formose_network.NetworkReactions:
		reactants = formose_network.NetworkReactions[r].Reactants
		products = formose_network.NetworkReactions[r].Products

		# remove reactions from the networks for which the reactants
		# are not detected experimentally.
		if not all(r in detected_compounds for r in reactants):
			remove_reactions.append(r)
		# remove all reactions in which the inputs are products
		for inp in carbon_inputs[e]:
			if inp in products:
				remove_reactions.append(r)

	# convert remove_reactions to a set to remove duplicate entries
	remove_reactions = set(remove_reactions)
	# remove the contents of remove_reactions from the network
	for r in remove_reactions:
		F.remove_node(r)

	# container for results
	R = nx.DiGraph()
	[R.add_node(header[i].split('/')[0]) for i in nonzero_idx_2]

	# first, find paths from inputs to compounds downstream
	for ci in carbon_inputs[e]:
		if ci == 'C=O':
			pass
		else:
			R.add_node(ci)
			for l in loc_header:
				edge_list = path_hlp.shortest_path_between_compounds(F, ci, l)
				if edge_list:
					R.add_edges_from(edge_list)

	# iterate down the amplitude list, find shortest pathways to every compound
	# downstream of the node.
	for s in range(0,len(loc_header)-1):
		for t in range(s+1,len(loc_header)-1):
			source = loc_header[s]
			target = loc_header[t]
			edge_list = path_hlp.shortest_path_between_compounds(F, source, target)
			if edge_list:
				R.add_edges_from(edge_list)

	# iterate through the graph and check that all products have edges
	# leading to them. If they do not, cycle backwards from them
	# in the sorted list until a connection is found

	# find products without input edges
	loose_ends = []
	for n in R.nodes:
		if len(R.in_edges(n)) == 0 and n not in carbon_inputs[e]\
		 		and n not in remove_nodes:
			loose_ends.append(n)

	# iterate over these loose ends and add connections
	for n in loose_ends:
		# if the compound is detected, try to connect to compounds upstream
		# from it
		if n in loc_header:
			reversed_list = list(reversed(loc_header[:loc_header.index(n)]))
			for r in reversed_list:
				edge_list = path_hlp.shortest_path_between_compounds(F, r, n)
				if edge_list:
					R.add_edges_from(edge_list)
		# otherwise, connect to inputs
		else:
			for r in carbon_inputs[e]:
				if r != 'C=O':
					edge_list = path_hlp.shortest_path_between_compounds(F, r, n)
					if edge_list:
						R.add_edges_from(edge_list)

	# insert some information for plotting
	for rn in R.nodes:
		if '>>' in rn:
			R.nodes[rn]['size'] = scale_factor*5e-5
		elif rn in loc_header or rn in carbon_inputs[e]:
			R.nodes[rn]['size'] = experiment_amplitudes[e][header.index(rn +'/ M')]*scale_factor
		else:
			R.nodes[rn]['size'] = scale_factor*1e-4

	# remove 'reagent nodes' (see remove_nodes above this loop)
	[R.remove_node(x) for x in remove_nodes if x in R.nodes]

	# copy the network into the result container.
	network_results[e] = R.copy()

# plotting the results
from helpers.network_plotting import plot_network
for n in network_results:
	if len(network_results[n].nodes) < 3:
		continue
	else:
		print(n, len(network_results[n].nodes), len(network_results[n].edges))
		fname = repository_dir/'PLOTS/network_plots/{}_network.png'.format(n)
		if len(network_results[n].edges) == 0:
			pass
		else:
			plot_network(network_results[n], fname, prog = 'neato')

# writing the resulting reaction lists to files
for n in network_results:
	output_list = []
	for node in network_results[n].nodes:
		if '>>' in node:
			output_list.append(node)

	with open(repository_dir/'REACTION_LISTS/{}_reaction_list.txt'.format(n), 'w') as f:
		for o in output_list:
			f.write(o)
			f.write('\n')
