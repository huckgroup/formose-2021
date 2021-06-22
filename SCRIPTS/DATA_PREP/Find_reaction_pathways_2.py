import numpy as np
import networkx as nx
from rdkit import Chem
from pathlib import Path
import matplotlib.pyplot as plt

from NorthNet.file_loads import data_loads, info_loads
from NorthNet import info_params

import __init__
from helpers.network_load_helper import load_reaction_list,convert_to_networkx
from helpers.network_load_helper import load_network_from_reaction_list
from helpers.loading_helper import carbon_inputs
from helpers import pathway_helpers as path_hlp

# base_directory = Path(r'C:\Users\willi\Documents')
base_directory = Path('/Users/williamrobinson/documents/nijmegen')
report_directory = base_directory/'safestore_DEP'

header = [x+'/ M' for x in info_params.smiles_to_names]
exp_info = info_loads.import_Experiment_information(report_directory/"Experiment_parameters.csv")
experiment_averages = data_loads.load_exp_compound_file('information_sources/AverageData.csv', header)
experiment_amplitudes = data_loads.load_exp_compound_file('information_sources/AmplitudeData.csv', header)

modifications = {x:[] for x in experiment_averages}
carbon_inputs = {x:[] for x in experiment_averages}
for v in experiment_averages:
	for p in exp_info[v].parameters:
		col_name = p
		if 'temperature' in p:
			continue
		if 'C' in p and not 'Ca' in p and exp_info[v].parameters[p] > 0.0:
			tag = p.split('/')[0][1:-1] + '/ M'
			if tag in header:
				i = header.index(tag)
				modifications[v].append(i)
				carbon_inputs[v].append(p.split('/')[0][1:-1])


scale_factor = 2e5
amplitude_filter = 1e-5
# enolates are the only species with C=C bonds
# considered within the scope of this study.
enol_patt = Chem.MolFromSmarts('C=C')

network_file = 'information_sources/FullFormoseReaction.txt'
formose_reactions = load_reaction_list(network_file)
formose_network = load_network_from_reaction_list(formose_reactions)

# conversion of formose reaction network to networkx
# DiGraph for further analysis
network = convert_to_networkx(formose_network)

# Define 'reagent' type nodes from the network
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

network_results = {}
null_networks = {}
for e in experiment_amplitudes:
	if np.sum(experiment_amplitudes[e]) == 0.0:
		print('no amps', e)
		continue

	# remove C7 compounds
	for c,p in enumerate(header):
		if p.count('C') == 7:
			experiment_amplitudes[e][c] = 0.0
			experiment_averages[e][c] = 0.0

	# sort amplitudes in decreasing magnitude
	idx = np.argsort(experiment_amplitudes[e])
	idx = np.flip(idx)
	sorted_amplitudes = experiment_amplitudes[e][idx]
	loc_header = [header[i] for i in idx]

	# filter amplitudes
	nonzero_idx_1 = np.argwhere(sorted_amplitudes > amplitude_filter).T[0]
	# loc_header becomes the ordered search list
	loc_header = [loc_header[i].split('/')[0] for i in nonzero_idx_1]

	# place the input reactants at the tope of the search list
	for ci in carbon_inputs[e]:
		if ci in loc_header:
			loc_header.remove(ci)
		else:
			pass

	# Get a list of all detected compounds
	nonzero_idx_2 = np.argwhere(experiment_averages[e] > 0.0).T[0]
	detected_compounds = [header[i].split('/')[0] for i in nonzero_idx_2]

	# include compounds which were not detected experimentally.
	# (if they could not be detected, they could still be present)
	detected_compounds.extend(undetected)

	# create a copy of the network to edit for pathway
	# searching.
	F = network.copy()

	# remove reactions from the networks for which the reactants
	# are not detected experimentally.
	# remove all reactions in which the inputs are products
	remove_reactions = []
	for r in formose_network.NetworkReactions:
		reactants = formose_network.NetworkReactions[r].Reactants
		products = formose_network.NetworkReactions[r].Products
		if not all(r in detected_compounds for r in reactants):
			remove_reactions.append(r)
		for inp in carbon_inputs[e]:
			if inp in products:
				remove_reactions.append(r)

	# convert remove_reactions to a set to remove duplicate entries
	remove_reactions = set(remove_reactions)
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

	# check through the graph and check that all products have edges
	# leading to them. If they do not, cycle backwards from them
	# in the sorted list until a connection is found
	loose_ends = []
	for n in R.nodes:
		if len(R.in_edges(n)) == 0 and n not in carbon_inputs[e]\
		 		and n not in remove_nodes:
			loose_ends.append(n)

	for n in loose_ends:
		if n in loc_header:
			reversed_list = list(reversed(loc_header[:loc_header.index(n)]))
			for r in reversed_list:
				edge_list = path_hlp.shortest_path_between_compounds(F, r, n)
				if edge_list:
					R.add_edges_from(edge_list)

		else:
			for r in carbon_inputs[e]:
				if r != 'C=O':
					edge_list = path_hlp.shortest_path_between_compounds(F, r, n)
					if edge_list:
						R.add_edges_from(edge_list)

	for rn in R.nodes:
		if '>>' in rn:
			R.nodes[rn]['size'] = scale_factor*5e-5
		elif rn in loc_header or rn in carbon_inputs[e]:
			R.nodes[rn]['size'] = experiment_amplitudes[e][header.index(rn +'/ M')]*scale_factor
		else:
			R.nodes[rn]['size'] = scale_factor*1e-4

	# remove nodes which clutter the layout
	[R.remove_node(x) for x in remove_nodes if x in R.nodes]

	network_results[e] = R.copy()

# plotting the results
from helpers.network_plotting import plot_network
for n in network_results:
	if len(network_results[n].nodes) < 3:
		continue
	else:
		print(n, len(network_results[n].nodes), len(network_results[n].edges))
		fname = 'plots_for_paper/network_plots/{}_network.png'.format(n)
		if len(network_results[n].edges) == 0:
			pass
		else:
			plot_network(network_results[n], fname, prog = 'neato')

# writing the resulting reaction lists
# to files
for n in network_results:
	output_list = []
	for node in network_results[n].nodes:
		if '>>' in node:
			output_list.append(node)

	with open('information_sources/network_lists/{}_reaction_list.txt'.format(n), 'w') as f:
		for o in output_list:
			f.write(o)
			f.write('\n')