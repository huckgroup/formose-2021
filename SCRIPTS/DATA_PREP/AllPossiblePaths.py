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
			  'O=C[C@@H](O)CO','OC1COC(O)CO1']

# find and append enolates to the list of undetected compounds
# enolates are the only species with C=C bonds
# considered within the scope of this study.
enol_patt = Chem.MolFromSmarts('C=C')
for c in formose_network.NetworkCompounds:
	compound = formose_network.NetworkCompounds[c]
	if compound.Mol.HasSubstructMatch(enol_patt):
		undetected.append(c)

# loop over amplitude data sets
e = 'FRN093B'

# remove C7 compounds: including C7 compounds
# in the search framework creates a huge jump
# in the size of the network.
for c,p in enumerate(header):
    if p.count('C') == 7:
        experiment_amplitudes[e][c] = 0.0
        experiment_averages[e][c] = 0.0

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

search_list = [d for d in detected_compounds if d not in remove_nodes]

# container for results
R = nx.DiGraph()
[R.add_node(d) for d in search_list]
# iterate through the list, finding shortest
# paths between each compound
for s in range(0,len(search_list)):
    for t in range(0,len(search_list)):
        source = search_list[s]
        target = search_list[t]
        if source == target:
            continue
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

# copy the network into the result container.
network_result = R.copy()

# writing the resulting reaction lists to files
output_list = []
for node in network_result.nodes:
    if '>>' in node:
        output_list.append(node)

with open(repository_dir/'REACTION_LISTS/{}_undirected_search_reaction_list.txt'.format(e), 'w') as f:
    for o in output_list:
        f.write(o)
        f.write('\n')
