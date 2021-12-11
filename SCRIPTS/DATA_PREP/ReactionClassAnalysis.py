'''
Processing to derive a measure of the relative reaction expression of reaction classes found in data-derived reaction networks.
'''
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
        lines = []
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

# normalise the the reaction class scores to the total number of reaction
# classes of the same type observed for all of the reaction systems analysed
reaction_numbers_normalised = np.zeros(reaction_numbers.shape)
for c,r in enumerate(class_names):
    reaction_numbers_normalised[:,c] = reaction_numbers[:,c]/observed_reaction_classes[r]

# get experiment labels
exp_names = [n.Name for n in networks]
exp_labels = [exp_info.loc[n.Name,'Experiment_entry'] for n in networks]

# organise zones of the reaction expression array by reaction expression
# load clusters upon which reaction ordering will be based.
with open(repository_dir/'RESOURCES/clusters.txt', 'r') as f:
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
        [f.write('{},'.format(reaction_numbers_r_order[c,x]))
                        for x in range(0,len(reaction_numbers_r_order[c]))]
        f.write('\n')

with open(repository_dir/'RESOURCES/reaction_expression_normalised.csv', 'w') as f:
    f.write(',')
    [f.write('{},'.format(rn)) for rn in class_names]
    f.write('\n')
    for c,v in enumerate(exp_names_r_order):
        f.write('{},'.format(v))
        [f.write('{},'.format(reaction_numbers_normalised_r_order[c,x]))
                        for x in range(0,len(reaction_numbers_normalised_r_order[c]))]
        f.write('\n')
