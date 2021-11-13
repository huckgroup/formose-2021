'''
Get the path between lyxose and dihydroxyacetone
with only proton transfer and formaldehyde addition 
reactions.
'''

import sys
import pickle
import copy
from pathlib import Path
import networkx as nx

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes
from helpers import chem_info as info_params
from helpers.network_load_helper import convert_to_networkx

data_folder = repository_dir/'DATA'
determined_params_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

# Paths will be found between dihydroxyacetone and lyxose
dihydroxyacetone = info_params.canonical_SMILES['dihydroxyacetone']
lyxose = info_params.canonical_SMILES['lyxose']

# loading in the formose reaction as a NorthNet Network Object
formose_file = repository_dir/'FORMOSE_REACTION/FormoseReactionNetwork.pickle'
with open(formose_file, 'rb') as f:
	FormoseNetwork = pickle.load(f)

# edit the network to remove sugar sugar reactions
sugar_sugar_reactions = []
formaldehyde_reactions = []
for r in FormoseNetwork.NetworkReactions:
    template = FormoseNetwork.get_reaction_template(r).Name
    if 'sugar' in template:
        sugar_sugar_reactions.append(FormoseNetwork.NetworkReactions[r])
    elif 'C=O' in template:
        formaldehyde_reactions.append(FormoseNetwork.NetworkReactions[r])

# Create a network from which sugar-sugar reactions will be removed
Network1 = copy.deepcopy(FormoseNetwork)
    
Network1.remove_reactions(sugar_sugar_reactions)

G1 = Network1.convert_to_networkx()

# Create a network from which formaldehyde addition reactions will be removed
Network2 = copy.deepcopy(FormoseNetwork)
    
Network2.remove_reactions(formaldehyde_reactions)

G2 = Network2.convert_to_networkx()

# a list of compounds which will be removed for pathway searches
secondary_reactants = ['C=O', '[OH-]', 'O']
[G1.remove_node(n) for n in secondary_reactants]
[G2.remove_node(n) for n in secondary_reactants]

# find paths from dihydroxyacetone to lyxose with only formaldehyde addition
# chain growth
formaldehyde_addition_path = nx.shortest_path(G1, source = dihydroxyacetone, target = lyxose)

formaldehyde_addition_reactions = [x for x in formaldehyde_addition_path if '>>' in x]

# Find paths from dihydroxyacetone to lyxose with only sugar sugar addition
# reactions

sugar_addition_path = nx.shortest_path(G2, source = dihydroxyacetone, target = lyxose)

sugar_addition_reactions = [x for x in sugar_addition_path if '>>' in x]

# results
'''
print(formaldehyde_addition_reactions)
# Formaldehyde addition reactions
O=C(CO)CO.[OH-]>>OC=C(O)CO.[OH-]
C=O.OC=C(O)CO>>O=C(CO)[C@H](O)CO
O=C(CO)[C@H](O)CO.[OH-]>>OC=C(O)[C@H](O)CO.[OH-]
C=O.OC=C(O)[C@H](O)CO>>O=C([C@@H](O)CO)[C@H](O)CO
O=C([C@@H](O)CO)[C@H](O)CO.[OH-]>>OCC(O)=C(O)[C@H](O)CO.[OH-]
O.OCC(O)=C(O)[C@H](O)CO>>O.O=C(CO)[C@@H](O)[C@H](O)CO 
O=C(CO)[C@@H](O)[C@H](O)CO.[OH-]>>OC=C(O)[C@@H](O)[C@H](O)CO.[OH-]
O.OC=C(O)[C@@H](O)[C@H](O)CO>>O.O=C[C@@H](O)[C@@H](O)[C@H](O)CO
'''

'''
print(sugar_addition_reactions)
# Sugar sugar addition reactions 
O=C(CO)CO.[OH-]>>OC=C(O)CO.[OH-]
O.OC=C(O)CO>>O.O=C[C@H](O)CO
O=C[C@H](O)CO.OC=CO>>O=C[C@@H](O)[C@@H](O)[C@H](O)CO
'''
