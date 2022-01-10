'''
Generates an reaction network for the formose reaction based on reaction rules.
'''
import sys
import pickle
from rdkit import Chem
from pathlib import Path

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)

from NorthNet import Classes
from NorthNet import network_generation as n_gen

from helpers import chem_info as info_params

# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

# Get reaction components

# There may be an error such as:
# 'mapped atoms in the reactants were not mapped in the products.
# unmapped numbers are: 5'
# the error is attributable to the Cannizzaro reaction SMARTS, in which
# O:5 is not mapped to the products. Apologies for this issue. To the best of
# my knowledge, it does not affect the results of the program as the reactions
# generated using this reaction SMARTS source appear to be correct.
reaction_rules = {}
for r in info_params.reaction_SMARTS:
    SMARTS = info_params.reaction_SMARTS[r]
    split_rxn_SMARTS = SMARTS.split('>>')
    reactant_SMARTS = split_rxn_SMARTS[0].split('.')
    product_SMARTS = split_rxn_SMARTS[1].split('.')
    reaction_rules[r] = Classes.ReactionTemplate(r, SMARTS,
                                            reactant_SMARTS,
                                            product_SMARTS)

# The following two lines offer a method to count the number of carbon atoms in
# a molecule.
C_patt = Chem.MolFromSmarts("[C]")
count_carbons = lambda x: x.GetSubstructMatches(C_patt)

# Name
generation_protocol_string = "C6_network_c2_c3_reactions_symmetric"

# Boundary conditions
# iterations overshoot for C6, but do so to get all reaction paths and 
# compounds possible up to C6 compounds
iterations = 6
# The initial set of compounds present in the reaction network
start_smiles = ['O','[OH-]','O=CCO', 'C=O', 'OC=CO','OC1COC(CO1)O']
initiator_species = [Classes.Compound(x) for x in start_smiles]

# Initialise an empty reaction network.
t_net = Classes.Network([],
                "{}".format(generation_protocol_string),
                generation_protocol_string)
# Add the initial set of compounds to the network.
t_net.add_compounds(initiator_species)

# Reactivity constructor
# The following lists emphasise the divisions of reaction types.
# They can be edited to remove the action of a particular reaction class upon
# the generated network. The names correspond to keys in the dictionary
# `reaction_rules`.
deprotonation_rules = ['deprotonation-t1','deprotonation-t2',
                       'deprotonation-t3','deprotonation-t4',
                       'deprotonation-t5','deprotonation-t6',
                       'deprotonation-t7','deprotonation-t8']

protonation_rules = ['protonation-t1-@','protonation-t1-@@',
                     'protonation-t2-@','protonation-t2-@@',
                     'protonation-t3-@','protonation-t3-@@',
                     'protonation-t4-@','protonation-t4-@@',
                     'protonation-t5-@','protonation-t5-@@',
                     'protonation-t6-@','protonation-t6-@@',
                     'protonation-t7-@','protonation-t7-@@',
                     'protonation-t8-@','protonation-t8-@@']

aldol_addition_rules = ['aldol_addition_fast_pt_C=O-t1-@',
                        'aldol_addition_fast_pt_C=O-t1-@@',
                        'aldol_addition_fast_pt_C=O-t2-@',
                        'aldol_addition_fast_pt_C=O-t2-@@',
                        'aldol_addition_fast_pt_C=O-t3-@',
                        'aldol_addition_fast_pt_C=O-t3-@@',
                        'aldol_addition_fast_pt_C=O-t4-@',
                        'aldol_addition_fast_pt_C=O-t4-@@',
                        'aldol_addition_fast_pt_C=O-t5-@',
                        'aldol_addition_fast_pt_C=O-t5-@@',
                        'aldol_addition_fast_pt_C=O-t6-@',
                        'aldol_addition_fast_pt_C=O-t6-@@',
                        'aldol_addition_fast_pt_C=O-t7-@',
                        'aldol_addition_fast_pt_C=O-t7-@@',
                        'aldol_addition_fast_pt_C=O-t8-@',
                        'aldol_addition_fast_pt_C=O-t8-@@']

sugar_aldol_rules = ['aldol_sugar_addition_fast_pt-t1-@@-@',
                     'aldol_sugar_addition_fast_pt-t2-@@-@',
                     'aldol_sugar_addition_fast_pt-t1-@@-@@',
                     'aldol_sugar_addition_fast_pt-t2-@@-@@',
                     'aldol_sugar_addition_fast_pt-t2-@-@',
                     'aldol_sugar_addition_fast_pt-t1-@-@',
                     'aldol_sugar_addition_fast_pt-t1-@-@@',
                     'aldol_sugar_addition_fast_pt-t2-@-@@']

others = ['retroaldol_to_enol-t1','retroaldol_to_enol-t2',
            'Cannizzaro', 'dimer_dissociation']

# The above reaction rule lists are combined into a single list.
reaction_pattern = deprotonation_rules + protonation_rules + \
                    aldol_addition_rules + sugar_aldol_rules + \
                    others

# Network generation operation
x = 0
while x < iterations:

    for task in reaction_pattern:
        n_gen.extend_network_task(t_net, reaction_rules[task])

    # Run the protonation/deprotonation rules repeatedly over the network
    # to generate all enol and epimer forms of the compounds present.
    n_gen.carbonyl_migration_isomers_multiclass(t_net,
            deprot_rules = [reaction_rules[d] for d in deprotonation_rules],
            prot_rules = [reaction_rules[p] for p in protonation_rules]
            )

    # Find and remove compounds containing more than 6 carbon atoms
    # In effect, this is equivalent to setting all chain-growing reaction
    # rules to not occur for C6 compounds
    # i.e. [$(C(O)=CO)!$(C(O)=C(O)C(O)C(O)C(O)CO)], etc.
    remove_compounds = []
    for c in t_net.NetworkCompounds:
        compound = t_net.NetworkCompounds[c] 
        if len(count_carbons(compound.Mol)) > 6:
            remove_compounds.append(compound)

    t_net.remove_compounds(remove_compounds)

    x+=1

# Saving the results
with open(repository_dir/'FORMOSE_REACTION/FullFormoseReaction.txt', 'w') as f:
    for r in t_net.NetworkReactions:
        f.write('{}\n'.format(r))

with open(repository_dir/'FORMOSE_REACTION/FormoseReactionNetwork.pickle','wb') as f:
    pickle.dump(t_net, f)
