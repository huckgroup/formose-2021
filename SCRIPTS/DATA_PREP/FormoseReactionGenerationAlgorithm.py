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
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes
from NorthNet import network_generation as n_gen

from helpers import chem_info as info_params

'''Get reaction components'''
# there may be an error like:
# 'mapped atoms in the reactants were not mapped in the products.
# unmapped numbers are: 5'
# the error is attributable to the Cannizzaro reaction SMARTS, in which
# O:5 is not mapped to the products
reactions = {}
for r in info_params.reaction_SMARTS:
    SMARTS = info_params.reaction_SMARTS[r]
    split_rxn_SMARTS = SMARTS.split('>>')
    reactant_SMARTS = split_rxn_SMARTS[0].split('.')
    product_SMARTS = split_rxn_SMARTS[1].split('.')
    reactions[r] = Classes.ReactionTemplate(r, SMARTS,
                                            reactant_SMARTS,
                                            product_SMARTS)

C_patt = Chem.MolFromSmarts("[C]")
count_carbons = lambda x: x.GetSubstructMatches(C_patt)

'''Name'''
generation_protocol_string = "C6_network_c2_c3_reactions_symmetric"

'''Boundary conditions'''
# iterations  overshoot for C6, but do so to get
# all reaction paths and compounds possible up
# to C6 compounds
iterations = 6
start_smiles = ['O','[OH-]','O=CCO', 'C=O', 'OC=CO','OC1COC(CO1)O']

initiator_species = [Classes.Compound(x) for x in start_smiles]
t_net = Classes.Network([],
                "{}".format(generation_protocol_string),
                generation_protocol_string)

t_net.add_compounds(initiator_species)

'''Reactivity constructor'''
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

reaction_pattern = deprotonation_rules + protonation_rules + \
                    aldol_addition_rules + sugar_aldol_rules + \
                    others

'''Expansion operation'''
x = 0
while x < iterations:
    for task in reaction_pattern:
        n_gen.extend_network_task(t_net, reactions[task])

    n_gen.carbonyl_migration_isomers_multiclass(t_net,
            deprot_rules = [reactions[d] for d in deprotonation_rules],
            prot_rules = [reactions[p] for p in protonation_rules]
            )

    '''Removing compounds > C6'''
    # In effect, this is equivalent to setting all chain-growing reaction
    # rules to not occur for C6 compounds
    # i.e. [$(C(O)=CO)!$(C(O)=C(O)C(O)C(O)C(O)CO)], etc.
    remove_compounds = [t_net.NetworkCompounds[c] for c in t_net.NetworkCompounds
        if len(count_carbons(t_net.NetworkCompounds[c].Mol)) > 6]
    t_net.remove_compounds(remove_compounds)

    x+=1

'''Saving the results'''
with open(repository_dir/'FORMOSE_REACTION/FullFormoseReaction.txt', 'w') as f:
    for r in t_net.NetworkReactions:
        f.write('{}\n'.format(r))

with open(repository_dir/'FORMOSE_REACTION/FormoseReactionNetwork.pickle','wb') as f:
    pickle.dump(t_net, f)
