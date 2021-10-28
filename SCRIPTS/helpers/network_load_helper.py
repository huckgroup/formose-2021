'''
Functions for loading and converting reaction network data.
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from NorthNet import Classes

def convert_to_networkx(network):
    '''
    Converts NorthNet network object to networkx object.

    Parameters
    ----------
    network: NorthNet Network object
        NorthNet Network to be converted to networkx network.
    save_images: Bool
        Whether to create images or not.
    Returns
    -------
    G: networkx DiGraph object
        Networkx version of the NorthNet network.
    '''
    import networkx as nx

    G = nx.DiGraph()

    for node in network.NetworkCompounds:
        G.add_node(node)

    for r in network.NetworkReactions:
        for sr in network.NetworkReactions[r].Reactants:
            for sp in network.NetworkReactions[r].Products:
                G.add_edge(sr,r)
                G.add_edge(r,sp)

    return G

def load_network_from_reaction_list(reaction_list):

    reactions = []
    for r in reaction_list:
        # Parse reaction SMILES
        split_reaction = r.split('>>')
        reactants = [Chem.MolFromSmiles(m) for m in split_reaction[0].split('.')]
        products = [Chem.MolFromSmiles(m) for m in split_reaction[1].split('.')]

        rxn = AllChem.ChemicalReaction() # Create an empty chemical reaction
        [rxn.AddReactantTemplate(r) for r in reactants]
        [rxn.AddProductTemplate(p) for p in products]
        reactions.append(Classes.Reaction(rxn))

    network = Classes.Network(reactions, '','')

    return network

def load_reaction_list(filename):
    from NorthNet import Classes

    reaction_list = []

    with open(filename, 'r') as f:
        for line in f:
            reaction_list.append(line.strip('\n'))

    return reaction_list

def load_from_edge_list(fname):
    import networkx as nx

    edge_list = []
    with open(fname, 'r') as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip('\n').split(',')
                edge_list.append(ins)

    G = nx.Graph()
    for e in edge_list:
        G.add_edge(e[0],e[1])

    return G

def load_coordinates_list(fname):
    pos = {}
    with open(fname, 'r') as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip('\n').split(',')
                pos[ins[0]] = [float(x) for x in ins[1:] if x != '']
    return pos
