from NorthNet import Classes

def load_text_networks(path):

    load_starts = []
    dset_names = []
    with open(path, 'r') as f:
        for c,line in enumerate(f):
            if 'Reaction_entry' in line:
                load_starts.append(c)
                dset_names.append(line.strip('\n').split('_')[-1])

    reaction_networks = {dset_names[c]:[] for c,x in enumerate(load_starts)}
    with open(path, 'r') as f:
        for c,line in enumerate(f):
            if c in load_starts:
                load_key = dset_names[load_starts.index(c)]
            else:
                reaction_networks[load_key].append(line.strip('/n'))
    return reaction_networks

def network_from_reaction_list(reaction_list, network_name = "", description = "", reaction_mapping = False):
    from NorthNet import Classes
    from rdkit import Chem
    from NorthNet.reaction_operations.attributes import reconstitute_reaction

    header = ["Reaction", "Description", "Reaction ID","References"]
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    ID_dict = {S:str(c) for c,S in enumerate(alphabet)}
    network = Classes.Network([],network_name, description)
    add_to_net = []
    for r in reaction_list:

        reactants = [x for x in r.split(">>")[0].split(".") if x != ""]
        products = [x for x in r.split(">>")[1].split(".") if x != ""]

        r_mols = [Chem.MolFromSmiles(x) for x in reactants]
        p_mols = [Chem.MolFromSmiles(x) for x in products]

        r_canon = [Chem.MolToSmiles(x, isomericSmiles = True) for x in r_mols]
        p_canon = [Chem.MolToSmiles(x , isomericSmiles = True) for x in p_mols]

        r_string = reconstitute_reaction(r_canon,p_canon)

        rID = "".join([ID_dict[s] for s in r_string if s in ID_dict])

        line = [r_string, "User input reaction.", rID, "None"]

        add_to_net.append(Classes.Reaction( Classes.Reaction_Database_Entry(header, line) ))

    network.add_reactions(add_to_net)

    return network

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

    rxns = []
    for r in reaction_list:
        rxns.append(Classes.Reaction(r))

    network = Classes.Network(rxns, '','')

    return network

def load_reaction_list(filename):
    from NorthNet import Classes

    reaction_list = []

    with open(filename, 'r') as f:
        for line in f:
            reaction_list.append(line.strip('\n'))

    return reaction_list

def add_flow_terms(network, inputs):
    '''
    Add flow terms into network (inputs and outputs are
    empty species and )
    '''

    add_reactions = []
    for i in inputs:
        r_obj = Classes.Reaction(">>{}".format(i))
        r_obj.Classification = 'flow_input'
        add_reactions.append(r_obj)
    for c in network.NetworkCompounds:
        r_obj = Classes.Reaction("{}>>".format(c))
        r_obj.Classification = 'flow_output'
        add_reactions.append(r_obj)

    network.add_reactions(add_reactions)

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
