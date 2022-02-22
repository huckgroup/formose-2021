'''
Functions for pathway searching.
'''
import networkx as nx

def shortest_path_between_compounds(G, source, target):
    '''
    Find the shortest pathway between the source and target in G.

    Parameters
    ----------
    G: networkx DiGraph
    source: str
        Source node in G.
    target: str
        Target node in G.

    Returns
    -------
    list or bool
    '''
    if nx.has_path(G,source,target):
        path = nx.shortest_path(G, source = source, target = target)
        # check for bimolecular reactions and make sure that
        # second reactants/products are included in the
        edge_list = [(a,b) for a,b in zip(path,path[1:])]
        for p in path:
            if '>>' in p:
                loose_ends = reactants_product_edges(p)
                edge_list.extend(loose_ends)

        return edge_list
    else:
        return False

def reactants_product_edges(reaction):
    '''
    Get the reactants and products of the reaction and create tuple pairs 
    (reactant, reaction) and (reaction, product).

    Parameters
    ----------
    reaction: str
        Reaction SMILES

    Returns
    -------
    edges: list of tuples
    '''

    reactants = reaction.split('>>')[0].split('.')
    products = reaction.split('>>')[1].split('.')

    edges = []
    for r in reactants:
        edges.append((r,reaction))
    for p in products:
        edges.append((reaction,p))

    return edges
