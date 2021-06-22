import networkx as nx

def reactants_product_edges(reaction):
    reactants = reaction.split('>>')[0].split('.')
    products = reaction.split('>>')[1].split('.')

    edges = []
    for r in reactants:
        edges.append((r,reaction))
    for p in products:
        edges.append((reaction,p))

    # G.add_edges_from([(0, 1), (1, 2)])
    return edges

def shortest_path_between_compounds(G, source, target):

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
