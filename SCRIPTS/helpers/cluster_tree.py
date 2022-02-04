'''
Functions for dealing with dendrogram data.
'''
import networkx as nx
from scipy.cluster import hierarchy

def graph_from_linkage(linkage_mat, id_modifier = ''):
    '''
    Create a networkx DiGraph from a scipy linkage matrix.

    linkage_mat: sparse array
        See
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
    id_modifier: str
        Modifier which will be added to node names.

    Returns
    -------
    G: networkx DiGraph
    '''
    rootnode, nodelist = hierarchy.to_tree(linkage_mat, rd= True)

    G = nx.DiGraph()

    for n in nodelist:
        source = '{}{}'.format(n.id, id_modifier)
        G.add_node(source, distance = n.dist, leaf = n.is_leaf(),
                    id = n.id)
        if n.count > 1:
            target_left = '{}{}'.format(n.left.id, id_modifier)

            G.add_node(target_left, distance = n.left.dist, leaf = n.left.is_leaf(),
                        id = n.left.id)
            G.add_edge(source, target_left, weight = abs(n.left.dist - n.dist),
                        dist = abs(n.left.dist - n.dist))

            target_right = '{}{}'.format(n.right.id, id_modifier)
            G.add_node(target_right, distance = n.right.dist, leaf = n.right.is_leaf(),
                        id = n.right.id)
            G.add_edge(source, target_right, weight = abs(n.right.dist - n.dist),
                        dist = abs(n.right.dist - n.dist))

    return G

def createNestedJSON(linkage_mat):
    '''
    Create a JSON string resembling the following:
    {
      name: "root",
      children: [
        {name: "child #1"},
        {
          name: "child #2",
          children: [
            {name: "grandchild #1"},
            {name: "grandchild #2"},
            {name: "grandchild #3"}
          ]
        }
      ]
    }

    Parameters
    ----------
    linkage_mat: scipy linkage matrix

    Returns
    -------
    json_string: str
        JSON string
    '''
    import json

    def connectNode(node, parent):
        new_entry = {"name":node.id, 'children':[]}
        parent["children"].append(new_entry)

        if node.left:
            connectNode(node.left, new_entry)
        if node.right:
            connectNode(node.right, new_entry)

    rootnode, nodelist = hierarchy.to_tree(linkage_mat, rd= True)

    root_name = rootnode.id

    json_build = {"name":root_name, 'children':[]}

    connectNode(rootnode.left, json_build)
    connectNode(rootnode.right, json_build)

    json_string = json.dumps(json_build)

    return json_string
