
import numpy as np
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cophenet, fcluster, dendrogram

from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

import __init__
from helpers.cluster_tree import graph_from_linkage
from helpers.loading_helper import data, amplitudes, conditions, exp_info

experiment_names = [*exp_info]
print(data.shape, len(experiment_names))
augmented_amplitude_matrix = np.hstack((data, amplitudes))
# Calculate pairwise euclidean distances in the data
augmented_amplitude_distances = pdist(augmented_amplitude_matrix, 'correlation')
# the linkage function should detect that a distance matrix is being passed to it.
augmented_amplitude_linkages = linkage(augmented_amplitude_distances, method='average',
                                metric='', optimal_ordering=False)
ci, coph_dists = cophenet(augmented_amplitude_linkages, augmented_amplitude_distances)

G = graph_from_linkage(augmented_amplitude_linkages, id_modifier = 'augmented_net')
pos = graphviz_layout(G, prog = 'neato', args='-GK=0.7')

c_ops.set_network_coords(G, pos)

with open('information_sources/dendrogram_coordinates.csv', 'w') as f:
    f.write('node,x,y\n')
    for n in G.nodes:
        if G.nodes[n]['leaf']:
            idx = experiment_names[G.nodes[n]['id']]
        else:
            idx = n
        f.write('{},{},{}\n'.format(idx,
                                    G.nodes[n]['pos'][0],
                                    G.nodes[n]['pos'][1]))

with open('information_sources/dendrogram_edgelist.csv', 'w') as f:
    f.write('source,target\n')
    for e in G.edges:
        if G.nodes[e[0]]['leaf']:
            source = experiment_names[G.nodes[e[0]]['id']]
        else:
            source = e[0]

        if G.nodes[e[1]]['leaf']:
            target = experiment_names[G.nodes[e[1]]['id']]
        else:
            target = e[1]
        f.write('{},{}\n'.format(source, target))
