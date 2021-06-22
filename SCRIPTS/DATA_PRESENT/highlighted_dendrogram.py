import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from networkx.drawing.nx_agraph import graphviz_layout

from scipy import cluster
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cophenet, fcluster, dendrogram

from NorthNet import info_params
from NorthNet import Classes
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops
from NorthNet.file_exports.plotting import network_plotting as northplot
from NorthNet.file_exports.plotting.misc import confidence_ellipse

import __init__
from helpers.loading_helper import exp_info
from helpers.loading_helper import data, amplitudes
from helpers.loading_helper import modified_averages, prime_header, clrs

from helpers.network_load_helper import load_reaction_list
from helpers.network_load_helper import convert_to_networkx

from helpers.cluster_tree import graph_from_linkage
from helpers.network_plotting import plot_network_in_axis

'''
Define groups of experiments to highlight
'''
groups = [['FRN090A', 'FRN093B', 'FRN093A','FRN090C','FRN093B','FRN060B','FRN060B',
            'FRN051A', 'FRN060B','FRN060A','FRN060A','FRN059B','FRN060A','FRN094C',
            'FRN060A','FRN094A'],
        ['FRN098A','FRN104B','FRN104B','FRN097F','FRN097F','FRN097K'],
        ['FRN098C','FRN103','FRN098C','FRN103','FRN097M','FRN097C','FRN098D',
        'FRN097C', 'FRN097I','FRN097E','FRN098H','FRN103'],
        ['FRN087A','FRN087E','FRN087C'],
        ['FRN071A','FRN071D', 'FRN071C'],
        ['FRN077A','FRN077B','FRN077C'],
        ['FRN099A','FRN099B','FRN099E'],
        ['FRN062C','FRN055B','FRN063E','FRN055B','FRN055B','FRN050A',
        'FRN061B','FRN055B','FRN055B','FRN054A','FRN054A','FRN054B',
        'FRN067C','FRN054A', 'FRN067A','FRN054A','FRN067B']

]

experiment_names = [*exp_info]
directory = 'information_sources/network_lists'

networks = {}
for f in os.listdir(directory):
    r_list = load_reaction_list('{}/{}'.format(directory,f))
    if len(r_list) == 0:
        pass
    else:
        rxns = [Classes.Reaction(r) for r in r_list]
        ntwrk = Classes.Network(rxns, f.split('_')[0],'')
        networks[f.split('_')[0]] = convert_to_networkx(ntwrk)

augmented_amplitude_matrix = np.hstack((data, amplitudes))
# Calculate pairwise euclidean distances in the data
augmented_amplitude_distances = pdist(augmented_amplitude_matrix, 'correlation')
# the linkage function should detect that a distance matrix is being passed to it.
augmented_amplitude_linkages = linkage(augmented_amplitude_distances, method='average',
                                metric='', optimal_ordering=False)
ci, coph_dists = cophenet(augmented_amplitude_linkages, augmented_amplitude_distances)
print('data + amplitude Cophenetic CC:', ci)


G3 = graph_from_linkage(augmented_amplitude_linkages, id_modifier = '')
pos = graphviz_layout(G3, prog = 'sfdp', args='-GK=0.7')

c_ops.set_network_coords(G3, pos)
c_ops.normalise_network_coordinates(G3)
net_coords = c_ops.get_network_scatter(G3)

fig,ax = plt.subplots(figsize = (5,5))
ax.scatter(net_coords[0],net_coords[1], c= 'k')
net_lines = c_ops.get_network_lineplot(G3)
ax.plot(net_lines[0], net_lines[1], c = '#000000', linewidth = 5)


subplot_width = 0.04
subplot_height = 0.04
for n in G3.nodes:
    if G3.nodes[n]['leaf']:
        pass
        # xy = G3.nodes[n]['pos']
        # L = xy[0] - subplot_width/2
        # B = xy[1] - subplot_height/2
        # W = subplot_width
        # H = subplot_height
        # axin = ax.inset_axes([L,B,W,H], transform=ax.transData)
        # if 'network' in n:
        #     idx = net_names[G3.nodes[n]['id']]
        # else:
        #     idx = experiment_names[G3.nodes[n]['id']]
        #
        # axin.pie(modified_averages[idx]/np.amax(modified_averages[idx]),
        #          colors = [info_params.colour_assignments[x.split('/')[0]]
        #          for x in prime_header])
ax.set_axis_off()

for g in groups:
    coords = []
    for exp in g:
        idx = experiment_names.index(exp)
        coords.append(G3.nodes[str(idx)]['pos'])
    coords = np.array(coords)

    xy = coords.T
    ax.plot(xy[0],xy[1], alpha = 1,
            c = 'w', linewidth = 20,
            solid_capstyle='round', zorder = 0,
            path_effects=[pe.Stroke(linewidth=25, foreground='#000000'),
            pe.Normal()])


plt.savefig('plots_for_paper/hierarchical_diagram_condition_highlights.png')
