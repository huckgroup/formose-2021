import sys
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cophenet, fcluster, dendrogram

from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

from helpers.cluster_tree import graph_from_linkage

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

exp_info = pd.read_csv(exp_info_dir, index_col = 0)

experiment_names = list(exp_info.index)

average_data = pd.read_csv(derived_parameters_dir/'AverageData.csv', index_col = 0)
# remove empty columns
average_data = average_data.dropna(axis = 1)
# remove columns containing only zeros
average_data = average_data.loc[:, (average_data != 0).any(axis=0)]

amplitude_data = pd.read_csv(derived_parameters_dir/'AmplitudeData.csv', index_col = 0)
# remove empty columns
amplitude_data = amplitude_data.dropna(axis = 1)
# remove columns containing only zeros
amplitude_data = amplitude_data.loc[:, (amplitude_data != 0).any(axis=0)]

data = average_data.to_numpy()
amplitudes = amplitude_data.to_numpy()

augmented_amplitude_matrix = np.hstack((data, amplitudes))
# Calculate pairwise euclidean distances in the data
augmented_amplitude_distances = pdist(augmented_amplitude_matrix, 'correlation')
# the linkage function should detect that a distance matrix is being passed to it.
augmented_amplitude_linkages = linkage(augmented_amplitude_distances, method='average',
                                metric='', optimal_ordering=False)

# get cophenetic correlation coefficient
ci, coph_dists = cophenet(augmented_amplitude_linkages, augmented_amplitude_distances)

G = graph_from_linkage(augmented_amplitude_linkages)
pos = graphviz_layout(G, prog = 'neato', args='-GK=0.7')

# add coordinates into nodes of G
for n in G.nodes:
    G.nodes[n]['pos'] = pos[n]

# write to file
with open(repository_dir/'RESOURCES/dendrogram_coordinates.csv', 'w') as f:
    f.write('node,x,y\n')
    for n in G.nodes:
        if G.nodes[n]['leaf']:
            idx = experiment_names[G.nodes[n]['id']]
        else:
            idx = n
        f.write('{},{},{}\n'.format(idx,
                                    G.nodes[n]['pos'][0],
                                    G.nodes[n]['pos'][1]))

with open(repository_dir/'RESOURCES/dendrogram_edgelist.csv', 'w') as f:
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
