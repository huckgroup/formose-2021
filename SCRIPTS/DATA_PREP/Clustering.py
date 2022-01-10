'''
Performs a clustering analysis to provide partitioned subsets of the data.
'''
import sys
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)

from NorthNet.network_visualisation import coordinates as c_ops

from helpers.network_load_helper import load_from_edge_list
from helpers.network_load_helper import load_coordinates_list

# get the repository directory for file output
repository_dir = Path(__file__).parents[2]
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

exp_info = pd.read_csv(exp_info_dir, index_col = 0)

experiment_names = list(exp_info.index)

average_data = pd.read_csv(
                            derived_parameters_dir/'AverageData.csv', 
                            index_col = 0
                            )
# remove empty columns
average_data = average_data.dropna(axis = 1)
# remove columns containing only zeros
average_data = average_data.loc[:, (average_data != 0).any(axis=0)]

amplitude_data = pd.read_csv(
                            derived_parameters_dir/'AmplitudeData.csv', 
                            index_col = 0
                            )
# remove empty columns
amplitude_data = amplitude_data.dropna(axis = 1)
# remove columns containing only zeros
amplitude_data = amplitude_data.loc[:, (amplitude_data != 0).any(axis=0)]

data = average_data.to_numpy()
amplitudes = amplitude_data.to_numpy()

augmented_amplitude_matrix = np.hstack((data, amplitudes))
# Calculate pairwise euclidean distances in the data
augmented_amplitude_distances = pdist(augmented_amplitude_matrix, 'correlation')
distances_squareform = squareform(augmented_amplitude_distances)
# the linkage function should detect that a distance matrix is being passed to it.
augmented_amplitude_linkages = linkage(
                                        augmented_amplitude_distances, 
                                        method='average',
                                        metric='', 
                                        optimal_ordering=False
                                        )

cut_level = 0.12
cluster_labels = fcluster(augmented_amplitude_linkages,
                          criterion = 'distance',
                          t = cut_level)
cluster_labels -= 1
n_clusters = len(np.unique(cluster_labels))

edge_list = repository_dir/'RESOURCES/dendrogram_edgelist.csv'
coord_list = repository_dir/'RESOURCES/dendrogram_coordinates.csv'
G3 = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G3, pos)
c_ops.normalise_network_coordinates(G3)
# c_ops.rotate_network(G, -np.pi/4 - np.pi/6)
lines = c_ops.get_network_lineplot(G3)
dots  = c_ops.get_network_scatter(G3)

kmeans = KMeans(n_clusters = 7)
kmeans.fit(distances_squareform)
cluster_labels = kmeans.labels_
n_clusters = len(np.unique(cluster_labels))
print('number of clusters', n_clusters)

fig, ax = plt.subplots(ncols = 2, figsize = (10,5))
ax[1].axhline(y = cut_level)
dendro = dendrogram(
                    augmented_amplitude_linkages, 
                    ax = ax[1], 
                    color_threshold = cut_level
                    )

ax[0].plot(
            lines[0], lines[1],
            c = '#000000', 
            linewidth = 2.5,
            zorder = 0
            )

for x in range(n_clusters):
    idx = np.where(cluster_labels == x)[0]
    names = [v for i,v in enumerate(exp_info.index) if i in idx]
    x_coords = []
    y_coords = []
    for n in names:
        if n in G3.nodes:
            xy = G3.nodes[n]['pos']
            x_coords.append(xy[0])
            y_coords.append(xy[1])

    ax[0].scatter(x_coords, y_coords, s = 10)

    ax[0].annotate(x+1, xy = (np.average(x_coords), np.average(y_coords)),
                    fontsize = 12, fontweight = 'bold')
# the y-axis is flipped upside down
# not necessary, but I forgot to remove it
# before making the final figures, so it
# gets left in!
ylm = ax[0].get_ylim()
ax[0].set_ylim(ylm[1],ylm[0])
ax[0].set_aspect('equal')
ax[0].set_axis_off()
plt.savefig(repository_dir/'RESOURCES/KMeans_cluster_positions.png', dpi = 600)

clusters = []
for x in range(n_clusters):
    clust = []
    idx = np.where(cluster_labels == x)[0]
    for i in idx:
        clust.append(experiment_names[i])
    clusters.append(clust)

with open(repository_dir/'RESOURCES/clusters.txt', 'w') as f:
    for c,v in enumerate(clusters):
        f.write('Cluster_{},'.format(c))
        [f.write('{},'.format(x)) for x in v]
        f.write('\n')

exp_list = [experiment_names[x] for x in dendro['leaves']]
with open(repository_dir/'RESOURCES/leaf_list.txt', 'w') as f:
    [f.write('{},'.format(x)) for x in exp_list]
