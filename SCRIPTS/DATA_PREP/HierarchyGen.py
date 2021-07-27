import sys
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cophenet, fcluster, dendrogram

from sklearn.cluster import KMeans

from NorthNet.network_visualisation import coordinates as c_ops

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
distances_squareform = squareform(augmented_amplitude_distances)
# the linkage function should detect that a distance matrix is being passed to it.
augmented_amplitude_linkages = linkage(augmented_amplitude_distances, method='average',
                                metric='', optimal_ordering=False)

from helpers import cluster_tree

json_string = cluster_tree.createNestedJSON(augmented_amplitude_linkages)

with open(repository_dir/'RESOURCES/dendrogramHierarchy.json', 'w') as f:
    f.write(json_string)
