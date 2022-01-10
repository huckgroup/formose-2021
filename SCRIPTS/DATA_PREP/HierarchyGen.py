'''
Creates a JSON file describing the hierarchical clustering of the ampltidue and
average compositional data.
'''
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)

from helpers import cluster_tree

# get the repository directory for file output
repository_dir = Path(__file__).parents[2]
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

exp_info = pd.read_csv(exp_info_dir, index_col = 0)

experiment_names = list(exp_info.index)

#############################################
# Load and pre-process compound averages data
#############################################
average_data = pd.read_csv(
                            derived_parameters_dir/'AverageData.csv', 
                            index_col = 0
                            )
# remove empty columns
average_data = average_data.dropna(axis = 1)
# remove columns containing only zeros
average_data = average_data.loc[:, (average_data != 0).any(axis=0)]

###############################################
# Load and pre-process compound amplitudes data
###############################################
amplitude_data = pd.read_csv(
                            derived_parameters_dir/'AmplitudeData.csv', 
                            index_col = 0
                            )
# remove empty columns
amplitude_data = amplitude_data.dropna(axis = 1)
# remove columns containing only zeros
amplitude_data = amplitude_data.loc[:, (amplitude_data != 0).any(axis=0)]

# Convert Pandas data frames to numpy arrays
data = average_data.to_numpy()
amplitudes = amplitude_data.to_numpy()

# Combine the averages and amplitudes arrays into a single array
# Clustering of the data is based on this array.
augmented_amplitude_matrix = np.hstack((data, amplitudes))
# Calculate pairwise distances between data entries
augmented_amplitude_distances = pdist(augmented_amplitude_matrix, 'correlation')
# the linkage function should detect that a distance matrix is being passed to it.
augmented_amplitude_linkages = linkage(
                                        augmented_amplitude_distances, 
                                        method='average',
                                        metric='', 
                                        optimal_ordering=False
                                        )

json_string = cluster_tree.createNestedJSON(augmented_amplitude_linkages)

with open(repository_dir/'RESOURCES/dendrogramHierarchy.json', 'w') as f:
    f.write(json_string)
