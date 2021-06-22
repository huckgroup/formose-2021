import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

import helpers.chem_info as info_params
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops
from helpers.network_load_helper import load_from_edge_list,load_coordinates_list

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

exp_info = pd.read_csv(exp_info_dir, index_col = 0)

edge_list = repository_dir/'RESOURCES/dendrogram_edgelist.csv'
coord_list = repository_dir/'RESOURCES/dendrogram_coordinates.csv'
G = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G,pos)
c_ops.normalise_network_coordinates(G)
# c_ops.rotate_network(G, -np.pi/2)
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)

figname = 'Figure2F'
series_sel = 'Ca_OH_grid_series'
condition_sel_1 = '[CaCl2]/ M'
condition_sel_2 = '[NaOH]/ M'
x_factor = 1000
y_factor = 1000

series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)
data_keys = series_seqs.loc[series_sel]
data_set_selections = list(data_keys.dropna())

node_path = [(a,b) for a,b in zip(data_set_selections,data_set_selections[1:])]

fig, ax = plt.subplots(figsize=(8.965/2.54,6.55/2.54))

ax.plot(lines[0],lines[1], zorder = 0, c = '#000000')
ax.scatter(dots[0],dots[1], zorder = 0, c = '#000000', s = 5)

for node_pair in node_path:
    arrow = FancyArrowPatch(G.nodes[node_pair[0]]['pos'],
                            G.nodes[node_pair[1]]['pos'],
                            arrowstyle='-',
                            path = None,
                            connectionstyle='Angle3',#'Angle3'
                            zorder = 1,
                            facecolor = '#2E2EFE',
                            edgecolor = '#2E2EFE',
                            linewidth = 1,
                            mutation_scale = 10,
                            shrinkA = 0,
                            shrinkB = 0,
                            alpha = 1)
    ax.add_patch(arrow)

ax.set_axis_off()
ylm = ax.get_ylim()
ax.set_ylim(ylm[1],ylm[0])
fig.tight_layout()
plt.savefig(repository_dir/'PLOTS/{}.png'.format(figname), dpi = 600)
plt.close()
