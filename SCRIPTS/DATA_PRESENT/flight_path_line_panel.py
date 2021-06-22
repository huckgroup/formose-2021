import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

from NorthNet import info_params
from NorthNet.calculations import calculations
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

import __init__
from helpers import load_series
from helpers.loading_helper import exp_info, data, header, modified_averages
from helpers.network_load_helper import load_from_edge_list,load_coordinates_list

edge_list = 'information_sources/dendrogram_edgelist.csv'
coord_list = 'information_sources/dendrogram_coordinates.csv'
G = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G,pos)
c_ops.normalise_network_coordinates(G)
c_ops.rotate_network(G, -np.pi/2)
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)

base_directory = Path('/Users/williamrobinson/Documents/Nijmegen')
# base_directory = Path(r'C:\Users\willi\Documents')

fname = base_directory/'safestore_DEP/Series_info.csv'
series_dict = load_series.load_series_sequences(fname)

series_sel = 'Formaldehyde_2_series'
condition_sel = '[C=O]/ M'
x_factor = 1000
y_factor = 1000

series_x_values = [x_factor*exp_info[x].parameters[condition_sel]
                        for x in series_dict[series_sel]]

idx = [[*exp_info].index(x) for x in series_dict[series_sel]]

series_stack = np.zeros((len(idx),len(data[0])))

for x in range(0,len(idx)):
    series_stack[x] = data[idx[x]]

series_progression = series_stack.T

node_path = [(a,b) for a,b in zip(series_dict[series_sel],series_dict[series_sel][1:])]

fig, ax = plt.subplots(figsize=(8.965/2.54,6.55/2.54))

ax.set_xlabel('[Formaldehyde]/ mM')
ax.set_ylabel('concentration/ mM')

ax.plot(lines[0],lines[1], zorder = 0, c = '#000000')
ax.scatter(dots[0],dots[1], zorder = 0, c = '#000000', s = 5)

for node_pair in node_path:
    arrow = FancyArrowPatch(G.nodes[node_pair[0]]['pos'],
                            G.nodes[node_pair[1]]['pos'],
                            arrowstyle='-|>',
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
plt.savefig('plots_for_paper/{}_flight_path_panel.png'.format(series_sel), dpi = 600)
plt.close()
