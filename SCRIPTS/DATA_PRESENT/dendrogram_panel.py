import numpy as np
from pathlib import Path
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from NorthNet import info_params
from NorthNet.calculations import calculations
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

import __init__
from helpers import load_series
from helpers.network_load_helper import load_from_edge_list
from helpers.network_load_helper import load_coordinates_list
from helpers.loading_helper import exp_info, modified_averages, prime_header

select_conditions = ['[C=O]/ M', '[CaCl2]/ M', '[NaOH]/ M','[O=C(CO)CO]/ M',
                    'residence_time/ s', 'temperature/ oC']
condition_names = ['[formaldehyde]/ mM', '[CaCl$_2$]/ mM',
                   '[NaOH]/ mM','[dihydroxyacetone]/ mM', 'Residence time/ s',
                   'Temperature/ $^\circ$C']
condition_dict = {}
for s in select_conditions:
    condition_dict[s] = {}
    for e in exp_info:
        condition_dict[s][e] = exp_info[e].parameters[s]

edge_list = 'information_sources/dendrogram_edgelist.csv'
coord_list = 'information_sources/dendrogram_coordinates.csv'
G = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G,pos)
c_ops.normalise_network_coordinates(G)
c_ops.rotate_network(G, -np.pi/2)
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)


fig, ax = plt.subplots(ncols = 1+len(select_conditions),
                        figsize = (13.44/2.54, 17.64/2.54))

for a in ax:
    a.plot(lines[0],lines[1],
           c = '#000000', linewidth = 2.5,
           zorder = 0)

subplot_width = 0.06
subplot_height = 0.06
for n in G.nodes:
    if n in exp_info:
        xy = G.nodes[n]['pos']
        L = xy[0] - subplot_width/2
        B = xy[1] - subplot_height/2
        W = subplot_width
        H = subplot_height
        axin = ax[0].inset_axes([L,B,W,H], transform=ax[0].transData)

        axin.pie(modified_averages[n]/np.amax(modified_averages[n]),
                 colors = [info_params.colour_assignments[x.split('/')[0]]
                 for x in prime_header])

factor = 1000
for c,v in enumerate(condition_dict,1):
    leaf_colours = []
    leaf_x = []
    leaf_y = []
    for n in G.nodes:
        if n in exp_info:
            xy = G.nodes[n]['pos']
            leaf_x.append(xy[0])
            leaf_y.append(xy[1])
            if 'residence' in v or 'temperature' in v:
                leaf_colours.append(condition_dict[v][n])
            else:
                leaf_colours.append(condition_dict[v][n]*factor)

    # ax[c].set_title(v)
    scattr = ax[c].scatter(leaf_x,leaf_y, s = 10, c = leaf_colours,
                            cmap = cm.coolwarm,
                            zorder = 1)

    new_cbar = fig.add_axes([0.025 + ((c-1)*(0.95/len(select_conditions))), 0.05,
                             1/(len(select_conditions)+1), 0.25])
    new_cbar.set_axis_off()
    cbar = plt.colorbar(scattr,
                        ax = new_cbar,
                        location="bottom")
    cbar.set_label(condition_names[c-1], fontsize = 6)

for a in ax:
    a.set_axis_off()
    ylm = a.get_ylim()
    a.set_ylim(ylm[1],ylm[0])

ax[0].set_position([0.1, 0.375, 0.85, 0.575])
for x in range(1,len(ax)):
    ax[x].set_position([0.025 + ((x-1)*(0.95/len(select_conditions))), 0.1,
                    1/(len(select_conditions)+1), 0.25])

plt.savefig('plots_for_paper/dendrogram_panel.png', dpi = 600)
plt.close()
