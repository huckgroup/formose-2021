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
                   'Temperature/ $^\circ$C', '[CaCl$_2$]/[NaOH]']
condition_dict = {}
for s in select_conditions:
    condition_dict[s] = {}
    for e in exp_info:
        condition_dict[s][e] = exp_info[e].parameters[s]

condition_dict['Ca_OH_ratio'] = {}
for e in exp_info:
    condition_dict['Ca_OH_ratio'][e] = exp_info[e].parameters['[CaCl2]/ M']/exp_info[e].parameters['[NaOH]/ M']/1000

edge_list = 'information_sources/dendrogram_edgelist.csv'
coord_list = 'information_sources/dendrogram_coordinates.csv'
G = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G,pos)
c_ops.normalise_network_coordinates(G)
c_ops.rotate_network(G, -np.pi/4 - np.pi/6)
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)

factor = 1000
for c,v in enumerate(condition_dict,0):
    fig, ax = plt.subplots(figsize = (3.13/2.54,2.81/2.54), frameon = False)
    ax.plot(lines[0],lines[1],
           c = '#000000', linewidth = 0.5,
           zorder = 0)
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
    scattr = ax.scatter(leaf_x,leaf_y, s = 2, c = leaf_colours,
                            cmap = cm.coolwarm,
                            zorder = 1)

    cbar = plt.colorbar(scattr,
                        ax = ax,
                        location="bottom")
    cbar.ax.tick_params(labelsize=6, length = 2, pad = 1)
    cbar.ax.set_position([0.0,0.25,0.95,0.05])
    ax.set_position([0.0,0.3,0.95,0.65])
    cbar.set_label(condition_names[c], fontsize = 6, labelpad = 2)
    ax.set_axis_off()
    ylm = ax.get_ylim()
    ax.set_ylim(ylm[1],ylm[0])
    plt.savefig('plots_for_paper/{}_dendrogam.png'.format(v.split('/')[0]), dpi = 600)
    plt.close()


# ax[0].set_position([0.1, 0.375, 0.85, 0.575])
# for x in range(1,len(ax)):
#     ax[x].set_position([0.025 + ((x-1)*(0.95/len(select_conditions))), 0.1,
#                     1/(len(select_conditions)+1), 0.25])
