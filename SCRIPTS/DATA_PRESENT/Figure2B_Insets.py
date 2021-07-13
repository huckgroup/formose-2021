import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet.network_visualisation import coordinates as c_ops

import helpers.chem_info as info_params
from helpers.network_load_helper import load_from_edge_list
from helpers.network_load_helper import load_coordinates_list

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

exp_info = pd.read_csv(exp_info_dir, index_col = 0)

select_conditions = ['[C=O]/ M', '[CaCl2]/ M', '[NaOH]/ M','[O=C(CO)CO]/ M',
                    'residence_time/ s', 'temperature/ oC']
condition_names = ['[formaldehyde]/ mM', '[CaCl$_2$]/ mM',
                   '[NaOH]/ mM','[dihydroxyacetone]/ mM', 'Residence time/ s',
                   'Temperature/ $^\circ$C', '[CaCl$_2$]/[NaOH]']

condition_dict = {}
for s in select_conditions:
    condition_dict[s] = {}
    for e in exp_info.index:
        condition_dict[s][e] = exp_info.loc[e,s]

condition_dict['Ca_OH_ratio'] = {}
for e in exp_info.index:
    condition_dict['Ca_OH_ratio'][e] = exp_info.loc[e,'[CaCl2]/ M']/exp_info.loc[e,'[NaOH]/ M']/1000

edge_list = repository_dir/'RESOURCES/dendrogram_edgelist.csv'
coord_list = repository_dir/'RESOURCES/dendrogram_coordinates.csv'
G = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G,pos)
c_ops.normalise_network_coordinates(G)
# c_ops.rotate_network(G, -np.pi/4 - np.pi/6)
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)

factor = 800
for c,v in enumerate(condition_dict,0):
    fig, ax = plt.subplots(figsize = (3.9/2.54,2.81/2.54), frameon = False)
    ax.plot(lines[0],lines[1],
           c = '#000000', linewidth = 0.5,
           zorder = 0)
    leaf_colours = []
    leaf_x = []
    leaf_y = []
    for n in G.nodes:
        if n in exp_info.index:
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
    cbar.ax.tick_params(labelsize=7, length = 2, pad = 1)
    cbar.ax.set_position([0.0,0.25,0.95,0.05])
    ax.set_position([0.0,0.3,0.95,0.65])
    cbar.set_label(condition_names[c], fontsize = 8, labelpad = 2)
    ax.set_axis_off()
    # the y-axis is flipped upside down
    # not necessary, but I forgot to remove it
    # before making the final figures, so it
    # gets left in!
    ylm = ax.get_ylim()
    ax.set_ylim(ylm[1],ylm[0])
    plt.savefig(repository_dir/'PLOTS/{}_Figure2B_inset.png'.format(v.split('/')[0]), dpi = 600)
    plt.savefig(repository_dir/'PLOTS/{}_Figure2B_inset.svg'.format(v.split('/')[0]))
    plt.close()
