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

from NorthNet import info_params
from NorthNet.calculations import calculations
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

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
                   'Temperature/ $^\circ$C']
condition_dict = {}
for s in select_conditions:
    condition_dict[s] = {}
    for e in exp_info.index:
        condition_dict[s][e] = exp_info.loc[e,s]

edge_list = repository_dir/'RESOURCES/dendrogram_edgelist.csv'
coord_list = repository_dir/'RESOURCES/dendrogram_coordinates.csv'
G = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G,pos)
c_ops.normalise_network_coordinates(G)
# c_ops.rotate_network(G, -np.pi/4 - np.pi/6)
lines = c_ops.get_network_lineplot(G)

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

plt.savefig(repository_dir/'PLOTS/dendrogram_panel.png', dpi = 600)
plt.close()
