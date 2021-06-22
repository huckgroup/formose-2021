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
c_ops.rotate_network(G, -np.pi/4 - np.pi/6)
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)

fig, ax = plt.subplots(figsize = (17.64/2.54, 13.44/2.54), frameon = False)

ax.plot(lines[0],lines[1],
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
        axin = ax.inset_axes([L,B,W,H], transform=ax.transData)

        axin.pie(modified_averages[n]/np.amax(modified_averages[n]),
                 colors = [info_params.colour_assignments[x.split('/')[0]]
                 for x in prime_header])
ylm = ax.get_ylim()
ax.set_ylim(ylm[1],ylm[0])
ax.set_position([0,0,1,1])
ax.set_axis_off()
plt.savefig(repository_dir/'PLOTS/dendrogram_panel.png', dpi = 600)
plt.close()
