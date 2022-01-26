'''
Illustration of how the concentration of sodium hydroxide and calcium chloride
induced compositional transitions accross the dendrogram. Figure 2F.
'''
import sys
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib.cm as cm

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet.network_visualisation import coordinates as c_ops
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
factor = 800
x_factor = 1000
y_factor = 1000

experiment_codes = exp_info.index
# Extract the required information for the plot
experiment_dict = {}
experiment_dict['Ca_OH_ratio'] = {}
for e in experiment_codes:
    experiment_dict[e] = exp_info.loc[e,'[CaCl2]/ M']/exp_info.loc[e,'[NaOH]/ M']/1000

series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)
data_keys = series_seqs.loc[series_sel]
data_set_selections = list(data_keys.dropna())

node_path = [(a,b) for a,b in zip(data_set_selections,data_set_selections[1:])]

# create leaf position scatter and colours
leaf_colours = []
leaf_x = []
leaf_y = []
for n in G.nodes:
    if n in experiment_codes:
        xy = G.nodes[n]['pos']
        leaf_x.append(xy[0])
        leaf_y.append(xy[1])
        leaf_colours.append(experiment_dict[n]*factor)

fig, ax = plt.subplots(figsize=(8.965/2.54,5.6/2.54))

scattr = ax.scatter(leaf_x,leaf_y, s = 10, c = leaf_colours,
                        cmap = cm.coolwarm,
                    zorder = 1)

ax.plot(lines[0],lines[1],linewidth = 1.5, 
    zorder = 0, c = '#000000')

# Colour bar for the nodes
cbar = plt.colorbar(scattr,
                    ax = ax,
                    location="bottom", aspect = 30)

for n in data_set_selections:
    pos = G.nodes[n]['pos']
    ax.scatter(pos[0],pos[1], c = None,
            facecolor = 'none',
            edgecolor = 'b',
            alpha = 0.5)

# add in the path formaldehyde induces across the dendrogram
arrow_colour = '#5B6F8D'
for node_pair in node_path:
    arrow = FancyArrowPatch(G.nodes[node_pair[0]]['pos'],
                            G.nodes[node_pair[1]]['pos'],
                            arrowstyle='-|>',
                            path = None,
                            connectionstyle='Angle3',
                            zorder = 1,
                            facecolor = arrow_colour,
                            edgecolor = arrow_colour ,
                            linewidth = 1,
                            mutation_scale = 10,
                            shrinkA = 3,
                            shrinkB = 3,
                            alpha = 1)
    ax.add_patch(arrow)

# plot settings 
ax.set_axis_off()
ylm = ax.get_ylim()

# Colourbar setting
cbar.ax.tick_params(labelsize=7, length = 2, pad = 1)
cbar.set_label('[CaCl$_2$]/[NaOH]', fontsize = 8, labelpad = 2)

# Colorbar to plot ratio
cbar_width = 0.9
cbar_height = 0.1
plot_width = 1
plot_height = 0.85
cbar.ax.set_position([0.05,0.075,cbar_width,cbar_height])
ax.set_position([0.0,0.15,plot_width,plot_height])

# Save output
plt.savefig(repository_dir/'PLOTS/{}.png'.format(figname), dpi = 600)
plt.savefig(repository_dir/'PLOTS/{}.svg'.format(figname), dpi = 600)
plt.close()
