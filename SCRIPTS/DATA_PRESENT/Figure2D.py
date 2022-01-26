'''
Dengrogrma heat map of formaldehyde concentration 
with formaldehyde induced compositional transitions 
overlayed.
'''
import sys
import pandas as pd
from pathlib import Path
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet.network_visualisation import coordinates as c_ops

from helpers.network_load_helper import load_from_edge_list
from helpers.network_load_helper import load_coordinates_list

# set up paths
data_folder = repository_dir/'DATA'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

# Load experiment conditions
exp_info = pd.read_csv(exp_info_dir, index_col = 0)

# Information to use in extracting information
formaldehyde_key = '[C=O]/ M'
condition_names = '[formaldehyde]/ mM'
figname = 'Figure2D'
series_sel = 'Formaldehyde_2_series'
condition_sel = '[C=O]/ M'

factor = 800
x_factor = 1000
y_factor = 1000

# Extract the required information for the plot
experiment_dict = {}
experiment_codes = exp_info.index
for e in experiment_codes:
    experiment_dict[e] = exp_info.loc[e,formaldehyde_key]

# Load information for the series sequences
series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)
data_keys = series_seqs.loc[series_sel]
data_set_selections = list(data_keys.dropna())

series_x_values = [x_factor*exp_info.loc[x,condition_sel]
                                            for x in data_set_selections]

node_path = [(a,b) for a,b in zip(data_set_selections,data_set_selections[1:])]

# Load in reaction network
edge_list = repository_dir/'RESOURCES/dendrogram_edgelist.csv'
coord_list = repository_dir/'RESOURCES/dendrogram_coordinates.csv'
G = load_from_edge_list(edge_list)
pos = load_coordinates_list(coord_list)
c_ops.set_network_coords(G,pos)
c_ops.normalise_network_coordinates(G)
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)

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

# Set colour for the arrows
arrow_colour = '#5B6F8D'

# create the plot
fig, ax = plt.subplots(figsize=(8.965/2.54,5.6/2.54))

# Plot lines for the dendrogram
ax.plot(lines[0],lines[1],linewidth = 1.5, 
        zorder = 0, c = '#000000')

for n in data_set_selections:
    pos = G.nodes[n]['pos']
    ax.scatter(pos[0],pos[1], c = None,
            facecolor = 'none',
            edgecolor = arrow_colour,
            alpha = 1.0)

# Plot the coloured leaf nodes 
scattr = ax.scatter(leaf_x,leaf_y, s = 10, c = leaf_colours,
                        cmap = cm.coolwarm,
                        zorder = 1)

# Colour bar for the nodes
cbar = plt.colorbar(scattr,
                    ax = ax,
                    location="bottom", aspect = 30)

# add in the path formaldehyde induces across the dendrogram
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
cbar.set_label(formaldehyde_key, fontsize = 8, labelpad = 2)

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
