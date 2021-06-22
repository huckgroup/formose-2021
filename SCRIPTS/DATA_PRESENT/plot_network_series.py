import numpy as np
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

from NorthNet import info_params
from NorthNet.file_exports.plotting import network_plotting as northplot
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

import __init__
from helpers import load_series
from helpers.loading_helper import header
from helpers.loading_helper import exp_info
from helpers.loading_helper import prime_header
from helpers.loading_helper import header
from helpers.loading_helper import data, experiment_amplitudes, experiment_averages
from helpers.network_load_helper import load_reaction_list
from helpers.network_load_helper import load_network_from_reaction_list, convert_to_networkx

base_directory = Path('/Users/williamrobinson/Documents/Nijmegen')
# base_directory = Path(r'C:\Users\willi\Documents')

fname = base_directory/'safestore_DEP/Series_info.csv'
network_folder = base_directory/'Packages/RPSO/FormosePaperAnalysis/Analysis2/information_sources/network_lists'
series_file = base_directory/'safestore_DEP/Series_info.csv'
series_dict = load_series.load_series_sequences(series_file)

series_sel = 'Ca_OH_grid_series'
condition_sel = '[C=O]/ M'
node_size_factor = 2.5e4

exp_list = series_dict[series_sel]

idx = [[*exp_info].index(x) for x in series_dict[series_sel]]
series_stack = np.zeros((len(exp_list),len(header)))
for x in range(0,len(idx)):
    series_stack[x] = data[idx[x]]

series_x_values = [exp_info[x].parameters[condition_sel]
                        for x in exp_info]

nets = []
for file in exp_list:
    r_list = load_reaction_list(network_folder/'{}_reaction_list.txt'.format(file))
    n_net = load_network_from_reaction_list(r_list)
    nets.append(convert_to_networkx(n_net))

F = nx.DiGraph()
for n in nets:
    [n.remove_node(x) for x in ['C=O','O','[OH-]'] if x in n.nodes]
    F = nx.compose(F,n)

pos = graphviz_layout(F, prog = 'neato', args='-GK=0.7')
c_ops.set_network_coords(F, pos)
c_ops.normalise_network_coordinates(F)
pos = {n:F.nodes[n]['pos'] for n in F.nodes}

for c,n in enumerate(nets):
    c_ops.set_network_coords(n, pos)
    for node in n.nodes:
        if '>>' in node:
            n.nodes[node]['color'] = '#000000'
            n.nodes[node]['size'] = node_size_factor*1e-6
        elif node+'/ M' in prime_header:
            n.nodes[node]['color'] = info_params.colour_assignments[node]
            c_idx = prime_header.index(node+'/ M')
            n.nodes[node]['size'] = experiment_averages[exp_list[c]][c_idx]*node_size_factor
        else:
            n.nodes[node]['color'] = info_params.colour_assignments[node]
            n.nodes[node]['size'] = node_size_factor*1e-6

fig, ax = plt.subplots(ncols = len(exp_list),
                       sharex = True,
                       sharey = True,
                       figsize = (16.85/2.54,16.85/2.54))
for c,n in enumerate(nets):
    x_dots = []
    y_dots = []
    clrs =[]
    sizes = []
    for node in n.nodes:
        x_dots.append(n.nodes[node]['pos'][0])
        y_dots.append(n.nodes[node]['pos'][1])
        clrs.append(n.nodes[node]['color'])
        sizes.append(n.nodes[node]['size'])

    northplot.draw_arrow_connectors(n,ax[c], color = '#000000',
                                    linew = 1.5,
                                    alpha = 0.8,
                                    zorder = 0,
                                    shrink_a = 0,
                                    shrink_b = 0)

    ax[c].scatter(x_dots, y_dots, c = clrs, s = sizes)
    ax[c].set_axis_off()
    ax[c].set_position([c*1/len(nets), 0.05, 1/len(nets), 0.9])
    # ax[c].set_xlim(0,1)
    # ax[c].set_ylim(0,1)

output_folder = Path('/Users/williamrobinson/Documents/Nijmegen/Packages/RPSO/FormosePaperAnalysis/Analysis2/plots_for_paper')
plt.savefig(output_folder/'{}_network_series.png'.format(series_sel), dpi = 600)
plt.close()
