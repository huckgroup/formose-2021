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

series_sel = 'Formaldehyde_paper_series'
condition_sel = '[C=O]/ M'
node_size_factor = 100

exp_list = series_dict[series_sel]

idx = [[*exp_info].index(x) for x in series_dict[series_sel]]
series_stack = np.zeros((len(exp_list),len(header)))
for x in range(0,len(idx)):
    series_stack[x] = data[idx[x]]

series_x_values = [exp_info[x].parameters[condition_sel]
                        for x in exp_list]
print(series_x_values)
nets = []
for file in exp_list:
    r_list = load_reaction_list(network_folder/'{}_reaction_list.txt'.format(file))
    n_net = load_network_from_reaction_list(r_list)
    nets.append(convert_to_networkx(n_net))

for c,n in enumerate(nets):
    for node in n.nodes:
        n.nodes[node]['formaldehyde'] = [series_x_values[c]]

F = nx.DiGraph()
for n in nets:
    [n.remove_node(x) for x in ['C=O','O','[OH-]'] if x in n.nodes]
    # F = nx.compose(F,n)

    for node in n.nodes:
        if node in F.nodes:
            F.nodes[node]['formaldehyde'].extend(n.nodes[node]['formaldehyde'])
        else:
            F.add_node(node, formaldehyde = n.nodes[node]['formaldehyde'])

    F.add_edges_from(n.edges)

pos = graphviz_layout(F, prog = 'neato', args='-GK=0.7')
c_ops.set_network_coords(F, pos)
c_ops.normalise_network_coordinates(F)
pos = {n:F.nodes[n]['pos'] for n in F.nodes}

for node in F.nodes:
    if '>>' in node:
        F.nodes[node]['color'] = '#000000'
        F.nodes[node]['size'] = node_size_factor/2
    elif node+'/ M' in prime_header:
        F.nodes[node]['color'] = info_params.colour_assignments[node]
        c_idx = prime_header.index(node+'/ M')
        F.nodes[node]['size'] = node_size_factor
    else:
        F.nodes[node]['color'] = info_params.colour_assignments[node]
        F.nodes[node]['size'] = node_size_factor


for e in F.edges:
    F.edges[e[0],e[1]]['weight'] = F.nodes[e[0]]['formaldehyde'][0]*100

fig, ax = plt.subplots(figsize = (16.85/2.54,16.85/2.54))
for c,n in enumerate(nets):
    x_dots = []
    y_dots = []
    clrs =[]
    sizes = []
    for node in F.nodes:
        if '>>' in node:
            continue
        else:
            x_dots.append(F.nodes[node]['pos'][0])
            y_dots.append(F.nodes[node]['pos'][1])
            clrs.append(F.nodes[node]['color'])
            sizes.append(F.nodes[node]['size'])

    northplot.draw_arrow_connectors(F,ax, color = '#000000',
                                    linew = 'edgewise',
                                    alpha = 0.8,
                                    zorder = 0,
                                    shrink_a = 0,
                                    shrink_b = 0)

    ax.scatter(x_dots, y_dots, c = clrs, s = sizes)
    ax.set_axis_off()
    # for n in F.nodes:
    #     ax.annotate(n, xy = F.nodes[n]['pos'], ha = 'center', fontsize = 6)

output_folder = Path('/Users/williamrobinson/Documents/Nijmegen/Packages/RPSO/FormosePaperAnalysis/Analysis2/plots_for_paper')
plt.savefig(output_folder/'{}_full_network.png'.format(series_sel), dpi = 600)
plt.close()

import pickle
with open('/Users/williamrobinson/documents/nijmegen/packages/rpso/formosepaperanalysis/analysis2/information_sources/FormoseReactionNetwork.pickle', 'rb') as f:
    FRN = pickle.load(f)


reaction_classes = [FRN.NetworkReactions[r].ReactionTemplate.Name for r in FRN.NetworkReactions]
reaction_classes.sort()
reaction_class_container = {r:[] for r in reaction_classes}
for n in F.nodes:
    if '>>' in n:
        class_key = FRN.NetworkReactions[n].ReactionTemplate.Name
        reaction_class_container[class_key].append(n)

with open(output_folder/'full_reaction_list.txt', 'w') as f:
    for r in reaction_class_container:
        f.write(r+'\n')
        for x in reaction_class_container[r]:
            f.write(x + '\n')

        f.write('\n')

for n in F.nodes:
    F.nodes[n]['size'] = ''
    F.nodes[n]['pos'] = ''
    F.nodes[n]['formaldehyde'] = ''
    print(F.nodes[n])
nx.write_gexf(F, output_folder/"formose_series_reactions.gexf")
