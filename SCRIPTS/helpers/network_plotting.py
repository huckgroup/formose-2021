'''
Functions for plotting graph representations.
'''
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt

from helpers import chem_info as info_params
from NorthNet.plotting import network_plotting as northplot
from NorthNet.network_visualisation import coordinates as c_ops

def plot_network(G, filename, prog = 'neato'):
    pos = graphviz_layout(G, prog = prog)

    c_ops.set_network_coords(G, pos)
    c_ops.normalise_network_coordinates(G)
    net_coords = c_ops.get_network_scatter(G)

    for n in G.nodes:
        if '>>' in n:
            G.nodes[n]['color'] = '#000000'
        elif n in info_params.colour_assignments:
            G.nodes[n]['color'] = info_params.colour_assignments[n]
        else:
            G.nodes[n]['color'] = '#EAD5A7'

    fig, ax = plt.subplots()
    northplot.plot_nodes(G, ax, color = 'nodewise', size = 'nodewise',
                         alpha = 1, zorder = 1)
    northplot.draw_arrow_connectors(G,ax, color = '#000000',
                            linew = 1, alpha = 1,
                             zorder = 1,
                             shrink_a = 0, shrink_b = 0)
    ax.set_axis_off()
    plt.savefig(filename, dpi = 600)
    plt.close()

def plot_network_in_axis(G, ax, layout_provided = False, prog = 'neato'):

    if not layout_provided:
        pos = graphviz_layout(G, prog = prog)
        c_ops.set_network_coords(G, pos)

    c_ops.normalise_network_coordinates(G)
    net_coords = c_ops.get_network_scatter(G)

    for n in G.nodes:
        if '>>' in n:
            G.nodes[n]['color'] = '#000000'
        elif n in info_params.colour_assignments:
            G.nodes[n]['color'] = info_params.colour_assignments[n]
        else:
            G.nodes[n]['color'] = '#EAD5A7'

    northplot.plot_nodes(G, ax, color = 'none', size = 'nodewise',
                         alpha = 1, zorder = 1)
    northplot.draw_arrow_connectors(G,ax, color = 'edgewise',
                            linew = 0.5, alpha = 1,
                             zorder = 1,
                             shrink_a = 0, shrink_b = 0)
    ax.set_axis_off()
