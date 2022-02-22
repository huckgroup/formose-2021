'''
Functions for plotting graph representations.
'''
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from networkx.drawing.nx_agraph import graphviz_layout

from helpers import chem_info as info_params
from NorthNet.network_visualisation import coordinates as c_ops

def plot_nodes(G, ax, color = '#000000', size = 2, alpha = 1, zorder = 1):
    '''
    Draw nodes as a scatter for a network
    
    Parameters
    ----------
    G: networkx Graph or DiGraph
        Graph containing nodes
    ax: matplotlib axis object
        axis in which nodes will be drawn.
    color: str
        if 'nodewise' colours nodes by 'color' attribute. Otherwise
        this parameter sets the colour for all nodes.
    size: float or 'nodewise'
        Node sizes. if 'nodewise' sizes are set according to
        the nodes' 'size' attribute. Otherwise uses single size.
    alpha: float
        Alpha values for points
    zorder: int
        Z-order of the points

    Returns
    -------
    None
    '''

    if color == 'nodewise':
        colors = [G.nodes[n]['color'] for n in G.nodes]
    else:
        colors = color

    if size == 'nodewise':
        sizes = [G.nodes[n]['size'] for n in G.nodes]
    else:
        sizes = size


    coords = c_ops.get_network_scatter(G)

    ax.scatter(coords[0], coords[1],
               s = sizes, c = colors, alpha = alpha,
               zorder = zorder)

def draw_arrow_connectors(G,ax, color = '#000000',linew = 'edgewise', alpha = 1,
                          zorder = 0, shrink_a = 0, shrink_b = 0):
    '''
    Draw arrow connectors for a directed graph network

    Parameters
    ----------
    G: networkx DiGraph with position info in nodes
    ax: matplotlib axis in which to place arrows.
    color: str or 'edgewise'
        if 'edgewise' colour edges according to their 'color' attribute
    linew: float or 'edgewise'
        if 'edgewise', linewidths are set according to the edge 'weight'
        attribute
    alpha: float
    zorder: int
    shrink_a: float
    shrink_b:float

    Returns
    -------
    None
    '''

    if color == 'edgewise':
        color_return = lambda x: G.edges[x]['color']
    else:
        color_return = lambda _: color

    if linew == 'edgewise':
        pass
    else:
        for e in G.edges:
            G.edges[e]['weight'] = linew

    for e in G.edges:
        begin = G.nodes[e[0]]['pos']
        end = G.nodes[e[1]]['pos']

        arrow = FancyArrowPatch(begin,end,arrowstyle='-|>', path = None,
                                connectionstyle='arc',#'Angle3'
                                zorder = zorder,
                                facecolor = color_return(e),
                                edgecolor = color_return(e),
                                linewidth = G.edges[e]['weight'],
                                mutation_scale = 5,
                                shrinkA = shrink_a,
                                shrinkB = shrink_b,
                                alpha = alpha)
        ax.add_patch(arrow)

def plot_network(G, filename, prog = 'neato'):
    '''
    Plot the network G.

    Parameters
    ----------
    G: networkx DiGraph
        Network to plot.
    filename: str
        Name for the file containing the graph.
    prog: str
        Assigment of the layout program to use.

    Returns
    -------
    None
    '''

    pos = graphviz_layout(G, prog = prog)

    c_ops.set_network_coords(G, pos)
    c_ops.normalise_network_coordinates(G)

    for n in G.nodes:
        if '>>' in n:
            G.nodes[n]['color'] = '#000000'
        elif n in info_params.colour_assignments:
            G.nodes[n]['color'] = info_params.colour_assignments[n]
        else:
            G.nodes[n]['color'] = '#EAD5A7'

    fig, ax = plt.subplots()

    plot_nodes(G, ax, color = 'nodewise', size = 'nodewise',
                         alpha = 1, zorder = 1)
    draw_arrow_connectors(G,ax, color = '#000000',

                            linew = 1, alpha = 1,
                             zorder = 1,
                             shrink_a = 0, shrink_b = 0)

    ax.set_axis_off()
    plt.savefig(filename, dpi = 600)
    plt.close()

def plot_network_in_axis(G, ax, layout_provided = False, prog = 'neato'):
    '''
    Add a plot of G into ax

    Parameters
    ----------
    G: networkx DiGraph
        Graph.
    ax: matplotlib axis
    layout_provided: bool
        If a layout for the graph has already been calculated.
    prog: str
        Assigment of the layout program to use.

    Returns
    -------
    None
    '''

    if not layout_provided:
        pos = graphviz_layout(G, prog = prog)
        c_ops.set_network_coords(G, pos)

    c_ops.normalise_network_coordinates(G)

    for n in G.nodes:
        if '>>' in n:
            G.nodes[n]['color'] = '#000000'
        elif n in info_params.colour_assignments:
            G.nodes[n]['color'] = info_params.colour_assignments[n]
        else:
            G.nodes[n]['color'] = '#EAD5A7'

    plot_nodes(G, ax, color = 'none', size = 'nodewise',
                         alpha = 1, zorder = 1)
    draw_arrow_connectors(G,ax, color = 'edgewise',
                            linew = 0.5, alpha = 1,
                             zorder = 1,
                             shrink_a = 0, shrink_b = 0)
    ax.set_axis_off()
