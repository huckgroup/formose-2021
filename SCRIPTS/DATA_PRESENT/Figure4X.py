import sys
import numpy as np
import pandas as pd
import networkx as nx
from rdkit import Chem
from pathlib import Path
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from helpers import chem_info as info_params
from helpers import pathway_helpers as path_hlp
from helpers.loading_helper import get_carbon_inputs
from helpers.network_load_helper import load_network_from_reaction_list
from helpers.network_load_helper import load_reaction_list,convert_to_networkx
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')
network_file = repository_dir/'FORMOSE_REACTION/FullFormoseReaction.txt'

series_sel = 'Formaldehyde_paper_series'
condition_sel_x = '[C=O]/ M'
factor = 1000



min_font = Global_formatting.min_font # should be around 1 mm on paper
col_width = Global_formatting.col_width
two_col_width = Global_formatting.two_col_width
full_page = Global_formatting.full_page
modest_display_item = Global_formatting.modest_display_item
composite_figure = Global_formatting.composite_figure
lines = Global_formatting.lines/2
ticklength = Global_formatting.ticklength
tickpad = Global_formatting.tickpad

tr_col_map = Global_formatting.cluster_colour_map
roman_numerals = ['I','II','III','IV','V','VI','VII','VIII','IX','X']

series_choice  ="Formaldehyde_paper_series"

directory = Path('/Users/williamrobinson/documents/nijmegen/Dynamic_environment_project')
plot_dir = Path(directory/'Paper_plots')
network_dir = Path(directory/'Networks')
temp_folder = Path('/Users/williamrobinson/Documents/Nijmegen/Packages/Scripts/RPSO_workflow/temp')
data_dir = directory/'DataAnalysis'

n_removals = ["C=O","O","[OH-]"]
base_linew = 0.6
clust_linew = 0.5
s_f = 9000
inclusion_thresh = 0.1
enolate_patt = Chem.MolFromSmarts('C=C')

home = os.getcwd()

exp_info = f_io.import_Experiment_information(directory/"Experiment_parameters.csv")
exp_info = {e:exp_info[e] for e in exp_info if 'FRN068' not in e}
exp_list = [*exp_info]
info_dict = f_io.get_Series_Information(directory/"Series_info.csv")
sets = info_dict[series_choice]
clusters = f_io.get_composition_base_vectors(data_dir/"cluster_results_all_points.csv")

roman_numerals = roman_numerals[:len(clusters)]
roman_numerals.reverse()

off_vec = f_io.get_vectors_from_file(data_dir/'AverageData.csv')
amp_vec = f_io.get_vectors_from_file(data_dir/'AmplitudeData.csv')
header = [x+' M' for x in info_params.smiles_to_names]

reactotypes = f_io.get_reaction_clusters(data_dir/"network_clusters.csv")

info_params.colour_assignments['OC[C@H](O)[C@@H](O)[C@H](O)CO'] = info_params.colour_assignments['OC[C@H](O)C(O)[C@H](O)CO']
for cnt, h in enumerate(header):
    if h == 'OC[C@H](O)C(O)[C@H](O)CO M':
        header[cnt] = 'OC[C@H](O)[C@@H](O)[C@H](O)CO M'

remove_inds = [header.index("C=O M"), header.index("O=C(CO)CO M"),
               header.index("O=CC(O)(CO)CO M"), header.index("O=CCO M"),
               header.index('O=C[C@H](O)CO M')]

for a in off_vec:
    for i in remove_inds:
        off_vec[a][i] = 0.0

# get reaction colours
r_colours = {}
with open(r'/Users/williamrobinson/documents/nijmegen/packages/info_files/reaction_colour_assignments.csv', 'r') as f:
    for c,line in enumerate(f):
        if c == 0:
            pass
        else:
            ins = line.strip('\n').split(',')
            r_colours[ins[0]] = ins[1]

DetectedCompoundNumbering = {}
with open(data_dir/'DetectedCompoundNumbering.txt', 'r') as f:
    for line in f:
        ins = line.strip('\n').split(',')
        DetectedCompoundNumbering[ins[0]] = ins[1]
DetectedCompoundNumbering['OC[C@H](O)[C@@H](O)[C@H](O)CO'] = DetectedCompoundNumbering['OC[C@H](O)C(O)[C@H](O)CO']

wghts = []
for a in sets:
    phenotype,resid = fitting_functions.fit_vectors_to_data(off_vec[a], clusters, header)
    wghts.append(phenotype)
wghts = np.array(wghts)

'''Arrange reactotypes'''
reactotype_nets = []
for r in reactotypes:
    n_net = n_gen.network_from_reaction_list(reactotypes[r],network_name = "{}".format(r))
    reactotype_nets.append(conv.convert_to_networkx(n_net))

F1 = nx.DiGraph()
for x in range(0,len(wghts)):
    for y in range(0,len(wghts[x])):
        if wghts[x,y]/np.amax(wghts[x]) > inclusion_thresh:
            n = reactotype_nets[y]
            for node in n_removals:
                if node in n.nodes:
                    n.remove_node(node)
            F1 = nx.compose(F1,n)

'''
Get all schemes
'''
schemes = f_io.load_schemes(network_dir)
F = nx.DiGraph()
nets = {}
northnets = {}
F = nx.compose(F,F1)
for r in sets:

    r_network = n_gen.network_from_reaction_list(schemes[r],
                                                 network_name = "{}".format(r))

    northnets[r] = r_network

    G = conv.convert_to_networkx(r_network)

    for node in n_removals:
        if node in G.nodes:
            G.remove_node(node)

    nets[r] = G

    F = nx.compose(F,G)

pos = graphviz_layout(F, prog="neato",
                      args="-Goverlap=scale,-GK=1, -Gnodesep={}".format(7000))

for node in F.nodes:
    F.nodes[node]['pos'] = pos[node]


for r in reactotype_nets:
    n_v.add_coords_to_network(r,pos)

base_net = n_v.get_network_coordinates(F)
x_center = np.nanmean(base_net[0])
y_center = np.nanmean(base_net[1])
net_width = np.nanmax(base_net[0])-np.nanmin(base_net[0])
net_height = np.nanmax(base_net[1])-np.nanmin(base_net[1])

for p in pos:
    b = (pos[p][0]-x_center)/net_width
    a = (pos[p][1]-y_center)/net_height
    pos[p] = (a,b)

for n in nets:
    for node in nets[n].nodes:
        nets[n].nodes[node]['pos'] = pos[node]
    for edge in nets[n].edges:
        for e in edge:
            if '>>' in e: col = r_colours[e]
        nets[n].edges[edge]['color'] = col

for node in F.nodes:
    F.nodes[node]['pos'] = pos[node]

for r in reactotype_nets:
    n_v.add_coords_to_network(r,pos)

base_net = n_v.get_network_coordinates(F)

chemo_sum = np.sum(wghts, axis = 0)
nets_to_include  = [x for x in range(0,len(chemo_sum)) if chemo_sum[x] > 0.1*np.sum(chemo_sum)]

'''Plotting series in four panels'''
fig,ax = plt.subplots(nrows = 2, ncols = 2, figsize = (7.91/2.54,8.42/2.54))
plt_pos_seq = [(0,0),(0,1),
               (1,0),(1,1)]

for c,v in enumerate(sets):
    ax[plt_pos_seq[c]].plot(base_net[0],base_net[1], c = '#acb5ad',
                            linewidth = base_linew,zorder =0)

    x,y = n_v.get_network_coordinates(nets[v])

    for r in northnets[v].NetworkReactions:
        ax[plt_pos_seq[c]].scatter(pos[r][0],pos[r][1], marker = 'd',
                                   zorder = 2, s = clust_linew*4, c = 'w',
                                   edgecolors = '#000000',
                                   facecolors = 'w')

    plotting_operations.draw_arrow_connectors(nets[v],ax[plt_pos_seq[c]],
                                              color = 'edgewise',
                                              linew = clust_linew*1.5,
                                              zorder = 1)

    for c2,h in enumerate(header):

        if off_vec[v][c2] > 0 and h[:-2] not in n_removals:

            ax[plt_pos_seq[c]].scatter(pos[h[:-2]][0],pos[h[:-2]][1],
                                        s =  off_vec[v][c2]*s_f,
                                        c = info_params.colour_assignments[h[:-2]],
                                        zorder = 2, edgecolors = 'None')

            ax[plt_pos_seq[c]].scatter(pos[h[:-2]][0],pos[h[:-2]][1],
                                        s =  (off_vec[v][c2]+amp_vec[v][c2])*s_f,
                                        edgecolors = info_params.colour_assignments[h[:-2]],
                                        zorder = 2,facecolors='None',linewidth = base_linew)
            #ax[plt_pos_seq[c]].annotate(DetectedCompoundNumbering[h[:-2]], xy = (pos[h[:-2]][0],pos[h[:-2]][1]), fontsize = 6)
    '''
    for en_n in nets[v].nodes:
        mol = Chem.MolFromSmiles(en_n)
        if mol == None:
            continue
        if mol.HasSubstructMatch(enolate_patt):
            ax[plt_pos_seq[c]].scatter(pos[en_n][0],pos[en_n][1], marker = 11,
                                       zorder = 2, s = clust_linew*2, c = '#E1D6B7')
    '''
    networks_to_plot = []
    sort_w = np.array([])
    for w in range(0,len(wghts[c])):
        if wghts[c,w]/np.amax(wghts[c]) > inclusion_thresh:
            networks_to_plot.append((reactotype_nets[w], tr_col_map[w], (wghts[c,w]/np.amax(wghts[c])), roman_numerals[w]))
            sort_w = np.hstack((sort_w,wghts[c,w]/np.amax(wghts[c]) ))

    networks_to_plot = sorted(networks_to_plot, key = lambda x:x[2])

    L = 0.8
    B = 0.75
    W = 0.2
    H = 0.2

    for rc,n in enumerate(networks_to_plot):
        G = n[0]
        if len(G.nodes) < 2:
            pass
        else:
            axin = ax[plt_pos_seq[c]].inset_axes([L,B,W,H])
            axin.set_xticks([])
            axin.set_yticks([])
            [axin.spines[i].set_linewidth(0.1) for i in axin.spines]
            B -= 0.26
            x,y = n_v.get_network_coordinates(G)

            axin.plot(x,y,linewidth = 1,color = n[1],
                                         zorder = 1, alpha = n[2],
                                         solid_capstyle='round')
            axin.plot(base_net[0],base_net[1], c = '#acb5ad',
                                    linewidth = base_linew/1.5,zorder = 0)

            axin.annotate(n[-1], xy = (0, -0.9),
                          annotation_clip=False, fontsize = min_font*2,
                          xycoords='data', horizontalalignment = 'center')



for x in plt_pos_seq:
    ax[x].set_axis_off()
    L = (x[1]/2)
    B = 0.5-x[0]/2
    W = 0.45
    H = 0.5
    ax[x].set_position([L,B,W,H])
    ylm = ax[x].get_ylim()
    ax[x].set_ylim(ylm[0]-0.1, ylm[1])

plt.savefig(plot_dir/'Figure3b.png'.format(series_choice), dpi = 800,
            transparent = False)
plt.close()

# make a scale bar for the compound nodes
fig,ax = plt.subplots(figsize = (8.39/2.54,0.93/2.54))
two_mM = s_f*2/1000
five_mM = s_f*5/1000
ten_mM = s_f*10/1000
twenty_mM = s_f*20/1000

ax.scatter(0,0,    s = two_mM,    c = '#acb5ad', edgecolors = 'None')
ax.scatter(0.25,0, s = five_mM,   c = '#acb5ad', edgecolors = 'None')
ax.scatter(0.5,0,  s = ten_mM,    c = '#acb5ad', edgecolors = 'None')
ax.scatter(0.75,0, s = twenty_mM, c = '#acb5ad', edgecolors = 'None')

ax.scatter(0,0, s = two_mM + five_mM,
                edgecolors = '#acb5ad',
                zorder = 2, facecolors='None', linewidth = base_linew)

ax.annotate('2 mM', xy = (0.05,0), fontsize = min_font*2,
                 verticalalignment = 'center', horizontalalignment = 'left')
ax.annotate('5 mM', xy = (0.3,0), fontsize = min_font*2,
                 verticalalignment = 'center', horizontalalignment = 'left')
ax.annotate('10 mM', xy = (0.55,0), fontsize = min_font*2,
                 verticalalignment = 'center', horizontalalignment = 'left')
ax.annotate('20 mM', xy = (0.8,0), fontsize = min_font*2,
                 verticalalignment = 'center', horizontalalignment = 'left')

ax.set_xlim(-0.1,1)
ax.set_position([0,0.1,1,0.9])
ax.set_axis_off()
plt.savefig(plot_dir/'Figure_3b_key.png', dpi = 600, transparent = True)
