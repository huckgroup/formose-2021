'''
A dendrogram with the compositions for each data set annotated on its leaves.
Figure 2B, main panel
'''
import sys
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet.network_visualisation import coordinates as c_ops

import helpers.chem_info as info_params
from helpers.loading_helper import get_carbon_inputs
from helpers.network_load_helper import load_from_edge_list
from helpers.network_load_helper import load_coordinates_list

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

exp_info = pd.read_csv(exp_info_dir, index_col = 0)

average_data = pd.read_csv(derived_parameters_dir/'AverageData.csv', index_col = 0)

# remove empty columns
average_data = average_data.dropna(axis = 1)
# remove columns containing only zeros
average_data = average_data.loc[:, (average_data != 0).any(axis=0)]

carbon_inputs = get_carbon_inputs(exp_info, average_data.columns)

# remove reactants
for c in carbon_inputs:
    for x in carbon_inputs[c]:
        if x in average_data.columns:
            average_data.loc[c,x] = 0.0

average_sum = average_data.sum(axis = 1)

compounds = [x.split('/')[0] for x in average_data.columns]
compound_clrs  = [info_params.colour_assignments[x] for x in compounds]

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
lines = c_ops.get_network_lineplot(G)
dots  = c_ops.get_network_scatter(G)

fig, ax = plt.subplots(figsize = (18.7/2.54, 8/2.54), frameon = False)

ax.plot(lines[0],lines[1],
       c = '#000000', linewidth = 2.5,
       zorder = 0,
       solid_capstyle="round")

subplot_width = 0.05
subplot_height = 0.05
for n in G.nodes:
    if n in exp_info.index:
        xy = G.nodes[n]['pos']
        L = xy[0] - subplot_width/2
        B = xy[1] - subplot_height/2
        W = subplot_width
        H = subplot_height
        axin = ax.inset_axes([L,B,W,H], transform=ax.transData)

        axin.pie(average_data.loc[n].to_numpy()/average_sum.loc[n],
                 colors = compound_clrs)
        # annotate experiment entry number
        # axin.annotate(exp_info.loc[n,'Experiment_entry'],
        #             xy = (0,0),
        #             ha = 'center', va = 'center')

        # axin.annotate(n,
        #             xy = (0,0),
        #             ha = 'center', va = 'center')

ax.set_position([0,0,1,1])
ax.set_axis_off()
plt.savefig(repository_dir/'PLOTS/Figure2B_Main.png', dpi = 600)
plt.savefig(repository_dir/'PLOTS/Figure2B_Main.svg')
plt.close()
