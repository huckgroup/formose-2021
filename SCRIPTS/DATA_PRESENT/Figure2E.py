'''
Subset of the data illustrating how the reaction composition changes with
varying sodium hydroxide and calcium chloride concentrations. Figure 2E.
'''
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

import helpers.chem_info as info_params
from helpers.loading_helper import get_carbon_inputs

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

exp_info = pd.read_csv(exp_info_dir, index_col = 0)

series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)

compound_numbering = info_params.compound_numbering

figname = 'Figure2E'
series_sel = 'Ca_OH_grid_series'
condition_sel_x = '[NaOH]/ M'
condition_sel_y = '[CaCl2]/ M'
factor = 1000

data_keys = series_seqs.loc[series_sel]
data_set_selections = list(data_keys.dropna())

average_data = pd.read_csv(derived_parameters_dir/'AverageData.csv', index_col = 0)
# remove empty columns
average_data = average_data.dropna(axis = 1)

carbon_inputs = get_carbon_inputs(exp_info, average_data.columns)

# remove reactants
for c in carbon_inputs:
    for x in carbon_inputs[c]:
        if x in average_data.columns:
            average_data.loc[c,x] = 0.0

exclusions = ['O=CCO/ M','O=C[C@H]C(O)CO/ M']
for x in exclusions:
    average_data.loc[:,x] = 0.0

sel = average_data.loc[data_set_selections]
# remove columns containing only zeros
sel = sel.loc[:, (sel != 0).any(axis=0)]

series_x_values = [factor*exp_info.loc[x,condition_sel_x]
                        for x in data_set_selections]

series_y_values = [factor*exp_info.loc[x,condition_sel_y]
                        for x in data_set_selections]

compounds = [x.split('/')[0] for x in sel.columns]
compound_clrs  = [info_params.colour_assignments[x] for x in compounds]

series_stack = sel.to_numpy()



subplot_width = 8
subplot_height = 8
fig, ax = plt.subplots(nrows = 1, figsize = (8.965/2.54,6.27/2.54))
ax.scatter(series_x_values, series_y_values, alpha = 0.0)
for x in range(0,len(series_stack)):
    L = series_x_values[x] - subplot_width/2
    B = series_y_values[x] - subplot_height/2
    W = subplot_width
    H = subplot_height
    axin = ax.inset_axes([L,B,W,H], transform=ax.transData)

    wedges, texts = axin.pie(series_stack[x]/series_stack[x].max(),
                        colors = compound_clrs,
                        radius = 1.9,
                        # labels = [compound_numbering[x.split('/')[0]] for x in
                        #         sel.columns],
                        # textprops={'fontsize': 3}
                        )

    for w in wedges:
        w.set_linewidth(0.1)
        w.set_edgecolor('#000000')


ylim = ax.get_ylim()
xlim = ax.get_xlim()
ax.set_position([0.2, 0.2, 0.7, 0.7])
ax.tick_params(which = 'both', axis = 'both', length = 2)
ax.set_ylim(ylim[0]-subplot_width/2, ylim[1]+subplot_width/2)
ax.set_xlim(xlim[0]-subplot_width/2, xlim[1]+subplot_width/2)
ax.set_xlabel('[NaOH]/ mM')
ax.set_ylabel('[CaCl$_2$]/ mM')

ratios = [a/b for a,b in zip(series_y_values, series_x_values)]
ratios = [50, 5, 2, 0.5]
points = [0.2, 2.5, max(series_x_values)]
line_func = lambda a,x: a*x
for r in ratios:
    ax.plot(points, [line_func(r,x) for x in points], '--',
            c = '#000000', alpha = 0.5, linewidth = 0.5)

plt.savefig(repository_dir/'PLOTS/{}.png'.format(figname), dpi = 600)
plt.savefig(repository_dir/'PLOTS/{}.svg'.format(figname))
plt.close()
