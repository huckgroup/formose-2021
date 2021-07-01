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

figname = 'ResidenceTimeSeries'
series_sel = 'Residence_time_2_series'
condition_sel = 'residence_time/ s'
x_name = 'residence time/ s'
x_factor = 1
y_factor = 1000

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

series_stack = sel.to_numpy()

series_progression = series_stack.T

series_x_values = [x_factor*exp_info.loc[x,condition_sel]
                                            for x in data_set_selections]

compounds = [x.split('/')[0] for x in sel.columns]
compound_clrs  = [info_params.colour_assignments[x] for x in compounds]

fig, ax = plt.subplots(figsize=(8.965/2.54,6.55/2.54))
ax.set_position([0.2, 0.2, 0.65, 0.7])

trace_labels = []
for x in range(0,len(series_progression)):
    if np.sum(series_progression[x]) == 0.0:
        continue
    species_name = compounds[x].split('/')[0]
    if species_name in info_params.colour_assignments:
        clr = info_params.colour_assignments[species_name]
    else:
        clr = 'k'
    ax.plot(series_x_values,series_progression[x]*y_factor, '-o',
            c = clr, markersize = 5, zorder = 100-species_name.count('C'))

    max_idx = np.where(series_progression[x] ==series_progression[x].max())[0]
    trace_labels.append((compound_numbering[species_name], clr))


    x_pos = series_x_values[max_idx[0]]
    y_pos = y_factor*series_progression[x,max_idx]
    # ax.annotate(compound_numbering[species_name],
    #             xy = (x_pos, y_pos),
    #             zorder = 1000)

ax.set_xlabel(x_name)
ax.set_ylabel('Concentration/ mM')
plt.savefig(repository_dir/'PLOTS/{}_series.png'.format(figname), dpi = 600)

label_x_positions = np.linspace(0,series_x_values[-1], num = len(trace_labels))
label_x_positions*=1.5
label_y_positions = np.full(len(trace_labels),y_factor*series_progression.max())

label_x_positions[13:] = label_x_positions[:len(label_x_positions)-13]
label_y_positions[13:] = y_factor*series_progression.max()*0.9
for c,t in enumerate(trace_labels):
    x_pos = label_x_positions[c]
    y_pos = label_y_positions[c]
    ax.annotate(t[0],
                xy = (x_pos, label_y_positions[c]),
                fontsize = 6,
                fontweight = 'bold',
                color = t[1],
                zorder = 1000,
                ha = 'center',
                va = 'center')
plt.savefig(repository_dir/'PLOTS/{}_series_annotated.png'.format(figname), dpi = 600)

plt.close()
