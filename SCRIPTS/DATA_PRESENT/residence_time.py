'''
A selection of data sets showing how the concentration of formaldehyde
affects the reaction composition. Figure 2C.
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

average_data = pd.read_csv(derived_parameters_dir/'AverageData.csv', index_col = 0)
# remove empty columns
average_data = average_data.dropna(axis = 1)

compound_numbering = info_params.compound_numbering

residence_time_experiments = ['FRN092B','FRN089B','FRN091','FRN092A']

x_axis = [exp_info['residence_time/ s'][r] for r in residence_time_experiments]

data_series = np.zeros((len(residence_time_experiments),len(average_data.columns)))
compound_list = average_data.columns
datasets = average_data.index.to_numpy()
for c,r in enumerate(residence_time_experiments):
    i = np.where(datasets == r)[0]
    data_series[c] = average_data.iloc[i]

y_factor = 1000
width = 10
height = 7
fig, ax = plt.subplots(figsize=(width/2.54,height/2.54))
ax.set_position([0.2, 0.2, 0.7, 0.7])
data_plot = data_series.T
for x in range(0,len(data_plot)):
    if np.sum(data_plot[x]) == 0.0:
        continue
    name = compound_list[x].split('/')[0]
    if name == 'C=O' or name == 'O=C(CO)CO':
        continue
    colour = info_params.colour_assignments[name]
    ax.plot(x_axis, data_plot[x]*y_factor, '-o', c = colour,
            markersize = 5, zorder = 100-name.count('C'))

ax.tick_params(which = 'both', axis = 'both', length = 2)
ax.set_xlabel('Residence time/ s')
ax.set_ylabel('Concentration/ mM')
figname = 'Residence_time'
plt.savefig(repository_dir/'PLOTS/{}_series.svg'.format(figname), dpi = 600)
plt.savefig(repository_dir/'PLOTS/{}_series.png'.format(figname), dpi = 600)
