import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

import helpers.chem_info as info_params

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

# read in experiment conditions
exp_info = pd.read_csv(exp_info_dir, index_col = 0)

# load in average concentration data
average_data = pd.read_csv(derived_parameters_dir/'AverageData.csv', index_col = 0)
# remove empty columns
average_data = average_data.dropna(axis = 1)
# create a list of compounds
compounds = list(average_data.columns)

# Choose a parameter from Experiment_parameters.csv
x_key = 'residence_time/ s'
y_key = '[CaCl2]/ M'
x_axis = exp_info.loc[:,x_key].to_numpy()
y_axis = exp_info.loc[:,y_key].to_numpy()

# 3D plot choice (2D: use x_axis)
threeD = False
if threeD:
    independents = [x_axis,y_axis]
    subplot_arg = {'projection':'3d'}
else:
    independents = [x_axis]
    subplot_arg = None

fig, ax = plt.subplots(subplot_kw = subplot_arg)
for c in compounds:
    clr = info_params.colour_assignments[c.split('/')[0]]
    if average_data.loc[:,c].sum() > 0.0:
        ax.plot(*independents, 1000*average_data.loc[:,c],
                '-o',
                linewidth = 1,
                c = clr)

ax.set_xlabel(x_key)
ax.set_ylabel('Concentration/ mM')
plt.show()
