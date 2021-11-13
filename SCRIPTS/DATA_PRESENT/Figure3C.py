'''
A heatmap illustrating the relative reaction class occurences in modulated
data sets. Figure 3C.
'''
import os
import sys
import pickle
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

from NorthNet import Classes

import helpers.chem_info as info_params
from helpers.network_load_helper import load_reaction_list
from helpers.network_load_helper import convert_to_networkx

# name for output files
figname = 'Figure3C'
# set paths to files
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

reaction_expression = pd.read_csv(
    repository_dir/'RESOURCES/reaction_expression_normalised.csv', index_col = 0
)
reaction_expression = reaction_expression.dropna(axis = 1)

class_names = reaction_expression.columns
experiment_names = reaction_expression.index
expression_array = reaction_expression.to_numpy()[:,:-1]

from mpl_toolkits.axes_grid1 import make_axes_locatable
fig, ax = plt.subplots(figsize = (20/2.54,9/2.54))

ax.tick_params(axis = 'both', which = 'both', length = 0)
im = ax.imshow(expression_array.T, cmap = 'cividis')

ax.set_xticks([])
# ax.set_xticks(np.arange(0,len(experiment_names),1))
# ax.set_xticklabels(experiment_names, fontsize = 3, rotation = 45)

ax.set_yticks(np.arange(0,len(class_names),1))
ax.set_yticklabels(class_names, fontsize = 6)

ax.set_xlabel('Experiment', fontsize = 9)
ax.set_ylabel('Reaction class', fontsize = 9)

divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="5%", pad=0.3)
cbar = plt.colorbar(im, cax=cax, orientation = 'horizontal')
cbar.set_label('Fractional expression',labelpad = 5, fontsize = 9)
cbar.ax.tick_params(labelsize= 6, length = 1)

fig.tight_layout()
plt.savefig(repository_dir/'PLOTS/{}.png'.format(figname), dpi = 600)
plt.savefig(repository_dir/'PLOTS/{}.svg'.format(figname))
plt.close()
