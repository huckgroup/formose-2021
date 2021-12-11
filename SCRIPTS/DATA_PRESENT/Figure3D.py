'''
A heatmap illustrating the relative reaction class occurences in modulated
data sets. Figure 3D.
'''
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import patches

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

# name for output files
figname = 'Figure3D'
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

series_seqs = pd.read_csv(repository_dir/'EXPERIMENT_INFO/Series_info.csv', index_col = 0)

#series_seqs = series_seqs.dropna(axis = 1)

formaldehyde_series = 'Formaldehyde_2_series'
formaldehyde_series_labels = series_seqs.loc[formaldehyde_series]
formaldehyde_series_labels = formaldehyde_series_labels.dropna()
CaOH2_series = 'Ca_OH_grid_series'
CaOH2_series_labels = series_seqs.loc[CaOH2_series]
CaOH2_series_labels = CaOH2_series_labels.dropna()
series_of_interest = pd.concat(
    [formaldehyde_series_labels, CaOH2_series_labels],
    axis = 0
    )
series_of_interest = series_of_interest.to_list()

'''
print(series_of_interest)
['FRN090B', 'FRN093A', 'FRN093B', 'FRN093C',
'FRN089A', 'FRN089B', 'FRN089C', 'FRN089D', 'FRN103',
'FRN104A', 'FRN104B', 'FRN105A', 'FRN105B', 'FRN088A',
'FRN088B', 'FRN089B']
'''

# break the data sets into branches manually
branches = [
    # branch I
    ['FRN090B', 'FRN093A', 'FRN093B', 'FRN093C'],
    # branch II
    ['FRN089A', 'FRN089B', 'FRN089C', 'FRN089D',
        'FRN088A', 'FRN088B'],
    # branch III
    ['FRN105A', 'FRN105B', 'FRN103'],
    # branch IV
    ['FRN104A', 'FRN104B'] 
]

branch_names = ['Branch I','Branch II',
            'Branch III','Branch IV']

# narrow down the data to the sets used in
# the branches described above
reaction_expression = reaction_expression.loc[series_of_interest]

# remove rows in which no reactions are counted
reaction_expression = reaction_expression.loc[(reaction_expression!=0).any(axis=1)]
#reaction_expression = reaction_expression.drop('Cannizzaro reaction', 1)
reaction_expression = reaction_expression.loc[:, (reaction_expression != 0).any(axis=0)]
class_names = reaction_expression.columns

fig, ax = plt.subplots(
                    figsize = (15/2.54,9/2.54),
                    )

min_val = reaction_expression.to_numpy().min()
max_val = reaction_expression.to_numpy().max()

plot_array = np.array([])
spacer = np.zeros((2,len(class_names)))
for c,b in enumerate(branches):
    selected_region = reaction_expression.loc[b,:]
    expression_array = selected_region.to_numpy()

    if len(plot_array) == 0:
        plot_array = expression_array
    else:
        plot_array = np.vstack((plot_array, expression_array))

    label_position = plot_array.shape[0] - expression_array.shape[0]/2
    ax.annotate(branch_names[c], xy = (label_position,-0.52),
        horizontalalignment = 'center', fontsize = 6)

    if c == len(branches)-1:
        break

    # add in dividers between each branch
    plot_array = np.vstack((plot_array, spacer))

    rect = patches.Rectangle([plot_array.shape[0]-2.5,-0.5], 2,
        len(class_names),
        color = 'w', angle=0.0, zorder = 1)

    ax.add_patch(rect)

im = ax.imshow(
                plot_array.T,
                cmap = 'cividis',
                vmin = min_val,
                vmax = max_val,
                zorder = 0
                )

ax.tick_params(axis = 'both', which = 'both', length = 0)
ax.set_xticks([])
ax.set_yticks(np.arange(0,len(class_names),1))
ax.set_yticklabels(class_names, fontsize = 6)

ax.set_xlabel('Experiment', fontsize = 9)
ax.set_ylabel('Reaction class', fontsize = 9)

cbar = fig.colorbar(im, orientation = 'vertical')
cbar.set_label('Fractional expression',labelpad = 5, fontsize = 9)
cbar.ax.tick_params(labelsize= 6, length = 1)

plt.savefig(repository_dir/'PLOTS/{}.png'.format(figname), dpi = 600)
plt.savefig(repository_dir/'PLOTS/{}.svg'.format(figname))
plt.close()
