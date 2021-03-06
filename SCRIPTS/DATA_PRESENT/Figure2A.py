'''
Example bar charts of compositional data. Figure 2A
'''
import sys
from pathlib import Path

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

import pandas as pd
import matplotlib.pyplot as plt

import helpers.chem_info as info_params

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'

data_set_selections = ['FRN094C', 'FRN089A', 'FRN090B','FRN077C']

df = pd.read_csv(derived_parameters_dir/'AverageData.csv', index_col = 0)
# remove empty columns
df = df.dropna(axis = 1)

# create a data frame with the selected compounds
sel = df.loc[data_set_selections]
# remove columns containing only zeros
sel = sel.loc[:, (sel != 0).any(axis=0)]

# convert to numpy
plot_stack = sel.to_numpy()

comp_ax = [x.split('/')[0] for x in sel.columns]
comp_clrs  = [info_params.colour_assignments[x] for x in comp_ax]
comp_numbers = [info_params.compound_numbering[x] for x in comp_ax]

# define plot parameters
class margins: left = 0.08;right = 0.01;top = 0.15;bottom = 0.025;
left = margins.left
width = 1-margins.left-margins.right
height = (1-margins.top-margins.bottom)/len(data_set_selections)
axis_label_font_size = 8

# create the plot
fig, ax = plt.subplots(nrows = len(data_set_selections),
                        figsize = (17.5/2.54,3.7/2.54))

for x in range(0,len(plot_stack)):
    ax[x].bar(comp_numbers, plot_stack[x],
                color = comp_clrs)
    if x > 0:
        ax[x].set_xticklabels([])
        ax[x].tick_params(axis = 'x', which = 'both',
                         length = 0.0)
    ax[x].set_yscale('log')
    box = ax[x].get_position([])
    bottom = margins.bottom + (x+1)*(1-margins.top-margins.bottom)/len(data_set_selections) - height/2
    ax[x].set_position([left, bottom, width, height])


for a in ax:
    a.tick_params(axis = 'y', which = 'major',
                length = 2,
                labelsize = 7,
                pad = 1,
                )
    a.tick_params(axis = 'both', which = 'minor',
                length = 1,
                labelsize = 7,
                pad = 1,
                )
    a.tick_params(axis = 'x', which = 'major',
                length = 1,
                labelsize = 7,
                pad = 1,
                )

    a.set_xticklabels([])

fig.text((margins.left - margins.right + 1)/2, 0.04,
         'Compound', ha='center',
         fontsize = axis_label_font_size)
fig.text(margins.left/3, 0.5, 'log$_{10}$(concentration/ M)',
            va='center', ha = 'center', rotation='vertical',
            fontsize = axis_label_font_size)

plt.savefig(plot_folder/'Figure2A.png', dpi = 600)
plt.savefig(plot_folder/'Figure2A.svg')
plt.close()

# Output source data
csv_stack = plot_stack.T
output_text = 'compound,'
for x in range(0,len(data_set_selections)):
    output_text += f"concentration (subplot {x+1})/ M,"
output_text += '\n'

for x in range(0,len(csv_stack)):
    output_text += f"{comp_ax[x]},"
    for y in range(0, len(csv_stack[x])):
        output_text += f"{csv_stack[x,y]},"
    output_text += "\n"

filename = repository_dir/"FIGURE_SOURCE_DATA/Figure2A_Source_Data.csv"
with open(filename, 'w') as file:
    file.write(output_text)

