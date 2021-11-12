'''
Insets for insertion into Figure 3
'''
import os
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from NorthNet import Classes

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

import helpers.chem_info as info_params

def processAxis(ax):
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(spline_width)

    ax.set_position([0.1, 0.3, 0.8, 0.7])
    ax.tick_params(axis='both', which='major',
            length = tick_length,
            labelsize=axis_label_size
    )
    #ax.set_xlabel('time/ ks', fontsize = 10)
    #ax.set_ylabel('Concentration/ mM', fontsize = 10)

plot_dir = repository_dir/'PLOTS'
data_dir = repository_dir/'DATA'/'DATA_REPORTS'

plot_width = 3.7/2.54
plot_height = 2.5/2.54
spline_width = 0.5
linewidth = 0.5
markersize = 2 
axis_label_size = 6
tick_length = 0.5
factor = 1000
time_factor = 0.0001

exp_code = 'FRN088C'
compound_number_sel = [3,9,14,20,18,12]

compound_numbering = {}
with open(repository_dir/'COMPOUND_INFO/compound_numbering.txt', 'r') as f:
    for line in f:
        ins = line.strip('\n').split(',')
        compound_numbering[ins[0]] = int(ins[1])

data = Classes.DataReport(
    file = data_dir/exp_code/f'{exp_code}_GCMS_concentration_report.csv'
)


for d in data.data:
    smiles = d.split('/')[0]
    if smiles in compound_numbering:
        compound_no = compound_numbering[smiles]
        if compound_numbering[smiles] in compound_number_sel:
            fig, ax = plt.subplots(figsize = (plot_width, plot_height))
            clr = info_params.colour_assignments[smiles]
            ax.plot(time_factor*data.series_values, factor*data.data[d],
                        '-o', markersize = markersize,
                        c = clr, linewidth = linewidth)
            processAxis(ax)
            plt.savefig(plot_dir/f'Figure3_inset_{compound_no}.svg')
