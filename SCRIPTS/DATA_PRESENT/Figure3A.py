# generate demonstrative amplitude diagrams for the formose reactions
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes

import helpers.chem_info as info_params

fontprops = fm.FontProperties(size=6)

experiment_code = 'FRN088C'
compound_number_sel = [3,9,14,20,18,12]
compound_numbering = {}
with open(repository_dir/'COMPOUND_INFO/compound_numbering.txt', 'r') as f:
    for line in f:
        ins = line.strip('\n').split(',')
        compound_numbering[int(ins[1])] = ins[0]

loaded_reports = []
path = 'DATA/DATA_REPORTS/{}'.format(experiment_code)
for file in os.listdir(path):
    if file.endswith('csv'):
        file_path = '{}/{}'.format(path,file)
        loaded_reports.append(Classes.DataReport(file = file_path))

data = Classes.DataSet(data_reports = loaded_reports)

fig,ax = plt.subplots(nrows =  len(compound_number_sel), ncols = 1, sharex = True)

class plot_parameters:
    plot_width = 7.93; plot_height = 9.95
    left = 0.12; right = 0.95
    top = 1; bottom = 0.1
    ax_height = 0.05

ax_width = plot_parameters.right-plot_parameters.left
fig.set_figheight(plot_parameters.plot_height/2.54)
fig.set_figwidth(plot_parameters.plot_width/2.54)
ax_height = (plot_parameters.top-plot_parameters.bottom)/len(ax)

for c,n in enumerate(compound_number_sel):
    name = compound_numbering[n]

    x, y = data.get_entry(name + '/ M')

    ax[c].plot(x,1000*y, '-o', markersize = 2.5, linewidth  = 1,
            c = info_params.colour_assignments[name])

for c2 in range(0,len(ax)):
    ax[c2].yaxis.set_major_locator(mticker.MaxNLocator(nbins=3, prune='both',
                                    min_n_ticks  = 3))
    ax[c2].tick_params(axis='both', which='major', labelsize = 4.5,
                       length = 2, pad = 2)

    ax_ypos = plot_parameters.top - plot_parameters.bottom - c2*ax_height

    ax_h = ax_height*0.8
    B = ax_ypos - 0.5 * ax_h

    ax[c2].set_position([plot_parameters.left,B,ax_width,ax_h])

    x_lims = ax[c2].get_xlim()
    y_lims = ax[c2].get_ylim()

fig.text(0.02, 0.55, "Concentration/ mM", va='center', rotation='vertical',
         fontsize = 9)
fig.text(0.55, 0.02, "time/ s", ha='center', fontsize = 9)

plt.savefig(repository_dir/'PLOTS/Figure_3A.png'.format(name), dpi = 600)
plt.close()
