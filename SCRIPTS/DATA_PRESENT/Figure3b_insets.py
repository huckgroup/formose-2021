'''
Example compound concentration timecourses. Figure 3A.
'''

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

experiment_code = 'FRN089B'
compound_number_sel = [3,9,13,21,18]
compound_numbering = {}
with open(repository_dir/'COMPOUND_INFO/compound_numbering.txt', 'r') as f:
    for line in f:
        ins = line.strip('\n').split(',')
        compound_numbering[int(ins[1])] = ins[0]

loaded_reports = []
path = repository_dir/'DATA/DATA_REPORTS/{}'.format(experiment_code)
for file in os.listdir(path):
    if file.endswith('csv'):
        file_path = '{}/{}'.format(path,file)
        loaded_reports.append(Classes.DataReport(file = file_path))

data = Classes.DataSet(data_reports = loaded_reports)

for c,n in enumerate(compound_number_sel):
    fig, ax = plt.subplots(figsize = (1.5/2.54,0.5/2.54))

    name = compound_numbering[n]

    x, y = data.get_entry(name + '/ M')

    trace = ax.plot(x,1000*y, '-o', markersize = 0.5, linewidth  = 0.5,
            c = info_params.colour_assignments[name])
    trace[0].set_clip_on(False)

    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3, prune='both',
                                    min_n_ticks  = 3))
    ax.tick_params(axis='both', which='major', labelsize = 3,
                       length = 2, pad = 2)

    #ax.set_xlabel('time/s', fontsize = 9)
    #ax.set_ylabel('Concentration/ mM', fontsize = 9)
    plt.savefig(repository_dir/f'PLOTS/Figure_3_{n}_inset.png', dpi = 600)
    plt.savefig(repository_dir/f'PLOTS/Figure_3_{n}_inset.svg')
    #plt.show()
    plt.close()

