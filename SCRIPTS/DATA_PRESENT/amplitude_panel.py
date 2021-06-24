# generate demonstrative amplitude diagrams for the formose reactions
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes

import helpers.chem_info as info_params
from helpers.loading_helper import load_data_set

from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
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


for n in compound_number_sel:
    name = compound_numbering[n]
    print(name)
    x, y = data.get_entry(name + '/ M')
    fig, ax = plt.subplots(figsize = (3/2.54,1.5/2.54))
    ax.plot(x,1000*y, '-o', markersize = 2.5, linewidth  = 1,
            c = info_params.colour_assignments[name])
    fig.tight_layout()
    ax.set_position([0.2, 0.1, 0.7, 0.8])
    ax.tick_params(which = 'both', axis = 'both',
                labelsize = 6, pad = 1, length = 1)
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_ylim(1000*(y.min() - y.min()*0.5), 1000*y.max()*1.1)

    # scalebar = AnchoredSizeBar(ax.transData,
    #                        6*60, '{} s'.format(6*60),
    #                        'lower right',
    #                        pad=0.1,
    #                        color='#000000',
    #                        frameon=False,
    #                        size_vertical=0.05,
    #                        fontproperties = fontprops)
    #
    # ax.add_artist(scalebar)
    ax.set_xticklabels([])
    plt.savefig(repository_dir/'PLOTS/{}_amps.png'.format(name), dpi = 600)
    plt.close()
