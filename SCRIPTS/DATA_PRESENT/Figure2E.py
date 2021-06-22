import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from NorthNet import info_params

import __init__
from helpers import load_series
from helpers.loading_helper import data
from helpers.loading_helper import header
from helpers.loading_helper import exp_info

base_directory = Path('/Users/williamrobinson/Documents/Nijmegen')
# base_directory = Path(r'C:\Users\willi\Documents')

with open('information_sources/compound_numbering.txt', 'r') as f:
    lines = f.readlines()

compound_numbering = {l.split(',')[0]:int(l.split(',')[1].strip('\n')) for l in lines}

fname = base_directory/'safestore_DEP/Series_info.csv'
series_dict = load_series.load_series_sequences(fname)

series_sel = 'Ca_OH_grid_series'
condition_sel_x = '[NaOH]/ M'
condition_sel_y = '[CaCl2]/ M'

idx = [[*exp_info].index(x) for x in series_dict[series_sel]]

series_stack = np.zeros((len(idx),len(data[0])))

for x in range(0,len(idx)):
    series_stack[x] = data[idx[x]]

factor = 1000
series_x_values = [factor*exp_info[x].parameters[condition_sel_x]
                        for x in series_dict[series_sel]]

series_y_values = [factor*exp_info[x].parameters[condition_sel_y]
                        for x in series_dict[series_sel]]

subplot_width = 8
subplot_height = 8
fig, ax = plt.subplots(nrows = 1, figsize = (8.965/2.54,6.55/2.54))
ax.scatter(series_x_values, series_y_values, alpha = 0.0)
for x in range(0,len(series_stack)):
    L = series_x_values[x] - subplot_width/2
    B = series_y_values[x] - subplot_height/2
    W = subplot_width
    H = subplot_height
    axin = ax.inset_axes([L,B,W,H], transform=ax.transData)

    axin.pie(series_stack[x]/series_stack[x].max(),
             colors = [info_params.colour_assignments[x.split('/')[0]]
             for x in header],
             radius = 1.9,
             # labels = [compound_numbering[x.split('/')[0]] for x in
             #        header]
             )

ylim = ax.get_ylim()
xlim = ax.get_xlim()
ax.set_position([0.2, 0.2, 0.7, 0.7])
ax.set_ylim(ylim[0]-subplot_width/2, ylim[1]+subplot_width/2)
ax.set_xlim(xlim[0]-subplot_width/2, xlim[1]+subplot_width/2)
ax.set_xlabel('[NaOH]/ mM')
ax.set_ylabel('[CaCl$_2$]/ mM')

plt.savefig('plots_for_paper/Ca_OH_grid.png', dpi = 600)
plt.close()
