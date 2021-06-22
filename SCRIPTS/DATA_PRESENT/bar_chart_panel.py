import numpy as np
import matplotlib.pyplot as plt

import __init__
from helpers.loading_helper import data, errors
from helpers.loading_helper import clrs
from helpers.loading_helper import header
from helpers.loading_helper import exp_info

exp_list = [*exp_info]
data_set_selections = ['FRN094C', 'FRN089A', 'FRN090B','FRN077C']
indices = [exp_list.index(s) for s in data_set_selections]

plot_stack = np.zeros((len(data_set_selections), len(data[0])))
error_stack = np.zeros((len(data_set_selections), len(data[0])))
for x,i in enumerate(indices):
    plot_stack[x] = data[i]
    error_stack[x] = errors[i]

idx = np.argwhere(np.all(plot_stack[..., :] == 0, axis=0))
plot_stack = np.delete(plot_stack, idx, axis = 1)
error_stack = np.delete(error_stack, idx, axis = 1)

comp_ax = [header[i] for i,c in enumerate(header) if i not in idx]
comp_clrs = [clrs[i] for i,c in enumerate(clrs) if i not in idx]
comp_numbers = [i for i,c in enumerate(header) if i not in idx]

class margins: left = 0.2;right = 0.05;top = 0.15;bottom = 0.13;
left = margins.left
width = 1-margins.left-margins.right
height = (1-margins.top-margins.bottom)/len(data_set_selections)
axis_label_font_size = 6
fig, ax = plt.subplots(nrows = len(data_set_selections),
                        figsize = (7.17/2.54,3.75/2.54))

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
    a.tick_params(axis = 'both', which = 'both', length = 2,
                        labelsize = 6)


fig.text((margins.left - margins.right + 1)/2, 0.04,
         'Compound number', ha='center',
         fontsize = axis_label_font_size)
fig.text(0.04, 0.53, 'log$_{10}$(concentration/ M)',
            va='center', rotation='vertical',
            fontsize = axis_label_font_size)

plt.savefig('plots_for_paper/Bar_charts.png', dpi = 600)
