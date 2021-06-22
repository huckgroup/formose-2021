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

fname = base_directory/'safestore_DEP/Series_info.csv'

series_sel = 'Formaldehyde_2_series'
condition_sel = 'residence_time/ s'

series_stack = np.zeros((len(exp_info),len(data[0])))

series_progression = data.T

series_x_values = [exp_info[x].parameters[condition_sel]
                        for x in exp_info]

fig, ax = plt.subplots(figsize=(8.71/2.54,6.55/2.54))
ax.set_position([0.2, 0.2, 0.7, 0.7])
for x in range(0,len(series_progression)):
    if np.sum(series_progression[x]) == 0.0:
        continue
    species_name = header[x].split('/')[0]
    if species_name in info_params.colour_assignments:
        clr = info_params.colour_assignments[species_name]
    else:
        clr = 'k'
    ax.plot(series_x_values,series_progression[x]*1000, 'o', c = clr)

for x in range(0,len(data)):
    ax.annotate([*exp_info][x], xy = (series_x_values[x], data[x].max()*1000), rotation = 90)

ax.set_xlabel('Residence time/ s')
ax.set_ylabel('Concentration/ mM')
plt.show()
# plt.savefig('plots_for_paper/{}_series.png'.format(series_sel), dpi = 600)
# plt.close()
