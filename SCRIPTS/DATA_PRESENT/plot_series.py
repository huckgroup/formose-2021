import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from NorthNet import info_params

import __init__
from helpers import load_series
from helpers.loading_helper import data
from helpers.loading_helper import exp_info
from helpers.loading_helper import header

base_directory = Path('/Users/williamrobinson/Documents/Nijmegen')
# base_directory = Path(r'C:\Users\willi\Documents')

fname = base_directory/'safestore_DEP/Series_info.csv'
series_dict = load_series.load_series_sequences(fname)

with open('information_sources/compound_numbering.txt', 'r') as f:
    lines = f.readlines()

compound_numbering = {l.split(',')[0]:int(l.split(',')[1].strip('\n')) for l in lines}

series_sel = 'Temperature_series'
condition_sel = 'temperature/ oC'
x_name = 'Temperature/ $^o$C'
x_factor = 1
y_factor = 1000

series_sel = 'Formaldehyde_2_series'
condition_sel = '[C=O]/ M'
x_name = '[formaldehyde]/ mM'
x_factor = 1000
y_factor = 1000


# series_sel = 'Residence_time_2_series'
# condition_sel = 'residence_time/ s'
# x_name = 'residence_time/ s'
# x_factor = 1
# y_factor = 1000

idx = [[*exp_info].index(x) for x in series_dict[series_sel]]

series_stack = np.zeros((len(idx),len(data[0])))

for x in range(0,len(idx)):
    series_stack[x] = data[idx[x]]

series_progression = series_stack.T

series_x_values = [x_factor*exp_info[x].parameters[condition_sel]
                        for x in series_dict[series_sel]]

fig, ax = plt.subplots(figsize=(8.965/2.54,6.55/2.54))
ax.set_position([0.2, 0.2, 0.7, 0.7])

trace_labels = []
for x in range(0,len(series_progression)):
    if np.sum(series_progression[x]) == 0.0:
        continue
    species_name = header[x].split('/')[0]
    if species_name in info_params.colour_assignments:
        clr = info_params.colour_assignments[species_name]
    else:
        clr = 'k'
    ax.plot(series_x_values,series_progression[x]*y_factor, '-o',
            c = clr, markersize = 5, zorder = 100-species_name.count('C'))

    max_idx = np.where(series_progression[x] ==series_progression[x].max())[0]
    trace_labels.append((compound_numbering[species_name], clr))


    x_pos = series_x_values[max_idx[0]]
    y_pos = y_factor*series_progression[x,max_idx]
    ax.annotate(compound_numbering[species_name],
                xy = (x_pos, y_pos),
                zorder = 1000)

print(len(trace_labels))
label_x_positions = np.linspace(0,series_x_values[-1], num = len(trace_labels))
label_x_positions*=1.5
label_y_positions = np.full(len(trace_labels),y_factor*series_progression.max())

label_x_positions[13:] = label_x_positions[:len(label_x_positions)-13]
label_y_positions[13:] = y_factor*series_progression.max()*0.9
for c,t in enumerate(trace_labels):
    x_pos = label_x_positions[c]
    y_pos = label_y_positions[c]
    ax.annotate(t[0],
                xy = (x_pos, label_y_positions[c]),
                fontsize = 6,
                fontweight = 'bold',
                color = t[1],
                zorder = 1000,
                ha = 'center',
                va = 'center')

ax.set_xlabel(x_name)
ax.set_ylabel('Concentration/ mM')
plt.savefig('plots_for_paper/{}_series.png'.format(series_sel), dpi = 600)

plt.close()
