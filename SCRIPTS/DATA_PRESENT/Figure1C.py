import os
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from ChromProcess import Classes
import matplotlib.patches as patches

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from ChromProcess import file_import
from ChromProcess import info_params
from ChromProcess import file_import as f_i

plot_folder = repository_dir/'PLOTS'

C_chain_regions = {"tetradecane" : [6.5,6.9],
                   "C3": [5.2,6.5],
                   "C4": [7.5,8.8],
                   "C5": [9.3,12],
                   "C6": [12,18]}

C_chain_colors = { "tetradecane" : "#000000",
                   "C3": "#f09c08",
                   "C4": "#2738e7",
                   "C5": "#cb340b",
                   "C6": "#30bd37",
                   "C7": "#592387"}

# Choose the experiment code (see Data_information.csv)
exp_name = 'FRN089B'


# importing information
storage_stem = Path(r'/Users/williamrobinson/documents/nijmegen/PrebioticDatabase')

calib_file_path = '/Users/williamrobinson/Documents/Nijmegen/PrebioticDatabase/Analysis_Information/GCMS/2020_03_16_GCMS_Calibrations.csv'
calib = Classes.Instrument_Calibration(file = calib_file_path)

Path_file = storage_stem/'Data_information/Data_information.csv'
exp_paths = Classes.DataPaths(Path_file)

figure_width = 4.7/2.54
figure_height = 5.67/2.54

experiment = exp_paths.exp_code_path[exp_name]
data_type = experiment.data_type
exp_path = experiment.path

mod_bd_file = storage_stem/'Data/GCMS/FRN/{}/{}_local_assignments.csv'.format(exp_name,exp_name)
modified_bounds = file_import.read_local_assignments(mod_bd_file)
calib.modify_boundaries(modified_bounds)

# State directory in which to store results
store_folder = Path(storage_stem/'Data'/data_type/'FRN'/exp_name/'Chromatograms')


fig,ax = plt.subplots(figsize = (9.05/2.54,6.1/2.54))
time, signal = f_i.load_chromatogram_csv(store_folder/'FRN089_036.csv')
ax.plot(time, signal, c = "#000000", linewidth = 1)

# Create a Rectangle patch
for c in C_chain_regions:
    xy = [C_chain_regions[c][0],0]
    width = C_chain_regions[c][1] - C_chain_regions[c][0]
    height = np.amax(signal)
    rect = patches.Rectangle(xy, width, height, facecolor = C_chain_colors[c],
                            alpha = 0.5, zorder = 0)
    # Add the patch to the Axes
    ax.add_patch(rect)

bbox_props = dict( fc="none", ec='none')
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

ax.set_xlim(5.2,17)
ax.set_yticklabels([])
ax.set_xlabel('Retention time/ min.', fontsize = 9)
ax.set_ylabel('Total ion count', fontsize = 9)
ax.tick_params(axis='both', which='major', labelsize = 7, length = 1.75)
ax.set_position([0.1,0.15,0.87,0.8])

plt.savefig(plot_folder/'Figure1C.png', dpi = 600)
plt.savefig(plot_folder/'Figure1C.svg', dpi = 600)
