'''
An example GC-MS chromatogram with annotation. Figure 1C of the main text.
'''
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from ChromProcess import Classes

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

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

# importing information
storage_stem = repository_dir/'DATA'/'CHROMATOGRAMS'

figure_width = 15/2.54
figure_height = 10/2.54
factor = 1e6

GCMS_chrom = Classes.Chromatogram(file = storage_stem/'FRN093_041_GCMS_chromatogram.csv')
HPLC_chrom = Classes.Chromatogram(file = storage_stem/'FRN093_041_HPLC_chromatogram.csv')
GCMS_inds = np.where((GCMS_chrom.time > 6) & (GCMS_chrom.time < 17))[0]
HPLC_inds = np.where((HPLC_chrom.time > 2) & (HPLC_chrom.time < 12))[0]

fig,ax = plt.subplots(nrows = 2, figsize = (figure_width, figure_height))

ax[0].plot(HPLC_chrom.time[HPLC_inds], HPLC_chrom.signal[HPLC_inds]/factor,
        c = "#000000", linewidth = 1.5)
ax[1].plot(GCMS_chrom.time[GCMS_inds], GCMS_chrom.signal[GCMS_inds]/factor,
        c = "#000000", linewidth = 1.5)

ax[0].set_ylabel('Intensity (364 nm)/ kV', fontsize = 10)
ax[1].set_xlabel('residence time/ min.', fontsize = 10)
ax[1].set_ylabel('Total ion counts/ 10$^6$', fontsize = 10)

ax[0].set_ylim(-0.01,0.2)
for a in ax:
    a.tick_params(labelsize = 8)

plt.savefig(plot_folder/'GCMS_HPLC_Comparison.png', dpi = 600)
plt.savefig(plot_folder/'GCMS_HPLC_Comparison.svg', dpi = 600)
plt.show()
