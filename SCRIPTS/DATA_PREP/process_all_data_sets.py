import sys
import numpy as np
import pandas as pd
from pathlib import Path

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from scipy.signal import find_peaks

from NorthNet import Classes
import helpers.chem_info as info_params

from helpers import parameter_determination_helper as params
from helpers.loading_helper import load_data_files_from_folder

# set paths to files
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

# load in experiment information.
exp_info = pd.read_csv(exp_info_dir, index_col = 0)
exp_list = list(exp_info.index)

# load in data sets
data_sets = {}
for e in exp_info.index:
    path = report_directory/e
    reports = load_data_files_from_folder(path)
    data_sets[e] = Classes.DataSet(data_reports = reports)

compound_header = [x+'/ M' for x in info_params.smiles_to_names]

# create storage arrays for the data
amps_out = np.zeros((len(exp_info),len(compound_header)))
averages_out = np.zeros((len(exp_info),len(compound_header)))
import matplotlib.pyplot as plt
for c1,d in enumerate(data_sets):

    if exp_info.loc[d, 'Modulated_component'] == 'None':
        amps_dict = {dep:0.0 for dep in data_sets[d].compounds}
        averages = {}
        for c in data_sets[d].compounds:
            comp_idx = compound_header.index(c)
            X,Y = data_sets[d].get_entry(c)
            averages_out[c1,comp_idx] = np.average(Y)
    else:
        print(d)
        period = exp_info.loc[d,'per[{}]/ s'.format(
                                        exp_info.loc[d,'Modulated_component'])]
        exp_freq = 1/period

        for c in data_sets[d].compounds:

            comp_idx = compound_header.index(c)
            X,Y = data_sets[d].get_entry(c)

            averages_out[c1,comp_idx] = np.average(Y)

            x_fourier,y_fourier = params.fourier_transform(X,Y)

            x_fourier = x_fourier[1:]
            y_fourier = y_fourier[1:]

            ida = np.where((x_fourier > exp_freq/2)
                            &(x_fourier < exp_freq*2))[0]

            x_fourier_section = x_fourier[ida]
            y_fourier_section = y_fourier[ida]

            idx, _ = find_peaks(y_fourier_section,distance = 10)

            if len(idx) > 0:
                amps_out[c1,comp_idx] = y_fourier[idx]

with open(derived_parameters_dir/'AverageData.csv', 'w') as f:
    f.write(',')
    [f.write('{},'.format(h)) for h in compound_header]
    f.write('\n')
    for x in range(0,len(averages_out)):
        f.write('{},'.format(exp_list[x]))
        [f.write('{},'.format(y)) for y in averages_out[x]]
        f.write('\n')

with open(derived_parameters_dir/'AmplitudeData.csv', 'w') as f:
    f.write(',')
    [f.write('{},'.format(h)) for h in compound_header]
    f.write('\n')
    for x in range(0,len(amps_out)):
        f.write('{},'.format(exp_list[x]))
        [f.write('{},'.format(y)) for y in amps_out[x]]
        f.write('\n')
