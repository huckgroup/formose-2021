import numpy as np
from pathlib import Path
from scipy.signal import find_peaks

from NorthNet import Classes
from NorthNet import info_params
from NorthNet.file_loads import info_loads

import __init__
from helpers import parameter_determination_helper as params
from helpers.loading_helper import load_data_files_from_folder

base_directory = Path('/Users/williamrobinson/Documents/Nijmegen')
#base_directory = Path(r'C:\Users\willi\Documents')
report_directory = base_directory/'safestore_DEP'
cleaned_output_directory = report_directory/'Cleaned_data_reports'
exp_info = info_loads.import_Experiment_information(report_directory/"Experiment_parameters.csv")
exp_list = [*exp_info]
header = [x+'/ M' for x in info_params.smiles_to_names]
# load in data sets
data_sets = {}
for e in exp_info:
    path = cleaned_output_directory/e
    reports = load_data_files_from_folder(path)
    data_sets[e] = Classes.DataSet(data_reports = reports)

amps_out = np.zeros((len(exp_info),len(header)))
tlags_out = np.zeros((len(exp_info),len(header)))
averages_out = np.zeros((len(exp_info),len(header)))

for c1,d in enumerate(data_sets):

    if exp_info[d].modulation == 'None':
        time_lag_dict = {dep:0.0 for dep in data_sets[d].compounds}
        amps_dict = {dep:0.0 for dep in data_sets[d].compounds}
        averages = {}
        for c in data_sets[d].compounds:
            X,Y = data_sets[d].get_entry(c)
            averages[c] = np.average(Y)
    else:
        period = exp_info[d].parameters['per[{}]/ s'.format(exp_info[d].modulation)]
        exp_freq = 1/period

        # find amplitudes
        # find averages
        averages = {}
        amps_dict  = {}
        for c in data_sets[d].compounds:

            X,Y = data_sets[d].get_entry(c)
            averages[c] = np.average(Y)

            x_fourier,y_fourier = params.fourier_transform(X,Y)

            x_fourier = x_fourier[1:]
            y_fourier = y_fourier[1:]

            ida = np.where((x_fourier > exp_freq/2)
                            &(x_fourier < exp_freq*2))[0]

            x_fourier_section = x_fourier[ida]
            y_fourier_section = y_fourier[ida]

            idx, _ = find_peaks(y_fourier_section,distance = 10)

            if len(idx) > 0:
                amps_dict[c] = y_fourier[idx]
            else:
                amps_dict[c] = 0.0

    for a in averages:
        idx = header.index(a)
        amps_out[c1,idx] = amps_dict[a]
        averages_out[c1,idx] = averages[a]

loc_header = [h for h in header]
with open('information_sources/AverageData.csv', 'w') as f:
    f.write(',')
    [f.write('{},'.format(h)) for h in loc_header]
    f.write('\n')
    for x in range(0,len(averages_out)):
        f.write('{},'.format(exp_list[x]))
        [f.write('{},'.format(y)) for y in averages_out[x]]
        f.write('\n')

with open('information_sources/AmplitudeData.csv', 'w') as f:
    f.write(',')
    [f.write('{},'.format(h)) for h in loc_header]
    f.write('\n')
    for x in range(0,len(amps_out)):
        f.write('{},'.format(exp_list[x]))
        [f.write('{},'.format(y)) for y in amps_out[x]]
        f.write('\n')
