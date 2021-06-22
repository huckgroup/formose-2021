import os
import copy
import numpy as np
from pathlib import Path

from NorthNet import Classes
from NorthNet import info_params
from NorthNet.file_loads import data_loads, info_loads
from NorthNet.data_processing import dataset_operations as d_ops

'''Loading data and defining parameters.'''
base_directory = Path('/Users/williamrobinson/documents/nijmegen')
#base_directory = Path(r'C:\Users\willi\Documents')
directory = base_directory/'Dynamic_environment_project'
plot_folder = directory/'Paper_plots'
save_data_folder = directory/'DataAnalysis'
save_data_folder_2 = directory/'DataAnalysisAllPoints'

def load_data_files_from_folder(path):
	loaded_reports = []
	for file in os.listdir(path):
		if file.endswith('csv'):
			file_path = '{}/{}'.format(path,file)
			loaded_reports.append(Classes.DataReport(file = file_path))
	return loaded_reports

def load_data_set(path):

	loaded_reports = load_data_files_from_folder(path)

	if 'time' in loaded_reports[0].series_unit:
		t_dep = True
	else:
		t_dep = False

	if len(loaded_reports) == 1:
		return loaded_reports[0]
	else:
		data_set = d_ops.combine_datasets(loaded_reports, time_dependent = t_dep)
		return data_set

def load_all_data_sets(exp_info):
	data_sets = {}
	for e in exp_info:
	    loaded_reports = []
	    for file in os.listdir(exp_info[e].path):
	        if file.endswith('csv'):
	            file_path = '{}/{}'.format(exp_info[e].path,file)
	            loaded_reports.append(Classes.DataReport(file = file_path))
	    if 'time' in loaded_reports[0].series_unit:
	        t_dep = True
	    else:
	        t_dep = False

	    data_sets[e] = d_ops.combine_datasets(loaded_reports, time_dependent = t_dep)

	return data_sets

# list of headers for the array columns
prime_header = [x+'/ M' for x in info_params.smiles_to_names]
# list of compound colours.
clrs  = [info_params.colour_assignments[x] for x in info_params.smiles_to_names]

report_directory = base_directory/'safestore_DEP'
# import experiment information
exp_info = info_loads.import_Experiment_information(report_directory/"Experiment_parameters.csv")

# import data
experiment_averages = data_loads.load_exp_compound_file('information_sources/AverageData.csv', prime_header)
averages_errors = data_loads.load_exp_compound_file('information_sources/Average_errors.csv', prime_header)
experiment_amplitudes = data_loads.load_exp_compound_file('information_sources/AmplitudeData.csv', prime_header)

# find the carbon-containing reactants for each experiment
# and remove them from the data.
modifications = {x:[] for x in experiment_averages}
carbon_inputs = {x:[] for x in experiment_averages}
for v in exp_info:
	for p in exp_info[v].parameters:
		col_name = p
		if 'temperature' in p:
			continue
		if 'C' in p and not 'Ca' in p and exp_info[v].parameters[p] > 0.0:
			tag = p.split('/')[0][1:-1] + '/ M'
			if tag in prime_header:
				i = prime_header.index(tag)
				modifications[v].append(i)
				carbon_inputs[v].append(p.split('/')[0][1:-1])

modified_averages = copy.deepcopy(experiment_averages)
for m in modifications:
	modified_averages[m][modifications[m]] = 0.0

'''Removing columns containing zeros and reactants from analysis matrix.'''
data = np.zeros((len(modified_averages), len(modified_averages[[*modified_averages][0]])))
for c,v in enumerate(modified_averages):
	data[c] = modified_averages[v]

errors = np.zeros((len(modified_averages), len(modified_averages[[*modified_averages][0]])))
for c,v in enumerate(modified_averages):
	errors[c] = averages_errors[v]

idx = np.argwhere(np.all(data[..., :] == 0, axis=0))
idx = np.vstack((idx,prime_header.index('O=CCO/ M')))
idx = np.vstack((idx,prime_header.index('O=C[C@H](O)CO/ M')))
data = np.delete(data, idx, axis = 1)
errors = np.delete(errors, idx, axis = 1)

# truncate the header to match the data matrix.
header = [h for i,h in enumerate(prime_header) if i not in idx]
# truncate the compound colour list to match the data matrix.
clrs = [info_params.colour_assignments[h.split('/')[0]] for h in header]
# Create a list of colours for cluster labelling
colourings = info_params.cluster_colour_map[:]
colourings.extend(info_params.cluster_grayscale)

'''Loading conditions.'''
condition_names = [*exp_info['FRN071A'].parameters]
conditions = np.zeros((len(exp_info), len(exp_info['FRN071A'].parameters)))
for x,e in enumerate(exp_info):
	for y,v in enumerate(exp_info[e].parameters):
		if v != 'Reaction_entry':
			conditions[x,y] = exp_info[e].parameters[v]

# removing input species from amplitudes
modified_amplitudes = copy.deepcopy(experiment_amplitudes)
for m in modifications:
	modified_amplitudes[m][modifications[m]] = 0.0

amplitudes = np.zeros((len(modified_averages), len(modified_amplitudes[[*modified_averages][0]])))
for c,v in enumerate(modified_amplitudes):
	amplitudes[c] = modified_amplitudes[v]

idx = np.argwhere(np.all(amplitudes[..., :] == 0, axis=0))
amplitudes = np.delete(amplitudes, idx, axis = 1)
