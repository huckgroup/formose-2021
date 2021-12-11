'''
Compiling the data sets ready for submission. Experiment data files are renamed
and recompiled to conform to an arabic numeral-based naming system.
'''
import sys
import pandas as pd
from pathlib import Path

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from helpers.loading_helper import load_data_files_from_folder
    
# set paths to files
data_folder = repository_dir/'DATA'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"

output_directory = repository_dir/'ConvertedData'

# load in experiment information.
exp_info = pd.read_csv(exp_info_dir, index_col = 0)
exp_list = list(exp_info.index)

# load in data sets
data_sets = []
for e in exp_info.index:
    path = report_directory/e
    reports = load_data_files_from_folder(path)
    data_sets.extend(reports)

output_lines = []
for d in data_sets:
    d.filename = d.filename.replace(d.experiment_code,
                        str(exp_info.loc[d.experiment_code,'Experiment_entry']))
    d.experiment_code = exp_info.loc[d.experiment_code,'Experiment_entry']

for d in data_sets:
    d.write_to_file(filename = d.filename, path = output_directory)
