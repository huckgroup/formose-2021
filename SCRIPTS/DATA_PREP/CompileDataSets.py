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

def write_to_file(dataset, filename = '', path = None):
    '''
    Parameters
    ----------
    filename: str
        name for file
    path: pathlib Path object
        Path to folder for file storage.
    '''

    import numpy as np

    if filename == '':
        filename = dataset.filename
    elif not filename.endswith('.csv'):
        filename = filename + '.csv'
    if path == None:
        fname = filename
    else:
        fname = path/filename

    with open(fname, 'w') as outfile:
        # writing experiment conditions to file
        dataset.write_conditions_header(outfile)
        # writing data
        sorted_keys = sorted([*dataset.data], key = lambda x:x.count('C'))

        outfile.write("start_data\n")

        p_header = [dataset.series_unit]

        out = np.array([dataset.series_values])

        for s in sorted_keys:
            p_header.append(s)
            out = np.vstack((out,dataset.data[s]))

        out = out.T
        [outfile.write("{},".format(x)) for x in p_header]

        outfile.write("\n")

        for x in range(0,len(out)):
            for y in range(0,len(out[x])):
                outfile.write("{},".format(out[x,y]))
            outfile.write("\n")

        outfile.write("end_data\n")

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
    write_to_file(d, path = output_directory)
