import os
import shutil
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from NorthNet import info_params
from NorthNet.file_loads import info_loads

import __init__
from helpers.loading_helper import load_all_data_sets
quit()
new_root = Path('/Users/williamrobinson/Documents/Nijmegen/safestore_DEP')

directory = Path(r'/Users/williamrobinson/Documents/Nijmegen/Dynamic_environment_project')
exp_info = info_loads.import_Experiment_information(directory/"Experiment_parameters.csv")
exp_info = {e:exp_info[e] for e in exp_info if 'FRN068' not in exp_info[e].name}

folder_translations = {e:exp_info[e].path.replace(str(directory)+'/','') for e in exp_info}

for e in exp_info:
    print(e)
    new_path = new_root/exp_info[e].path.replace(str(directory)+'/','').replace(folder_translations[e], e)
    print(new_path)
    file_to_copy = ''
    for file in os.listdir(exp_info[e].path):
        if file.endswith('.csv'):
            fname = '{}/{}'.format(exp_info[e].path,file)
            target_name = new_path/file
            with open(fname, 'r', encoding = 'latin-1') as f:
                type = ['','']
                for line in f:
                    if 'Chromatography_method' in line:
                        type = line.strip('/n').split(',')
                        break
            if type[1] == 'GCMS':
                file_to_copy = fname
                destination = target_name
                break
            elif type[1] == '' and 'GCMS' in file:
                file_to_copy = fname
                destination = target_name
                break

    if file_to_copy == '':
        print('passed {}'.format(e))
    else:
        if not os.path.exists(destination):
            print(file_to_copy, '-->', destination)
            print()
            shutil.copyfile(file_to_copy, destination)
