'''
Function for loading files.
'''
import os
import sys
import copy
import numpy as np
from pathlib import Path

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

import helpers.chem_info as info_params

from NorthNet import Classes

def load_data_files_from_folder(path):
    '''
    Load all data reports in a folder
    '''

    loaded_reports = []
    for file in os.listdir(path):
        if file.endswith('csv'):
            file_path = '{}/{}'.format(path,file)
            loaded_reports.append(Classes.DataReport(file = file_path))
    return loaded_reports

def get_carbon_inputs(exp_info, compounds_from_data):
    '''
    exp_info: NorthNet ExperimentInformation

    carbon_inputs: dict
    '''

    carbon_inputs = {x:[] for x in exp_info.index}
    for v in exp_info.index:
        for p in exp_info.columns:
            col_name = p
            if 'temperature' in p:
                continue
            if 'C' in p and not 'Ca' in p and exp_info.loc[v,p] > 0.0:
                tag = p.split('/')[0][1:-1] + '/ M'
                if tag in compounds_from_data:
                    carbon_inputs[v].append(tag)
    return carbon_inputs
