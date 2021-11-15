'''
Program for printing condition information
given an experiment code.
'''
import os
import sys
import copy
import pandas as pd
from pathlib import Path

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

conditions_file = repository_dir/'EXPERIMENT_INFO'/'Experiment_parameters.csv'

conditions = pd.read_csv(conditions_file, index_col=0)

# remove empty columns
conditions = conditions.dropna(axis = 1)

condition_names = conditions.columns
experiment_names = list(conditions.index)

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print('Provide experiment code')
        quit()

    exp_code = sys.argv[1]
    if exp_code in experiment_names:
        for c in condition_names:
            print(f'{c}: {conditions.loc[exp_code,c]}')
    else:
        print(f'{exp_code}: Invalid experiment code.')
        quit()
