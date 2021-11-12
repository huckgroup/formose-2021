'''

'''
import sys
import pickle
import pandas as pd
import networkx as nx
from pathlib import Path
import matplotlib.pyplot as plt

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from NorthNet import Classes
from helpers import chem_info as info_params
from helpers.network_load_helper import convert_to_networkx

data_folder = repository_dir/'DATA'
determined_params_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')

exp_code = 'FRN088B'
reaction_file = reaction_list_directory/f'FRN088B_reaction_list.txt'

with open(reaction_file, 'r') as f:
    for line in f:
        print(line)
    lines = f.readlines()

reactions = [x.strip('/n') for x in lines]

print(lines)
