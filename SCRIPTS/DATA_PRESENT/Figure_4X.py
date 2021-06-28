import sys
import numpy as np
import pandas as pd
import networkx as nx
from rdkit import Chem
from pathlib import Path
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]

from helpers import chem_info as info_params
from helpers import pathway_helpers as path_hlp
from helpers.loading_helper import get_carbon_inputs
from helpers.network_load_helper import load_network_from_reaction_list
from helpers.network_load_helper import load_reaction_list,convert_to_networkx
from NorthNet.network_manipulations.networkx_ops import coordinates as c_ops

data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')
network_file = repository_dir/'FORMOSE_REACTION/FullFormoseReaction.txt'

series_sel = 'Formaldehyde_paper_series'
condition_sel_x = '[C=O]/ M'
factor = 1000
