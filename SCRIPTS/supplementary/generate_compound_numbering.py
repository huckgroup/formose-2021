import NorthNet
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

import sys
from pathlib import Path
# add the SCRIPTS directory to the system path
# so that its contents can be imported
script_dir = Path(__file__).parents[1].as_posix()
sys.path.append(script_dir)
# get the repository directory for file output
repository_dir = Path(__file__).parents[2]
# set paths to files
data_folder = repository_dir/'DATA'
derived_parameters_dir = data_folder/'DERIVED_PARAMETERS'
plot_folder = repository_dir/'PLOTS'
report_directory = data_folder/'DATA_REPORTS'
exp_info_dir = repository_dir/"EXPERIMENT_INFO/Experiment_parameters.csv"
reaction_list_directory = Path(repository_dir/'REACTION_LISTS')


df = pd.read_csv(derived_parameters_dir/'AverageData.csv', index_col = 0)
df = df.dropna(axis = 1)

df = df.loc[:, (df != 0).any(axis=0)]

compound_numbering = {str(x.split('/')[0]):i for i,x in enumerate(df.columns,1)}


with open(repository_dir/'COMPOUND_INFO/compound_numbering.txt', 'w') as f:
    for c in compound_numbering:
        f.write('{},{}\n'.format(c,compound_numbering[c]))

molecules = []
for x in compound_numbering:
    mol = Chem.MolFromSmiles(x)
    if x == 'C=O':
        mol = Chem.AddHs(mol)

    molecules.append(mol)

for m in molecules:
    tmp=AllChem.Compute2DCoords(m)

img=Draw.MolsToGridImage(molecules,molsPerRow=10,
                        subImgSize=(100,100),
                        legends=[str(compound_numbering[c]) for c in compound_numbering],
                        useSVG = False)

svg=Draw.MolsToGridImage(molecules,molsPerRow=10,
                        subImgSize=(300,300),
                        legends=[str(compound_numbering[c]) for c in compound_numbering],
                        useSVG = True)
img.save(repository_dir/'PLOTS/compound_numbering.png')

with open(repository_dir/'PLOTS/compound_numbering.svg', 'w') as f:
    f.write(svg)


from NorthNet import info_params
import matplotlib.pyplot as plt
print(len(info_params.colour_assignments))
fig,ax = plt.subplots(ncols = 10, nrows = 6, figsize = (10,10))
ax = ax.flatten()
for i,x in enumerate(compound_numbering):
    ax[i].scatter(1,1, c = info_params.colour_assignments[x], s = 10000)
    ax[i].annotate(info_params.smiles_to_names[x], xy = (1,1),
                ha = 'center', va = 'center', fontsize = 7,
                rotation = 60, fontweight = 'bold', color = 'k')
    ax[i].set_title(compound_numbering[x])
for a in ax:
    a.set_axis_off()

plt.savefig(repository_dir/'PLOTS/colour_key.png')
plt.close()
