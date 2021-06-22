import pandas as pd
import NorthNet
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

df = pd.read_csv('information_sources/AverageData.csv', index_col = 0)
df = df.dropna(axis = 1)

df = df.loc[:, (df != 0).any(axis=0)]

compound_numbering = {str(x.split('/')[0]):i for i,x in enumerate(df.columns,1)}


with open('information_sources/compound_numbering.txt', 'w') as f:
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
img.save('plots_for_paper/compound_numbering.png')

with open('plots_for_paper/compound_numbering.svg', 'w') as f:
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

plt.savefig('plots_for_paper/colour_key.png')
