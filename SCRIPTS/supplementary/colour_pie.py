import pandas as pd
import NorthNet
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np

df = pd.read_csv('information_sources/AverageData.csv', index_col = 0)
df = df.dropna(axis = 1)

df = df.loc[:, (df != 0).any(axis=0)]

compound_numbering = {str(x.split('/')[0]):i for i,x in enumerate(df.columns,1)}

from NorthNet import info_params
import matplotlib.pyplot as plt


vals = np.full(len(compound_numbering),1)
clrs = [info_params.colour_assignments[x] for x in compound_numbering]
clrs.reverse()
fig,ax = plt.subplots(figsize = (10,10))
wedges, texts = ax.pie(vals/vals.max(), colors = clrs,
                    startangle=90,
                    wedgeprops={"edgecolor":"k",
                                'linewidth': 0.5,
                                'antialiased': True},
                                )
fig.tight_layout()
plt.savefig('plots_for_paper/pie_key.png', dpi = 600)
plt.close()
