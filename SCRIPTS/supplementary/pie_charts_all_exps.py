from helpers.loading_helper import modified_averages, prime_header, clrs
from NorthNet import info_params
import numpy as np
import matplotlib.pyplot as plt

compound_numbers = {p.split('/ M')[0]:i for i,p in enumerate(prime_header)}
labels = np.arange(1,len(compound_numbers)+1, 1)

for idx in modified_averages:
    fig, ax = plt.subplots()

    inds = np.where(modified_averages[idx] > 0.0)[0]
    loc_labels = [labels[i] for i in inds]
    wedges, texts = ax.pie(modified_averages[idx][inds]/np.amax(modified_averages[idx][inds]),
        colors = [info_params.colour_assignments[x.split('/')[0]]
                for c,x in enumerate(prime_header) if c in inds],
        startangle=-40,
        labels = loc_labels,
        wedgeprops={"edgecolor":"k",
                    'linewidth': 1,
                    'antialiased': True}
        )

    ax.set_title(idx)
    plt.savefig('pie_plots/{}.png'.format(idx),dpi = 300)
    plt.close()


with open('pie_plots/compound_numbering.csv', 'w') as f:
    f.write('comopounds SMILES, number\n')
    for c in compound_numbers:
        f.write('{},{}\n'.format(c+1,compound_numbers[c]))
