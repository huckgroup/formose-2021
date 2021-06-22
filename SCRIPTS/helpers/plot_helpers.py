def PCA_loading_plot(X, clrs, ax, components = [0,1]):
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt

    pca = PCA()
    pca.fit(X)

    component_loadings = pca.components_[components,:].T
    for x in range(0,len(X)):
        ax.arrow(0,0,component_loadings[x,0],component_loadings[x,1],
                 color = clrs[x],
                 alpha = 0.9,
                 linewidth = 2)

    circ = plt.Circle((0, 0), radius=0.5,
                    edgecolor='#000000', facecolor='None', alpha = 0.2)
    ax.add_patch(circ)
    ax.set_xlabel('Component {}'.format(components[0]+1))
    ax.set_ylabel('Component {}'.format(components[1]+1))
    ax.set_aspect('equal')

def plot_cluster_dash(data, clrs, fname, components = 10):
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.neighbors import KernelDensity

    corr_coeffs = np.corrcoef(data.T)

    fig = plt.figure(constrained_layout=True)
    outer_grid = fig.add_gridspec(2, 2, wspace=0, hspace=0)
    square = math.ceil(np.sqrt(len(data)))
    inner_grid = outer_grid[0, 0].subgridspec(int(square),int(square),
    										  wspace=0, hspace=0)
    axs = inner_grid.subplots()  # Create all subplots for the inner grid.
    if square == 1:
        for i in range(0,len(data)):
        	axs.pie(np.nan_to_num(data[i]/np.amax(data[i])), colors = clrs)
    else:
        ax_set = axs.flatten()
        [a.set_axis_off() for a in ax_set]
        for i in range(0,len(data)):
        	ax_set[i].pie(np.nan_to_num(data[i]/np.amax(data[i])), colors = clrs)
        	ax_set[i].set_axis_on()

    pca = PCA()
    pca.fit(np.nan_to_num(corr_coeffs))

    ax = fig.add_subplot(outer_grid[1,1])
    ax.plot(np.arange(1,components+1,1),np.cumsum(pca.explained_variance_ratio_[:components]),'o-' ,c='#000000')
    ax.set_axis_on()
    ax.set_ylabel('Explained variance')
    ax.set_xlabel('Number of components')
    ax.set_title('Explained variance')

    ax = fig.add_subplot(outer_grid[1,0])
    plot_specs = data.T
    for z in range(0,len(data[0])):
    	idx_2 = np.where(plot_specs[z] > 0.0)[0]
    	if len(idx_2) > 0:
    		bin_width = 0.0005
    		if len(idx_2) == 1:
    			bins = 1
    		else:
    			bins = np.arange(np.amin(plot_specs[z,idx_2]), np.amax(plot_specs[z,idx_2])+bin_width, bin_width)

    		x_ax = np.arange(0.0, data.max(), 0.00005)
    		kde = KernelDensity(kernel='gaussian',
    						bandwidth=0.0005).fit(plot_specs[z,idx_2].reshape(-1,1))
    		log_density = kde.score_samples(x_ax.reshape(-1,1))
    		ax.plot(x_ax, np.exp(log_density), c = clrs[z], zorder = 1)

    ax.set_ylabel('Probability')
    ax.set_xlabel('concentration/ M')
    ax.set_title('Kernel density estimation')
    ax.set_xlim(0,np.amax(data))

    ax = fig.add_subplot(outer_grid[0,1])

    PCA_loading_plot(np.nan_to_num(corr_coeffs), clrs, ax)
    ax.set_title('Loading plot')

    plt.savefig(fname, dpi = 600)
    plt.close()
