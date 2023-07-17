# Plotting functions for use with scanpy
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from anndata import AnnData
# import anndata as ad

# Knee plot to quality check UMI counts for single-cell data
def knee_plot(
    ADATA,
    x_lim=[0, 20000],
    line_width=2,
    line_color="b",
    title="Knee plot",
    verbose=False
):
    import matplotlib.pyplot as plt
    expected_num_cells = 10000

    knee = np.sort((np.array(ADATA.X.sum(axis=1))).flatten())[::-1]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(
        knee,
        range(len(knee)),
        linewidth=line_width,
        color=line_color
    )
#     ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
#     ax.axhline(y=expected_num_cells, linewidth=3, color="k")

    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    ax.set_title(title)

    plt.xlim(x_lim)
    plt.grid(True, which="both")
    plt.show()

# Faceted plot for any embedding
# scanpy github issue reference- https://github.com/scverse/scanpy/issues/955
def facet_embedding(
    adata, 
    clust_key, 
    basis, 
    size=60, 
    frameon=False, 
    legend_loc=None, 
    **kwargs
):
    # import scanpy as sc

    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.embedding(
        tmp, 
        groups=tmp.obs[clust].cat.categories[1:].values, 
        color=adata.obs[clust_key].cat.categories.tolist(), 
        basis=basis,
        size=size, frameon=frameon, legend_loc=legend_loc, 
        **kwargs
    )


# scanpy version of seuListPlot - grids of plots for a single embedding
def plot_grid_of_embeddings(
        adata_dict, 
        color, 
        ncols=None, 
        figsize=(10, 10), 
        min_value=None, 
        na_color='lightgrey', 
        # cmap="plasma", 
        same_scale=False, 
        **kwargs
    ):
    """
    Plot a grid of embeddings for multiple Anndata objects.

    Parameters:
        adata_dict (dict): Dictionary of Anndata objects. The keys represent the plot titles and the values are the corresponding Anndata objects.
        color (list): List of feature names to use for coloring the embeddings.
        ncols (int): Number of columns in the grid. If None, it is set to the number of entries in adata_dict.
        figsize (tuple): Figure size (width, height).
        same_scale (bool): If True, use the same color scale for plots showing the same feature.
        **kwargs: Additional parameters to be passed to scanpy.pl.embedding.

    """

    if ncols is None:
        ncols = len(adata_dict)

    nrows = len(list(color))

    if len(color) >= 1:
        fig, axes = plt.subplots(
            nrows=nrows, 
            ncols=ncols, 
            figsize=figsize,
            # layout="constrained",
            sharex=True, 
            sharey='row'
        )
        axes = axes.flatten()
    else:
        print("Need to specify `color`")

    for i, feat in enumerate(color):
        if same_scale:
            # Find color scale range for current `feat`
            color_min = np.inf
            color_max = -np.inf
            for j, (title, adata) in enumerate(adata_dict.items()):
                # Find the min and max values of feat within each adata object
                if feat in adata.obs_names and adata.obs[feat].dtype.name != 'category': # metadata variables
                    feat_min = np.min(adata.obs[feat])
                    feat_max = np.max(adata.obs[feat])
                elif feat in adata.var_names: # gene expression
                    feat_min = np.min(adata[:, feat].X)
                    feat_max = np.max(adata[:, feat].X)
                else: # categorical variables
                    feat_min = 0 #None
                    feat_max = 42 #None
                color_min = min(color_min, feat_min)
                color_max = max(color_max, feat_max)

        for j, (title, adata) in enumerate(adata_dict.items()):
            ax = axes[i * ncols + j]
            ax.set_aspect('equal')

            # Set color limits for the plot
            vmin = color_min if same_scale else None
            vmax = color_max if same_scale else None
            
            # Set 
            if kwargs["cmap"] is not None and type(kwargs["cmap"]) == str:
                kwargs["cmap"] = cm.get_cmap(kwargs["cmap"])
                kwargs["cmap"].set_under(na_color)
                
            if min_value != None:
                vmin = min_value

            sc.pl.embedding(
                adata,
                color=feat,
                title=None,
                show=False,
                ax=ax,
                vmin=vmin,
                vmax=vmax,
                **kwargs
            )

            # y-axis
            if j == 0:
                ax.set_ylabel(feat, fontsize='medium', fontstyle='italic')
            else:
                ax.set_ylabel(None)

            # x-axis
            ax.set_xlabel(None)

            # plot titles
            if i == 0:
                ax.set_title(title)
            else:
                ax.set_title(None)

    plt.tight_layout()
    # plt.legend(loc='right')
    plt.show()


