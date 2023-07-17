# Utility functions for use with scanpy
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib.path import Path
from scipy.sparse import issparse
from typing import Union
import scipy.sparse as sp


def npcs(
        ADATA, 
        var_perc=0.95, 
        reduction="pca"
    ):
    """
    Calculate the number of Principal Components (PCs) that contain a certain proportion of the variance.
    
    Parameters:
    -------
    ADATA -- An AnnData object. Must have already been processed with PCA and contain a 'pca' entry in its 'obsm' field.
    var_perc -- A float indicating the proportion of variance to be covered by the selected PCs. Default is 0.95.
    reduction -- A string indicating the type of dimensionality reduction to use. Default is 'pca'.
    
    Returns:
    -------
    n_pcs -- The number of PCs needed to cover the specified proportion of the variance. If the specified 'reduction' is not found, returns None.
    """
    from numpy import sum, var
    get_var = lambda i: var(ADATA.obsm[reduction][:,i])

    if ADATA.obsm[reduction] is None:
        print(f"Reduction '{reduction}', not found!")
        return None
    else:
        var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
        var_cut = var_perc * sum(var_tmp)
        n_pcs = 0
        var_sum = 0
        while var_sum<var_cut and n_pcs<ADATA.obsm[reduction].shape[1]-1:
            var_sum = var_sum + var_tmp[n_pcs]
            n_pcs = n_pcs + 1

        return(n_pcs)


# Reorder a reduction by decreasing % variance
def reorder_reduction(
        ADATA, 
        reduction="pca",
        verbose=False
    ):
    """
    Re-order a dimensions of a reduction by decreasing % variance.
    
    Parameters:
    -------
    ADATA -- An AnnData object. Must have already been processed with PCA and contain a 'pca' entry in its 'obsm' field.
    reduction -- A string indicating the type of dimensionality reduction to use. Default is 'pca'.
    verbose -- A boolean to indicate whether to print the variance for each dimension. Default is False.
    
    This function doesn't return anything, but it modifies the AnnData object in place, re-ordering the dimensions
    of the specified reduction in the 'obsm' field based on their variance (in decreasing order).
    """
    from numpy import var, argsort

    if reduction in ADATA.obsm:
        get_var = lambda i: var(ADATA.obsm[reduction][:,i])
        var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
        if verbose:
            print("Reduction variance by dimension:")
            print(var_tmp)

        pc_order = argsort(var_tmp)[::-1]
        ADATA.obsm[reduction] = ADATA.obsm[reduction][:,pc_order]
    else:
        print(f"The reduction '{reduction}' was not found...")


# Read in a list of gene lists from .csv (each column is a gene list)
def read_csv_to_dict(
        filename, 
        names2check=""
    ):
    """
    Read in a list of gene lists from .csv (each column is a gene list).
    
    Parameters:
    -------
    filename -- A string specifying the location of the csv file.
    names2check -- A list of gene names to filter the dictionary by. If not specified, all gene names are included. 

    Returns:
    -------
    dict_out -- A dictionary with column headers from the csv file as keys, and lists of genes as values. If names2check is specified, only genes in names2check are included in the lists.
    """
    import csv

    # Open the CSV file
    with open(filename, 'r') as file:

        # Create a CSV reader object
        reader = csv.reader(file)

        # Read the first row as header
        header = next(reader)

        # Create an empty dictionary to store the columns
        dict_out = {col: [] for col in header}

        # Loop through each row in the CSV file
        for row in reader:

            # Loop through each column in the row
            for col, value in zip(header, row):

                # Add the value to the corresponding column in the dictionary
                if value: # skip empty strings
                    dict_out[col].append(value)

    # Filter out unwanted entries based on the list in `names2check`
    if len(names2check) > 1:
        for KEY in dict_out.keys():
            dict_out[KEY] = [k for k in dict_out[KEY] if k in names2check]

    # Return the dictionary
    return dict_out


# Function to export DGEA results to a .csv file
def export_dgea_to_csv(
    adata: ad.AnnData,
    dgea_name,
    n_features,
    csv_out,
    axis=0,
    wide=False
):
    """
    Function to export DGEA (Differential Gene Expression Analysis) results to a .csv file.
    
    Parameters:
    -------
    adata -- An AnnData object which stores the gene expression data and metadata.
    dgea_name -- A string specifying the name of the DGEA results in adata.uns[].
    n_features -- An integer specifying the number of top features to be included in the exported .csv file.
    csv_out -- A string specifying the path to the output .csv file.
    axis -- An integer specifying how to write results for each group. If 1, results are written horizontally. If 0, results are written vertically. Default is 0.
    wide -- A boolean specifying the format of the .csv file. If True, the .csv file will have one row per feature, and each column will be a group. If False, the .csv file will have one row per group-feature combination. Default is False.

    Returns:
    -------
    The function doesn't return anything but writes the DGEA results to a .csv file.
    """
    import pandas as pd
    import scanpy as sc

    result = adata.uns[dgea_name]
    groups = result['names'].dtype.names

    if wide:
        celltype_markers = pd.DataFrame(
            {group + '_' + key[:-1]: result[key][group]
            for group in groups for key in ['names', 'logfoldchanges','pvals']}).head(n_features)
        celltype_markers.to_csv(csv_out, index=False)
    else:
        marker_list = list()
        for group in adata.uns[dgea_name]['names'].dtype.names:
            markers = sc.get.rank_genes_groups_df(adata, key=dgea_name, group = group).head(n_features)
            markers['celltypes'] = group
            marker_list.append(markers)
        
        celltype_markers = pd.concat(
            marker_list, 
            axis=axis
        )
        celltype_markers.to_csv(
            csv_out,
            index=False
        )


# Function to convert feature names 
def convert_feature_names(
        adata: ad.AnnData,
        gtf_info: pd.DataFrame, 
        from_col: str='GENEID',
        to_col: str='GeneSymbol',
        inplace: bool=True,
        verbose: bool=True
) -> ad.AnnData:
    """
    Function to convert feature names in an AnnData object using mapping provided in a DataFrame.
    
    Parameters:
    -------
    adata       -- An AnnData object which stores the gene expression data and metadata.
    gtf_info    -- A DataFrame containing the mapping from one set of feature names to another.
    from_col    -- A string specifying the column in gtf_info to be mapped from. Default is 'GENEID'.
    to_col      -- A string specifying the column in gtf_info to be mapped to. Default is 'GeneSymbol'.
    inplace     -- A boolean specifying whether to perform the conversion inplace or return a new AnnData object. Default is True.
    verbose     -- A boolean specifying whether to print progress information. Default is True.

    Returns:
    -------
    If inplace is False, returns a new AnnData object with converted feature names.
    """

    if not inplace:
        adata = adata.copy()

    # Filter gtf_info to keep only the gene names found in the anndata object
    if from_col not in gtf_info.columns:
        raise ValueError(f"Column {from_col} not found in gtf_info")
    else:
        gtf_info_filtered = gtf_info[gtf_info[from_col].isin(adata.var_names)]
    
    if verbose:
        num_found = len(gtf_info_filtered)
        num_total = len(adata.var_names)
        # fraction_found = num_found / num_total
        print(f"Fraction of adata.var_names found in gtf_info[{from_col}]: {num_found} out of {num_total}")
    
    gene_name_mapping = dict(zip(
        gtf_info[from_col], 
        gtf_info[to_col]
    ))
    
    adata.var[from_col] = adata.var_names
    adata.var[to_col] = adata.var[from_col].map(gene_name_mapping)

    adata.var.dropna(subset=[to_col], inplace=True)
    adata.var.reset_index(drop=True, inplace=True)

    # mask = adata.var_names.isin(adata.var[to_col].values)
    # adata = adata[:, mask].copy()
    adata.var_names = adata.var[to_col]
    adata.var_names_make_unique()

    if not inplace:
        return adata


# Remove cells with fewer than K neighbors within a distance D
# def spatial_singlet_filter(
#         adata: ad.AnnData,
#         basis="spatial",
#         D: int = 100, 
#         K: int = 10
#     ) -> ad.AnnData:
#     """
#     This function calculates the Euclidean distances between spatial coordinates, 
#     finds the number of spatial neighbors within a given distance (D), and filters out cells 
#     with fewer than a given number of neighbors (K).

#     Parameters
#     ----------
#     adata: sc.AnnData
#         An AnnData object. The function expects that this object contains spatial coordinates in adata.obsm['spatial'].
        
#     D: int, default=100
#         The maximum Euclidean distance to consider when determining spatial neighbors.
        
#     K: int, default=10
#         The minimum number of spatial neighbors required to keep a cell in the output AnnData object.

#     Returns
#     -------
#     adata: ad.AnnData
#         An updated AnnData object which only contains cells with more than K neighbors within a distance D. 
#         Also, this object contains a new field in .obsm: "spatial_distances", which contains the Euclidean 
#         distances between spatial coordinates, and a new field in .obs: "spatial_neighbors_{D}_true", which contains
#         the number of spatial neighbors for each cell within a distance D.
#     """
#     import scipy.spatial as scisp
#     import numpy as np
#     import anndata as ad

#     # Calculate Euclidean distances
#     adata.obsm[f"{basis}_distances"] = scisp.distance.squareform(scisp.distance.pdist(adata.obsm[basis]))
    
#     # Calculate number of spatial neighbors within distance D for each cell
#     adata.obs[f"{basis}_neighbors_{D}"] = np.sum(adata.obsm[f"{basis}_distances"] < D, axis=0)
    
#     # Filter out cells with fewer than K spatial neighbors within distance D
#     adata = adata[adata.obs[f"{basis}_neighbors_{D}"] > K,]
    
#     return adata


# Function to segment tissues based on spatial information
def segment_tissues(adata, threshold='auto', num_tissues=None, inplace=True, verbose=True):
    """
    Segment tissues based on spatial information.

    Parameters:
        adata (Anndata): Annotated data object containing spatial coordinates.
        threshold (float or str): Threshold distance for tissue segmentation. If 'auto', the threshold will be calculated based on the expected number of tissues.
        num_tissues (int): Expected number of tissues. Required if threshold='auto'.
        inplace (bool): If True, the tissue labels will be added as adata.obs[f"segment_{threshold}"]. If False, the function will return a modified copy of the input Anndata object.
        verbose (bool): If True, print progress and summary messages.

    Returns:
        Anndata: Annotated data object with tissue labels added as adata.obs[f"segment_{threshold}"].

    """
    if not inplace:
        adata = adata.copy()

    spatial_coords = adata.obsm['spatial']
    num_cells = len(spatial_coords)

    # Calculate the average distance between cells
    avg_distance = np.mean(np.linalg.norm(spatial_coords - np.mean(spatial_coords, axis=0), axis=1))

    # Calculate the threshold based on the expected number of tissues if 'auto' is selected
    if threshold == 'auto':
        if num_tissues is None:
            raise ValueError("num_tissues must be provided when threshold='auto'.")
        threshold = num_tissues * avg_distance

    tissue_labels = np.arange(num_cells)
    for i in range(num_cells):
        if verbose and i % (num_cells // 4) == 0:
            print(f"Processing cell {i+1}/{num_cells}")

        distances = np.linalg.norm(spatial_coords - spatial_coords[i], axis=1)
        similar_cells = np.where(distances <= threshold)[0]
        tissue_labels[similar_cells] = tissue_labels[i]

    unique_labels, label_counts = np.unique(tissue_labels, return_counts=True)
    label_mapping = {label: str(i) for i, label in enumerate(unique_labels)}
    tissue_labels = pd.Categorical(tissue_labels.astype(str))
    tissue_labels = tissue_labels.rename_categories(label_mapping)
    adata.obs[f"segment_{threshold}"] = tissue_labels

    if verbose:
        unique_labels, label_counts = np.unique(tissue_labels, return_counts=True)
        print(f"Segmentation completed using threshold={threshold}")
        print(f"Number of identified tissues: {len(unique_labels)}")
        print("Tissue labels summary:")
        for label, count in zip(unique_labels, label_counts):
            print(f"Tissue {label}: {count} cells")

    if not inplace:
        return adata


# Function to add biotype % values to AnnData object
def add_biotypes_pct(
    adata: ad.AnnData,
    biomart: Union[None, pd.DataFrame] = None, # DataFrame containing gene biotypes
    gene_colname: str = "GeneSymbol",
    biotype_colname: str = "Biotype",
    add_as: str = "obs", # how percent features should be added
    prefix: str = "pct.",
    scale: int = 100,
    verbose: bool = True
) -> ad.AnnData:
    """
    This function adds gene biotype percentage values to an AnnData object.
    
    Args:
        adata (AnnData): The AnnData object containing gene expression data.
        biomart (pd.DataFrame, optional): A DataFrame containing gene biotypes.
        gene_colname (str, optional): Column name in biomart DataFrame for gene identifiers. Default is "GeneSymbol".
        biotype_colname (str, optional): Column name in biomart DataFrame for biotype. Default is "Biotype".
        add_as (str, optional): Determines how percent features should be added. Default is "obs".
        prefix (str, optional): Prefix for column names added to the AnnData object. Default is "pct.".
        scale (int, optional): Determines the scaling for the percentage. Default is 100.
        verbose (bool, optional): Determines whether to print messages during function execution. Default is True.
        
    Returns:
        AnnData: The original AnnData object with added gene biotype percentage values.
    """

    if biomart is None:
        if verbose: print("Need a list of gene biotypes! Nothing done.")
        return adata

    if add_as == "var":
        if verbose: print("add_as='var' is not yet implemented")
        return adata

    if verbose: print(f"Adding gene biotype percentage values as {add_as} ...")

    biotypes = biomart[biotype_colname].unique()

    for biotype in biotypes:
        # Subset out current gene biotype
        tmp_mart = biomart[biomart[biotype_colname] == biotype]

        # Get unique gene names which are present in the adata object
        tmp_feat = tmp_mart[gene_colname][tmp_mart[gene_colname].isin(adata.var_names)].unique()

        if len(tmp_feat) == 0:
            if verbose: print(f"  No {biotype} genes found...")
        else:
            if add_as == "obs":
                col_name = prefix + biotype
                
                # Calculate the percentage
                if issparse(adata.X):
                    gene_pct = adata[:, tmp_feat].X.sum(axis=1) / adata.X.sum(axis=1)
                else:
                    gene_pct = adata[:, tmp_feat].X.sum(axis=1) / adata.X.sum(axis=1)
                
                # Add the percentage to the AnnData object
                adata.obs[col_name] = np.asarray(gene_pct).flatten()
                
                # Scale the data if necessary
                if scale == 1: # [0,1]
                    adata.obs[col_name] /= 100
                elif scale != 100: # [0,100]
                    if verbose: print("Given scale was not found. Scaling to 100...")
            else:
                if verbose: print("`add_as` option not found... Try again.")

    return adata


# Function to label cells within a region of interest (ROI) polygon
def label_roi_polygon(adata, roi_dict, metadata_col_name='roi', as_string=False):
    """
    Add a column to adata.obs indicating whether each cell is within the defined region of interest (ROI) polygon.

    Parameters:
        adata (AnnData): Anndata object containing spatial coordinates in adata.obsm['spatial'].
        roi_dict (dict): Dictionary containing the vertices of the ROI polygon.
                         Format: {'x': [x1, x2, x3, ...], 'y': [y1, y2, y3, ...]}
        metadata_col_name (str): Name of the new metadata column. Default is 'roi'.
        as_string (bool): If True, store the column values as strings ('True' or 'False').
                          If False, store the column values as booleans (True or False). Default is False.

    Returns:
        None
    """
    roi_path = Path(list(zip(roi_dict['x'], roi_dict['y'])))
    spatial_df = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'])
    roi_mask = roi_path.contains_points(spatial_df[['x', 'y']].values)

    if as_string:
        roi_mask = pd.Series(roi_mask).map({True: 'True', False: 'False'})

    adata.obs[metadata_col_name] = roi_mask.values


# Function to add a matrix as a new layer to an AnnData object
def add_mtx_as_layer(
    adata: ad.AnnData,
    mtx_path: str,
    row_names: np.array,
    col_names: np.array,
    layer_name: str,
    intersect: bool = False,
    inplace: bool = True,
):
    """
    Read an mtx file and add it as a new layer to an AnnData object.

    Args:
        adata (AnnData): AnnData object to which the matrix should be added as a new layer.
        mtx_path (str): Path to the mtx file.
        row_names (np.array): Array of row names in the mtx file.
        col_names (np.array): Array of column names in the mtx file.
        layer_name (str): Name of the layer to be added.
        intersect (bool, optional): Whether to consider only the intersection of observations.
            If False, the union of observations will be used. Default is False.
        inplace (bool, optional): Whether to modify the input AnnData object in place or create a new copy.
            If True, modifications are made to the input AnnData object. If False, a new AnnData object is created
            with the added layer. Default is True.

    Returns:
        None if inplace=True. Returns a new AnnData object with the added layer if inplace=False.
    """
    if not inplace:
        # Create a copy of the AnnData object
        adata = adata.copy()

    # Read the matrix file
    matrix = sp.load_npz(mtx_path)

    # Filter column names based on intersection or union with existing obs_names
    existing_obs_names = adata.obs_names
    if intersect:
        col_names = np.intersect1d(col_names, existing_obs_names)
    else:
        col_names = np.union1d(col_names, existing_obs_names)

    # Find the column indices corresponding to col_names in the matrix
    col_indices = np.where(np.isin(existing_obs_names, col_names))[0]

    # Filter the matrix and column names based on the selected columns
    matrix = matrix[:, col_indices]
    col_names = existing_obs_names[col_indices]

    # Create a new layer in the AnnData object
    adata.layers[layer_name] = matrix

    # Assign row and column names to the AnnData object
    adata.var_names = col_names
    adata.obs_names = row_names

    if not inplace:
        return adata

def top_n_genes(adata, n):
    """
    Identify the top N genes with the highest sum of expression values in an AnnData object.

    Parameters:
    adata (AnnData): An AnnData object storing the gene expression data and metadata.
    n (int): The number of top genes to return.

    Returns:
    list: A list of top N genes sorted by their sum of expression values.
    """
    sorted_genes = adata.var_names[np.argsort(adata.X.sum(axis=0))[::-1]]
    return sorted_genes[:n].tolist()