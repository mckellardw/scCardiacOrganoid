# Utility functions for use with scanpy

# Calculate the number of PCs that contain some proportion (default is 95%) of the variance
def npcs(
  ADATA,
  var_perc=0.95,
  reduction="pca"
):
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


# Re-order a dimensions of a reduction by decreasing % variance
def reorder_reduction(
    ADATA,
    reduction="pca",
    verbose=False
):
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
## Useful for plotting
def read_csv_to_dict(filename):
    import csv

    # Open the CSV file
    with open(filename, 'r') as file:

        # Create a CSV reader object
        reader = csv.reader(file)

        # Read the first row as header
        header = next(reader)

        # Create an empty dictionary to store the columns
        columns = {col: [] for col in header}

        # Loop through each row in the CSV file
        for row in reader:

            # Loop through each column in the row
            for col, value in zip(header, row):

                # Add the value to the corresponding column in the dictionary
                if value: # skip empty strings
                    columns[col].append(value)

    # Return the dictionary
    return columns


# Function to export DGEA results to a .csv file
def export_dgea_to_csv(
    adata, # adata
    dgea_name, # 'rank_gene_groups_clusters_leiden_0.5', name in adata.uns[]
    n_features,
    csv_out,
    wide=False
):
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
            markers["celltypes"] = group
            marker_list.append(markers)
        
        celltype_markers = pd.concat(marker_list, axis=1)
        celltype_markers.to_csv(csv_out, index=False)