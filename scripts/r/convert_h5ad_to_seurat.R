# Script to convert the fully processed single-cell data from .h5ad to Seurat

# Set workdir ----
setwd("/workdir/dwm269/scCardiacOrganoid/")

# Load lib(s) -----
library(glue)

# Custom function built on `anndata-R`----
source("~/DWM_utils/sc_utils/adata2seu.R")

tmp.name <- "scCO_v1c"
if(!exists(glue("data/pyobjs/{tmp.name}.h5ad"))){
  scco.seu <- adata2seu(
    adata_path=glue("data/pyobjs/{tmp.name}.h5ad"), #7-4
    counts = "counts",
    data = NULL,
    project = "scCardiacOrganoid",
    meta.data = "obs",
    reductions = "obsm",
    graphs = NULL,#"obsp",
    save.rdata = glue("data/robjs/{tmp.name}.rds"),
    verbose=T
  )
}

# Cardiac/mesoderm subset
tmp.name <- "mesoderm_2a"
if(!exists(glue("data/pyobjs/{tmp.name}.h5ad"))){
  meso.seu <- adata2seu(
    adata_path=glue("data/pyobjs/{tmp.name}.h5ad"), #7-4
    counts = "counts",
    data = NULL,
    project = "scCardiacOrganoid",
    meta.data = "obs",
    reductions = "obsm",
    graphs = NULL,#"obsp",
    save.rdata = glue("data/robjs/{tmp.name}.rds"),
    verbose=T
  )
}

# Endoderm subset
tmp.name <- "endoderm_2a"
# if(!exists(glue("data/pyobjs/{tmp.name}.h5ad"))){
#   scco.seu <- adata2seu(
#     adata_path=glue("data/pyobjs/{tmp.name}.h5ad"), #7-4
#     counts = "counts",
#     data = NULL,
#     project = "scCardiacOrganoid",
#     meta.data = "obs",
#     reductions = "obsm",
#     graphs = NULL,#"obsp",
#     save.rdata = glue("data/robjs/{tmp.name}.rds"),
#     verbose=T
#   )
# }

#
# With SeuratDisk (VERY BUGGY, DON'T RECOMMEND) -----
# Convert AnnData to Seurat (.h5ad -> .h5seurat)
## https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
#Load libraries
if(FALSE){ #switch for this code block...
  library(Seurat, quietly = T)
  library(SeuratDisk, quietly = T)
  
  # Filename setup...
  RData_file = "data/robjs/ma_iPSC_v6names"
  h5ad_file = "data/pyobjs/ma_iPSC_v6names2"
  
  # RData_file = "robjs/scCO_v7"
  # h5ad_file = "pyobjs/scCO_v7"
  
  if(!file.exists(paste0(RData_file,".RData"))){
    if(!file.exists(paste0(h5ad_file,".h5seurat"))){
      Convert(
        source=paste0(h5ad_file,".h5ad"),
        assay = "RNA",
        dest = "h5seurat",
        overwrite = TRUE
      )
    }
    
    # Load in Seurat object
    scco.seu <- LoadH5Seurat(
      file = paste0(h5ad_file,".h5seurat"),
      assays=c("RNA")
    )
    
    save(
      scco.seu,
      file=paste0(RData_file,".RData")
    )
    
  }else(
    message(paste0(h5ad_file,".h5ad has already been converted!"))
  )
}
# anndata -> SingleCellExperiment -> Seurat ----
# https://theislab.github.io/zellkonverter/reference/index.html
if(FALSE){
  #TODO, but also probably not because this would just further the portability issues of single-cell data between R and python.
}
