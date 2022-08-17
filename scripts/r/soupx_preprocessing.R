

# Load libraries  & helper functions ----
library(Matrix)
library(dplyr)
library(Seurat)
library(future)
library(parallel)
library(SoupX)


source("~/DWM_utils/sc_utils/seurat_helpers/seutils.R")

# Load metadata -----

meta <- read.csv("resources/metadata.csv")

# Helper function(s) ----
getClusterIDs <- function(toc, verbose=F){
  seu <- CreateSeuratObject(toc)
  seu <- seu %>%
    NormalizeData(verbose=verbose)%>%
    ScaleData(verbose=verbose) %>%
    FindVariableFeatures(verbose=verbose)%>%
    RunPCA(verbose=verbose) %>%
    FindNeighbors(verbose=verbose) %>%
    FindClusters(verbose=verbose)

  return(
    setNames(Idents(seu), Cells(seu))
  )
}

#     SoupX ----
# https://github.com/constantAmateur/SoupX
NCORES<-min(c(nrow(meta),16))

cl <- makeCluster(NCORES)

#load souplist objects from STARsolo outputs
soup.list <- mclapply(
  paste0(meta$data.dir,"/Solo.out/GeneFull/"),
  FUN = function(DIR){
    tod = Seurat::Read10X(paste0(DIR,'raw')) #droplets
    toc = Seurat::Read10X(paste0(DIR,'filtered')) # cells

    return(
      SoupChannel(
        tod=tod,
        toc=toc
      )
    )
  },
  mc.cores=NCORES
)

# set cluster IDs for cells before soup estimations
soup.list <- mclapply(
  soup.list,
  FUN = function(sc){
    # quick preprocessing/clustering to get cluster IDs for cells
    tmp.clusters <- getClusterIDs(toc=sc$toc)

    return(tryCatch(
      setClusters(sc,tmp.clusters),
      error=function(e) NULL
    ))
  },
  mc.cores=NCORES
)

soup.list.est <- mclapply(
  soup.list,
  FUN = function(sc){
    return(tryCatch(autoEstCont(sc), error=function(e) NULL))
  },
  mc.cores=NCORES
)

adj.mat.list <- mclapply(
  soup.list.est,
  FUN = function(sc){
    return(tryCatch(adjustCounts(sc), error=function(e) NULL))
  },
  mc.cores=NCORES
)
stopCluster(cl)

# Save adjusted matrices to disk
for(i in 1:length(adj.mat.list)){
  if(!is.null(adj.mat.list[[i]]) & !file.exists(paste0(meta$data.dir[i],"/Solo.out/GeneFull/soupx/matrix.mtx.gz"))){
    cat(paste0("Writing matrix for ", meta$sample[i],"...\n"))
    write_sparse(
      path=paste0(meta$data.dir[i],"/Solo.out/GeneFull/soupx"), # name of new directory
      x=adj.mat.list[[i]], # matrix to write as sparse
      barcodes=NULL, # cell IDs, colnames
      features=NULL, # gene IDs, rownames
      overwrite=F,
      verbose=T
    )
    # DropletUtils:::write10xCounts(
    #   path=paste0(meta$data.dir[i],"/soupx"), #path to each sample's STARsolo output
    #   adj.mat.list[[i]]
    # )
  }else{
    message(paste0("Skipping ", meta$sample[i],"...\n"))
  }
}

# save Rho values
rhos <- list()
for(i in 1:length(soup.list.est)){
  rhos[[i]] <- mean(soup.list.est[[i]]$metaData$rho)
}
rhos <- do.call(rbind,rhos)

#TODO- update and save metadata

#     Read in the count matrices ----

#TODO

# mat.list <- list()
# soupx.used <- list()
# for(i in 1:length(meta$data.dir)){
#   # only us SoupX values for samples from whole muscle (i.e., not FACS-sorted)
#   if(file.exists(paste0(meta$data.dir[i], '/soupx')) & meta$is.whole.muscle[i]){
#     cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
#     mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i],"/soupx"))
#     soupx.used[[i]] <- T
#   }else if(file.exists(paste0(meta$data.dir[i], '/outs/filtered_feature_bc_matrix'))){ # use raw counts for FACS-sorted samples
#     cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
#     mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i], '/outs/filtered_feature_bc_matrix'))
#     soupx.used[[i]] <- F
#   }else{
#     cat("Data not found for # ", i, " (", meta$data.dir[i], ")", "\n")
#     soupx.used[[i]] <- NULL
#   }
# }
#
# cat(sum(unlist(lapply(mat.list, ncol))),"cells (total) loaded...\n")
