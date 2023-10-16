# Script information -------------------------------------------------------------------------------

# Title: Create 10x reference atlases
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will create reference atlases from our 10x data. With these altases
# we will project public data to compare and contrast our clustering and determine if any genes are
# differentially expressed between health and disease. As our dataset has already been processed 
# appropriately, we will will largely replicate the ProjecTILs vignette, with slight modifications:
# https://carmonalab.github.io/ProjecTILs.demo/build_ref_atlas.html

# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '3. Create 10x reference atlases' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, Matrix, ProjecTILs, openxlsx, umap)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load 10x data
all_list <- readRDS("Data/R_out/10x/all_list.rds")

# Define function to set Idents to serurat_clusters
fix_ident <- function(input_seurat){
  Idents(input_seurat) <- "seurat_clusters"
  return(input_seurat)
}

# Create character vectors of sample names
sample_names <- c("all_10x", "BM_10x", "PB_10x")


# 3. Pre-processing --------------------------------------------------------------------------------

# Set seed to same value used for Seurat clustering - meaning of life I guess?
seed <- 42
set.seed(42)

# Set the number of dimensions to use in PCA and UMAP
ndim_pca <- 10
ndim_umap <- 10

# Set assay to use integrated 
which.assay <- "integrated"

# Set UMAP config to match as closely to Seurat as we can
# NB: We will calculate 50 PCs but use on 10 going forward
# NB: Seurat has switched to using "uwot" - this may be contributing to slight differences 
# observable in clusters - try to switch to Seurat RunUMAP - compare outputs.
umap.config <- umap.defaults
umap.config$n_neighbors <- 30
umap.config$min_dist <- 0.3
umap.config$metric <- "cosine"
umap.config$n_components <- 2
umap.config$random_state <- seed
umap.config$transform_state <- seed
umap.config$umap.method = "uwot"


# 4. Create reference atlases ----------------------------------------------------------------------

# Create empty lists to store calculated PCA/UMAP and modified Seurat objects
all_PCA <- list()
all_UMAP <- list()
all_refs <- list()

#Run a loop
for (i in names(all_list)){
  data.seurat <- all_list[[i]]
  
  # Extract variable features from the integrated object and move into an ordered dataframe
  varfeat <- data.seurat@assays[[which.assay]]@var.features 
  refdata <- data.frame(t(data.seurat@assays[[which.assay]]@data[varfeat,]))
  refdata <- refdata[, sort(colnames(refdata))]
  
  # Calculate new PC and rotate data
  ref.pca <- prcomp(refdata, rank. = 50, scale. = TRUE, center = TRUE, retx = TRUE)
  ref.pca$rotation[1:5,1:5]
  
  # Perform nonlinear dimensionality reduction want this stored in a list for all 3 as above
  ref.umap <- umap(ref.pca$x[,1:ndim_umap], config = umap.config)
  
  # Copy embeddings into the atlas object
  colnames(ref.umap$layout) <- c("UMAP_1","UMAP_2")
  data.seurat@reductions$umap@cell.embeddings <- ref.umap$layout
  data.seurat@reductions$pca@cell.embeddings <- ref.pca$x
  data.seurat@reductions$pca@feature.loadings <- ref.pca$rotation
  colnames(data.seurat@reductions$pca@cell.embeddings) <- 
    gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl=TRUE)
  colnames(data.seurat@reductions$pca@feature.loadings) <- 
    gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl=TRUE)
  #Store the complete PCA and UMAP object in @misc
  data.seurat@misc$pca_object <- ref.pca
  data.seurat@misc$umap_object <- ref.umap
  data.seurat@misc$projecTILs = "custom_atlas"
  
  # Perform FindNeighbors from PCA reduction
  DefaultAssay(data.seurat) <- "integrated"
  data.seurat <- FindNeighbors(data.seurat, reduction = "pca", dims = 1:ndim_pca)
  
  # Perform FindClusters
  data.seurat <- 
    FindClusters(data.seurat, 
                 resolution = c(0.2, 0.25, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.45, 0.5))
  
  # Store the custom reference atlas
  all_refs[[i]] <- data.seurat
  # Store the calculated PCA
  all_PCA[[i]] <- ref.pca
  # Store the calculated UMAP
  all_UMAP[[i]] <- ref.umap
}

# Change the cluster resolution to appropriate values to ensure 7 clusters
all_refs$all_10x[["seurat_clusters"]] <- 
  all_refs$all_10x[["integrated_snn_res.0.42"]]
all_refs$BM_10x[["seurat_clusters"]] <- 
  all_refs$BM_10x[["integrated_snn_res.0.42"]]
all_refs$PB_10x[["seurat_clusters"]] <- 
  all_refs$PB_10x[["integrated_snn_res.0.38"]]

# Fix the identity of the reference atlases back to "seurat_clusters"
all_refs <- lapply(all_refs, fix_ident)


# 5. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_PCA, "Data/R_out/ProjecTILs/all_PCA.rds")
saveRDS(all_UMAP, "Data/R_out/ProjecTILs/all_UMAP.rds")
saveRDS(all_refs, "Data/R_out/ProjecTILs/all_refs.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finshed running '3. Create 10x reference atlases' script")