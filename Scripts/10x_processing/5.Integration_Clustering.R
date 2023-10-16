# Script information -------------------------------------------------------------------------------

# Title: Integrate and cluster data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we prepare our SCTransformed data and run a standard dimensionality
# reduction workflow to cluster the integrated data. We will first run the script targeting the
# first 10 principal components (PC). We will calculate all combinations of our data at 10 PC
# NB: Seurat advise is to run all dimensionality reduction to same level - i.e. UMAP = 1:10
# NB: For ease of downstream use and due to the biology of CD8+ T-cells, we will initially target 7 
# clusters across all samples by altering the resolution parameter.

# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '5. Integration and clustering' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire)
})


# 2. Load data and define common functions ---------------------------------------------------------

# Load the SCTransformed, filtered 10x data
all_seurat_filtered_SCT <- readRDS("Data/R_out/10x/all_seurat_filtered_SCT.rds")

# Get all combinations of PB/BM and 13/31/43/63 & create a character vector for downsampled objects
all_samples <- apply(expand.grid(c("BM", "PB"), c("13","31","43","63")), 1, paste, collapse="")
all_samples_ds <- apply(expand.grid(c("BM", "PB"), c("13_ds","31_ds","43_ds","63_ds")), 
                        1, paste, collapse="")

# Define function to set Idents to serurat_clusters
fix_ident <- function(input_seurat){
  Idents(input_seurat) <- "seurat_clusters"
  return(input_seurat)
}

# Define function to fix clonal bin order
fix_clones <- function(input_seurat){
  slot(input_seurat, "meta.data")$cloneType =
    factor(slot(input_seurat, "meta.data")$cloneType, 
           levels = c("Expanded (0.1 < X <= 1)", 
                      "Large (0.01 < X <= 0.1)", 
                      "Medium (0.001 < X <= 0.01)", 
                      "Small (0 < X <= 0.001)", NA))
  return(input_seurat)
}

# Define function to define clonal bins
# NB: This is no longer required as clonal bins are now defined based on proportion
#fix_clones <- function(input_seurat){
#  slot(input_seurat, "meta.data")$cloneType <- factor(slot(input_seurat, "meta.data")$cloneType, 
#                                           levels = c("Expanded (628 < X <= 6275)", 
#                                                      "Large (63 < X <= 628)", 
#                                                      "Medium (6 < X <= 63)", 
#                                                      "Small (1 < X <= 6)", 
#                                                      "Single (0 < X <= 1)", NA))
#  return(input_seurat)
#}


# 3. Prepare data for integration  -----------------------------------------------------------------

# NB: 2021-02-15: Decision to exclude BM13 and BM31 due to quality concerns.
# NB: To allow clonal binning we will reduce the number of cells for analysis to the lowest common
# denominator (PB63: 6275 cells). This is necessary to allow proportions on a single cell level, not
# a proportion (i.e. you can't have half a cell). Using head will remove the randomness associated 
# with random downsampling with sample, allowing reproducibility but will reduce the number of cells
# available to analyse by ~10%. Alternatively - could use sample with a seed...
# Addendum 2022-02-22: This is no longer required as scRepertorire now allows clonal bining based
# on frequency. We will still run the code, as binning by frequency will not allow identification of
# single clonotype and will bin them with "small" expansions. Options are good...

# Create empty list for downsampled data
all_seurat_filtered_SCT_ds <- list()

# Downsample all cells, excluding BM31 (3rd item of list) (only 2429 cells recovered)
for (i in 1:length(all_samples)) {
  if (i == c(3)) next
  temp <- subset(all_seurat_filtered_SCT[[i]], 
                 cells = head(Cells(all_seurat_filtered_SCT[[i]]), 6275))
  all_seurat_filtered_SCT_ds[[i]] <- temp
}

# Replace the 3rd item in the list with the actual sample
all_seurat_filtered_SCT_ds[3] <- all_seurat_filtered_SCT[3]

# Rename the items in the list
names(all_seurat_filtered_SCT_ds) <- all_samples_ds

# Create an empty list of lists
all_lists <- list()

# Create lists for integration of both raw and downsampled data
# NB: This could likely be cleaned up a bit but it is functional
for (patient in c(13, 31, 43, 63)){
  all_lists[[paste0("pt",patient,".list")]] <- list(
    all_seurat_filtered_SCT[[paste0("BM", patient)]],
    all_seurat_filtered_SCT[[paste0("PB", patient)]])
  names(all_lists[[paste0("pt",patient,".list")]]) <- 
    c(paste0("BM", patient), paste0("PB", patient))
}

all_lists$BM.list <- 
  list(BM43 = all_seurat_filtered_SCT$BM43, BM63 = all_seurat_filtered_SCT$BM63)
all_lists$PB.list <- 
  list(PB13 = all_seurat_filtered_SCT$PB13, PB31 = all_seurat_filtered_SCT$PB31, 
       PB43 = all_seurat_filtered_SCT$PB43, PB63 = all_seurat_filtered_SCT$PB63)
all_lists$all.list <- 
  list(BM43 = all_seurat_filtered_SCT$BM43, PB43 = all_seurat_filtered_SCT$PB43, 
       BM63 = all_seurat_filtered_SCT$BM63, PB63 = all_seurat_filtered_SCT$PB63)

for (patient in c(13, 31, 43, 63)){
  all_lists[[paste0("pt",patient,"_ds.list")]] <- list(
    all_seurat_filtered_SCT_ds[[paste0("BM", patient,"_ds")]],
    all_seurat_filtered_SCT_ds[[paste0("PB", patient,"_ds")]])
  names(all_lists[[paste0("pt",patient,"_ds.list")]]) <- 
    c(paste0("BM", patient), paste0("PB", patient))
}

all_lists$BM_ds.list <- 
  list(BM43_ds = all_seurat_filtered_SCT_ds$BM43_ds, 
       BM63_ds = all_seurat_filtered_SCT_ds$BM63_ds)
all_lists$PB_ds.list <- 
  list(PB13_ds = all_seurat_filtered_SCT_ds$PB13_ds, 
       PB31_ds = all_seurat_filtered_SCT_ds$PB31_ds, 
       PB43_ds = all_seurat_filtered_SCT_ds$PB43_ds, 
       PB63_ds = all_seurat_filtered_SCT_ds$PB63_ds)
all_lists$all_ds.list <- 
  list(BM43_ds = all_seurat_filtered_SCT_ds$BM43_ds, 
       PB43_ds = all_seurat_filtered_SCT_ds$PB43_ds, 
       BM63_ds = all_seurat_filtered_SCT_ds$BM63_ds, 
       PB63_ds = all_seurat_filtered_SCT_ds$PB63_ds)

# Create emptylist to store Integration features
features_list <- list()

# Run loop to select Integration features
for (i in names(all_lists)) {
  temp <- SelectIntegrationFeatures(all_lists[[i]], nfeatures = 3000)
  name <- gsub("\\.list", ".features", i)
  features_list[[name]] <- temp
}

# Run Prep SCT Integration
all_lists <- map2(all_lists, features_list,
                  function(x, y){
                    PrepSCTIntegration(
                      object.list = x, 
                      anchor.features = y
                    )
                  }
)

# Find Integration Anchors
anchors_list <- map2(all_lists, features_list,
                     function(x, y){
                       FindIntegrationAnchors(
                         object.list = x,
                         normalization.method = "SCT",
                         anchor.features = y
                       )
                     }
)

# Rename items in the anchors list to .anchors
anchor_names <- names(anchors_list)
anchor_names <- gsub("\\.list", ".anchors", anchor_names)
names(anchors_list) <- anchor_names


# 5. Integrate data and run dimensionality reduction -----------------------------------------------

# Create empty list to store Integrated data
integrated_list <- list()

# Integrate data and run dimensionality reduction on 10 principal components and find clusters
# NB: Previous analysis showed that bulk of variance is captured in the first 10 PC for all samples
# NB: Seurat advise is to run all dimensionality reduction to same level - i.e. UMAP = 1:10
# NB: For ease of downstream use and due to the biology of CD8+ T-cells, we will target 7 clusters
# across all samples by altering the resolution
for (i in names(anchors_list)) {
  temp <- IntegrateData(anchors_list[[i]], normalization.method = "SCT") %>%
    RunPCA(dims = 1:10) %>% RunUMAP(dims = 1:10) %>%
    FindNeighbors(reduction = "pca", dims = 1:10) %>% 
    FindClusters(resolution = c(0.12, 0.14, 0.145, 0.15, 0.155, 0.16, 0.178, 0.18, 0.19, 0.2, 0.202, 
                                0.21, 0.22, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 
                                1.0))
  integrated_list[[i]] <- temp
}

# Rename items in the integrated_list to .integrated
integrated_names <- names(integrated_list)
integrated_names <- gsub("\\.anchors", ".integrated", integrated_names)
names(integrated_list) <- integrated_names

# Change the cluster resolution to appropriate values to ensure 7 clusters
integrated_list$pt13.integrated[["seurat_clusters"]] <- 
  integrated_list$pt13.integrated[["integrated_snn_res.0.22"]]
integrated_list$pt31.integrated[["seurat_clusters"]] <- 
  integrated_list$pt31.integrated[["integrated_snn_res.0.26"]]
integrated_list$pt43.integrated[["seurat_clusters"]] <- 
  integrated_list$pt43.integrated[["integrated_snn_res.0.5"]]
integrated_list$pt63.integrated[["seurat_clusters"]] <- 
  integrated_list$pt63.integrated[["integrated_snn_res.0.26"]]
integrated_list$BM.integrated[["seurat_clusters"]] <- 
  integrated_list$BM.integrated[["integrated_snn_res.0.35"]]
integrated_list$PB.integrated[["seurat_clusters"]] <- 
  integrated_list$PB.integrated[["integrated_snn_res.0.26"]]
integrated_list$all.integrated[["seurat_clusters"]] <- 
  integrated_list$all.integrated[["integrated_snn_res.0.22"]]
integrated_list$pt13_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$pt13_ds.integrated[["integrated_snn_res.0.2"]]
integrated_list$pt31_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$pt31_ds.integrated[["integrated_snn_res.0.28"]]
integrated_list$pt43_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$pt43_ds.integrated[["integrated_snn_res.0.45"]]
integrated_list$pt63_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$pt63_ds.integrated[["integrated_snn_res.0.3"]]
integrated_list$BM_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$BM_ds.integrated[["integrated_snn_res.0.4"]]
integrated_list$PB_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$PB_ds.integrated[["integrated_snn_res.0.22"]]
integrated_list$all_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$all_ds.integrated[["integrated_snn_res.0.22"]]

# Set Idents to serurat_clusters
integrated_list <- lapply(integrated_list, fix_ident)

# Fix clonal bin order
integrated_list <- lapply(integrated_list, fix_clones)

# Check clustering
DimPlot(integrated_list$all_ds.integrated)


# 6. Save out data ---------------------------------------------------------------------------------

saveRDS(anchors_list, "Data/R_out/10x/anchors_list.rds")
saveRDS(integrated_list, "Data/R_out/10x/integrated_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '5. Integration and clustering' script")