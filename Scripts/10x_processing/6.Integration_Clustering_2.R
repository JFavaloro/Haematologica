# Script information -------------------------------------------------------------------------------

# Title: Integrate and Cluster data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we repeat our integration and clustering workflow. Although previous 
# analyses demonstrated the bulk of data variance is captured in the first 10 PC, we will calculate 
# the first 50 PC, then run the workflow using 10, 12 and 15 PC to ensure we are not missing 
# information by setting the dimensions too low. We will focus only on three combinations of data 
# moving forward; the standardised combined BM and PB from 2x NDMM patients, BM only from 2x NDMM 
# patients and PB only from 4x NDMM patients. To speed up computation, we will import the 
# previously calculated anchors_list and cull it to the relevant objects.
# NB: Seurat advise is to run all dimensionality reduction to same level - i.e. UMAP = 1:10
# NB: For ease of downstream use and due to the biology of CD8+ T-cells, we will initially target 7 
# clusters across all samples by altering the resolution parameter.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '6. Integration and clustering 2' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire)
})


# 2. Load data and define common functions ---------------------------------------------------------

# Load the SCTransformed, filtered 10x data
all_seurat_filtered_SCT <- readRDS("Data/R_out/10x/all_seurat_filtered_SCT.rds")

# Load the anchors_list and cull to the three objects of interest
anchors_list <- readRDS("Data/R_out/10x/anchors_list.rds")
anchors_list <- anchors_list[-c(1:11)]

# Get all combinations of PB/BM and 13/31/43/63 & create a character vector for downsampled objects
all_samples <- apply(expand.grid(c("BM", "PB"), c("13","31","43","63")), 1, paste, collapse="")
all_samples_ds <- apply(expand.grid(c("BM", "PB"), c("13_ds","31_ds","43_ds","63_ds")), 
                        1, paste, collapse="")

# Define function to set Idents to serurat_clusters
fix_ident <- function(input_seurat){
  Idents(input_seurat) <- "seurat_clusters"
  return(input_seurat)
}


# 3. Prepare data for integration  -----------------------------------------------------------------

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

# Create lists for downsampled data
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

# Create empty list to store Integration features
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


# 4. Integrate data and run dimensionality reduction on 10 PC --------------------------------------

# Create empty list to store Integrated data
integrated_list <- list()

# Integrate data and run dimensionality reduction/clustering on 10 principal components. 
for (i in names(anchors_list)) {
  temp <- IntegrateData(anchors_list[[i]], normalization.method = "SCT") %>%
    RunPCA(dims = 1:50) %>% RunUMAP(dims = 1:10) %>%
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
integrated_list$BM_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$BM_ds.integrated[["integrated_snn_res.0.4"]]
integrated_list$PB_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$PB_ds.integrated[["integrated_snn_res.0.22"]]
integrated_list$all_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$all_ds.integrated[["integrated_snn_res.0.22"]]

# Set Idents to serurat_clusters
integrated_list <- lapply(integrated_list, fix_ident)

# Check clustering
DimPlot(integrated_list$all_ds.integrated)

# Rename the integrated list to allow export
integrated_list_10_PC <- integrated_list 

# 5. Integrate data and run dimensionality reduction on 12 PC --------------------------------------

# Empty the Integrated data list
integrated_list <- list()

# Integrate data and run dimensionality reduction/clustering on 12 principal components. 
for (i in names(anchors_list)) {
  temp <- IntegrateData(anchors_list[[i]], normalization.method = "SCT") %>%
    RunPCA(dims = 1:50) %>% RunUMAP(dims = 1:12) %>%
    FindNeighbors(reduction = "pca", dims = 1:12) %>% 
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
integrated_list$BM_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$BM_ds.integrated[["integrated_snn_res.0.45"]]
integrated_list$PB_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$PB_ds.integrated[["integrated_snn_res.0.2"]]
integrated_list$all_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$all_ds.integrated[["integrated_snn_res.0.24"]]

# Set Idents to serurat_clusters
integrated_list <- lapply(integrated_list, fix_ident)

# Check clustering
DimPlot(integrated_list$all_ds.integrated)

# Rename the integrated_list to allow export
integrated_list_12_PC <- integrated_list


# 6. Integrate data and run dimensionality reduction on 15 PC --------------------------------------

# Empty the Integrated data list
integrated_list <- list()

# Integrate data and run dimensionality reduction/clustering on 15 principal components. 
for (i in names(anchors_list)) {
  temp <- IntegrateData(anchors_list[[i]], normalization.method = "SCT") %>%
    RunPCA(dims = 1:50) %>% RunUMAP(dims = 1:15) %>%
    FindNeighbors(reduction = "pca", dims = 1:15) %>% 
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
integrated_list$BM_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$BM_ds.integrated[["integrated_snn_res.0.35"]]
integrated_list$PB_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$PB_ds.integrated[["integrated_snn_res.0.2"]]
integrated_list$all_ds.integrated[["seurat_clusters"]] <- 
  integrated_list$all_ds.integrated[["integrated_snn_res.0.2"]]

# Set Idents to serurat_clusters
integrated_list <- lapply(integrated_list, fix_ident)

# Check clustering
DimPlot(integrated_list$all_ds.integrated)

# Rename the integrated_list to allow export
integrated_list_15_PC <- integrated_list


# 7. Save out data ---------------------------------------------------------------------------------

saveRDS(integrated_list_10_PC, "Data/R_out/10x/integrated_list_10_PC.rds")
saveRDS(integrated_list_12_PC, "Data/R_out/10x/integrated_list_12_PC.rds")
saveRDS(integrated_list_15_PC, "Data/R_out/10x/integrated_list_15_PC.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '6. Integration and clustering 2' script")