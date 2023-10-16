# Script information -------------------------------------------------------------------------------

# Title: Control dataset #1 - Integrate and cluster data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we prepare our SCTransformed data and run a standard dimensionality
# reduction workflow to cluster the integrated data.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '3. Control dataset #1 - Integrate and cluster data' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load filtered, SCTransformed data
all_seurat_controls_filtered_SCT <- 
  readRDS("Data/R_out/Controls_1/all_seurat_controls_filtered_SCT.rds")

# Create a character vector for all controls
controls <- list.files("Data/Public_data/Controls_1", pattern = "*.gz", full.names = TRUE)
controls <- gsub("Data/Public_data/Controls_1/", "", gsub(".txt.gz", "", controls))


# 3. Prepare data for integration  -----------------------------------------------------------------

# 3.1 Create lists of combinations we're interested in assessing -----------------------------------
# NB: This could likely be cleaned up a bit but it is functional

# Create an empty list of lists
all_lists <- list()

# Combinations are BM, PB, all, rest and act

# Create BM list
all_lists$BM.list <- list(D1_BM_act = all_seurat_controls_filtered_SCT$D1_BM_act, 
                          D1_BM_rest = all_seurat_controls_filtered_SCT$D1_BM_rest, 
                          D2_BM_act = all_seurat_controls_filtered_SCT$D2_BM_act, 
                          D2_BM_rest = all_seurat_controls_filtered_SCT$D2_BM_rest)

# Create PB list
all_lists$PB.list <- list(D1_PB_act = all_seurat_controls_filtered_SCT$D1_PB_act, 
                          D1_PB_rest = all_seurat_controls_filtered_SCT$D1_PB_rest, 
                          D2_PB_act = all_seurat_controls_filtered_SCT$D2_PB_act, 
                          D2_PB_rest = all_seurat_controls_filtered_SCT$D2_PB_rest)

# Create act list
all_lists$act.list <- list(D1_BM_act = all_seurat_controls_filtered_SCT$D1_BM_act, 
                           D2_BM_act = all_seurat_controls_filtered_SCT$D2_BM_act, 
                           D1_PB_act = all_seurat_controls_filtered_SCT$D1_PB_act, 
                           D2_PB_act = all_seurat_controls_filtered_SCT$D2_PB_act)

# Create rest list
all_lists$rest.list <- list(D1_BM_rest = all_seurat_controls_filtered_SCT$D1_BM_rest, 
                            D2_BM_rest = all_seurat_controls_filtered_SCT$D2_BM_rest, 
                            D1_PB_rest = all_seurat_controls_filtered_SCT$D1_PB_rest, 
                            D2_PB_rest = all_seurat_controls_filtered_SCT$D2_PB_rest)

# Create all list
all_lists$all.list <- list(D1_BM_act = all_seurat_controls_filtered_SCT$D1_BM_act, 
                           D1_BM_rest = all_seurat_controls_filtered_SCT$D1_BM_rest, 
                           D1_PB_act = all_seurat_controls_filtered_SCT$D1_PB_act, 
                           D1_PB_rest = all_seurat_controls_filtered_SCT$D1_PB_rest, 
                           D2_BM_act = all_seurat_controls_filtered_SCT$D2_BM_act, 
                           D2_BM_rest = all_seurat_controls_filtered_SCT$D2_BM_rest, 
                           D2_PB_act = all_seurat_controls_filtered_SCT$D2_PB_act, 
                           D2_PB_rest = all_seurat_controls_filtered_SCT$D2_PB_rest) 


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


# 4. Integrate data and run dimensionality reduction -----------------------------------------------

# Create empty list to store Integrated data
integrated_list <- list()

# Integrate data and run dimensionality reduction on 10 principal components and find clusters
# NB: Previous analysis showed that bulk of variance is captured in the first 10 PC for all samples
# NB: Seurat advise is to run all dimensionality reduction to same level - i.e. UMAP = 1:10
for (i in names(anchors_list)) {
  temp = IntegrateData(anchors_list[[i]], normalization.method = "SCT") %>%
    RunPCA(dims = 1:10) %>% RunUMAP(dims = 1:10) %>%
    FindNeighbors(reduction = "pca", dims = 1:10) %>% 
    FindClusters(resolution = c(0.3, 0.4, 0.415, 0.45, 0.5, 0.6, 0.7, 0.8))
  integrated_list[[i]] <- temp
}

# Rename items in the integrated_list to .integrated
integrated_names <- names(integrated_list)
integrated_names <- gsub("\\.anchors", ".integrated", integrated_names)
names(integrated_list) <- integrated_names

# Change the cluster resolution to appropriate values to ensure 7 clusters

# Define function to set default cluster resolution to 0.4
fix_res <- function(input_seurat){
  input_seurat[["seurat_clusters"]] <- input_seurat[["integrated_snn_res.0.4"]]
  return(input_seurat)
}

# Apply to all Seurat objects
integrated_list <- lapply(integrated_list, fix_res)

# Set resolution to 0.45 for All object
integrated_list$all.integrated$seurat_clusters <- 
  integrated_list$all.integrated$integrated_snn_res.0.45
Idents(integrated_list$all.integrated) <- "seurat_clusters"

# Set resolution to 0.7 for BM object
integrated_list$BM.integrated$seurat_clusters <- 
  integrated_list$BM.integrated$integrated_snn_res.0.7
Idents(integrated_list$BM.integrated) <- "seurat_clusters"

# Set resolution to 0.415 for PB object
integrated_list$PB.integrated$seurat_clusters <-
  integrated_list$PB.integrated$integrated_snn_res.0.415
Idents(integrated_list$PB.integrated) <- "seurat_clusters"

# Define function to set Idents to serurat_clusters
fix_ident <- function(input_seurat){
  Idents(input_seurat) <- "seurat_clusters"
  return(input_seurat)
}

# Apply to all Seurat objects
integrated_list <- lapply(integrated_list, fix_ident)

# Check clustering
DimPlot(integrated_list$all.integrated)
DimPlot(integrated_list$BM.integrated)
DimPlot(integrated_list$PB.integrated)


# 5. Save out data and print stats -----------------------------------------------------------------

saveRDS(anchors_list, "Data/R_out/Controls_1/anchors_list.rds")
saveRDS(integrated_list, "Data/R_out/Controls_1/integrated_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '3. Control dataset #1 - Integrate and cluster data' script")