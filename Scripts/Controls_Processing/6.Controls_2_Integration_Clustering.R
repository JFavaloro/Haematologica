# Script information -------------------------------------------------------------------------------

# Title: Control dataset #2 - Integrate and cluster data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we prepare our SCTransformed data and run a standard dimensionality
# reduction workflow to cluster the integrated data.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '6. Control dataset #2 - Integrate and cluster data' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load filtered, SCTransformed CD8+ T-cells
Controls_2_SCT <- readRDS("Data/R_out/Controls_2/Controls_2_SCT.rds")

# Create a character vector for all samples IDs
all_samples_controls <- apply(expand.grid(c("HD", "MGUS", "MM", "SMM"), 1:11), 
                              1, paste0, collapse=" ") %>% str_replace_all(fixed(" "), "")
discard <- c("HD10", "HD11", 
             "MGUS6", "MGUS7", "MGUS8", "MGUS9", "MGUS10", "MGUS11", 
             "MM8", "MM9", "MM10", "MM11")
all_samples_controls <- setdiff(all_samples_controls, discard)
rm(discard)


# 3. Prepare data for integration  -----------------------------------------------------------------
# NB: As there were insufficient cells to process samples individually, 'lists' of common 'Disease'
# status have been prepared. We will therefore remove the Controls_2 object from Controls_2_SCT and 
# proceed with the remaining workflow.

# Remove the Controls_2 object from the Controls_2_SCT list and rename it all_lists to allow utility
# of code previously written for 10x
all_lists <- Controls_2_SCT[-1]

# Move the Seurat objects one layer down the all_list
all_lists$all.list <- list("HD" = all_lists$HD, 
                           "MGUS" = all_lists$MGUS, 
                           "MM" = all_lists$MM, 
                           "SMM" = all_lists$SMM)

# Remove the unnecessary objects from the all_lists list
all_lists <- all_lists[-c(1:4)]

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


# 5. Integrate data and run dimensionality reduction -----------------------------------------------

# Create empty list to store Integrated data
integrated_list <- list()

# Integrate data and run dimensionality reduction on 10 principal components and find clusters
# NB: Previous analysis showed that bulk of variance is captured in the first 10 PC for all samples
# NB: Seurat advise is to run all dimensionality reduction to same level - i.e. UMAP = 1:10
# NB: For this integrated object resolution of 0.19 will yield 5 clusters, 0.195 will yield 6.
for (i in names(anchors_list)) {
  temp <- IntegrateData(anchors_list[[i]], normalization.method = "SCT") %>%
    RunPCA(dims = 1:50) %>% RunUMAP(dims = 1:10) %>%
    FindNeighbors(reduction = "pca", dims = 1:10) %>% 
    FindClusters(resolution = c(0.19, 0.195))
  integrated_list[[i]] <- temp
}

# Rename items in the integrated_list to .integrated
integrated_names <- names(integrated_list)
integrated_names <- gsub("\\.anchors", ".integrated", integrated_names)
names(integrated_list) <- integrated_names

# Define function to set default cluster resolution to 0.195 - 6 clusters
fix_res <- function(input_seurat){
  input_seurat[["seurat_clusters"]] <- input_seurat[["integrated_snn_res.0.195"]]
  return(input_seurat)
}

# Apply to all Seurat objects
integrated_list <- lapply(integrated_list, fix_res)

# Define function to set Idents to serurat_clusters
fix_ident <- function(input_seurat){
  Idents(input_seurat) <- "seurat_clusters"
  return(input_seurat)
}

# Apply to all Seurat objects
integrated_list <- lapply(integrated_list, fix_ident)

# Check clustering
DimPlot(integrated_list$all)


# 6. Save out data and print stats -----------------------------------------------------------------

saveRDS(anchors_list, "Data/R_out/Controls_2/anchors_list.rds")
saveRDS(integrated_list, "Data/R_out/Controls_2/integrated_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '6. Control dataset #2 - Integrate and cluster data' script")