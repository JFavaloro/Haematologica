# Script information -------------------------------------------------------------------------------

# Title: Further QC
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will perform some further QC and determine genes that may
# contribute to unhelpful variance using the vignettes:
# https://github.com/hbctraining/scRNA-seq/blob/master/lessons/06_SC_SCT_and_integration.md
# https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html
# We will then normalise the data using SCTransform for downstream processing.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '4. Further QC and SCTransform' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire, SingleCellExperiment, scater)
})


# 2. Load data, define functions and create character vectors --------------------------------------

all_seurat_filtered <- readRDS("Data/R_out/10x/all_seurat_filtered.rds")
all_gex <- readRDS("Data/R_out/10x/all_gex.rds")

# Load cell cycle genes
load("Data/Public_data/cycle.rda")

# Get all combinations of PB/BM and 13/31/43/63 & create a character vector
all_samples <- apply(expand.grid(c("BM", "PB"), c("13","31","43","63")), 1, paste, collapse="")


# 3. Normalise RNA data and add cell cycle metrics -------------------------------------------------

# Define function to run NormalizeData
# NB: This is fine to do here before SCT - have compared pre and post with same results
norm_data <- function(input_seurat){
  input_seurat <- NormalizeData(input_seurat)
  return(input_seurat)
}

# Apply to filtered Seurat objects
all_seurat_filtered <- lapply(all_seurat_filtered, norm_data)

# Define function to add the cell cycle stage
# NB: Have checked through each sample to see if cell cycle needs to be regressed - doesn't appear
# to; several genes not found in our data.
add_cycle <- function(input_seurat){
  input_seurat <- CellCycleScoring(input_seurat, g2m.features = g2m_genes, s.features = s_genes)
  return(input_seurat)
}

# Apply to all SCTransformed Seurat objects
all_seurat_filtered <- lapply(all_seurat_filtered, add_cycle)


# 4. Determine unhelpful variance in data ----------------------------------------------------------

# Create a temp list of Seurat objects to QC
temp_seurat <- all_seurat_filtered

# Process data the "old way" using FindVariableFeatures/ScaleData/RunPCA
for (i in 1:length(temp_seurat)) {
  x <- FindVariableFeatures(temp_seurat[[i]], nfeatures = 2000) %>% ScaleData() %>% RunPCA()
  temp_seurat[[i]] <- x
}

# Set names of temp list
names(temp_seurat) <- names(all_seurat_filtered)

# 4.1 Determine if cell cycle genes contribute variance to data ------------------------------------

# Plot the PCA colored by cell cycle phase across all samples
for (i in 1:length(temp_seurat)) {
  x <- DimPlot(temp_seurat[[i]], reduction = "pca", group.by= "Phase", split.by = "Phase")
  ggsave(plot = x, filename = paste0("Output/QC/10x/Phase/Phase_", all_samples[i], ".png"))
}


# 4.2 Determine if mitochondrial genes contribute variance to data ---------------------------------

# Create an empty list to hold break points for plotting mitochondrial data
mito_data <- list()

# Determine the breaks used for plotting data
for (i in 1:length(temp_seurat)) {
  x <- summary(temp_seurat[[i]]@meta.data[["percent.mt"]])
  mito_data[[i]] <- x
}

# Rename objects in the list
names(mito_data) <- names(temp_seurat)

# Add mitochondrial fraction data to temp_seurat objects
for (i in 1:length(temp_seurat)) {
  x <- cut(temp_seurat[[i]]@meta.data[["percent.mt"]],
           breaks=mito_data[[i]], 
           labels=c("Low","Low-mid","Mid","Mid-high","High"))
  temp_seurat[[i]]$mitoFr <- x
}

# Plot the PCA colored by cell cycle phase across all samples
for (i in 1:length(temp_seurat)) {
  x <- DimPlot(temp_seurat[[i]], reduction = "pca", group.by= "mitoFr", split.by = "mitoFr")
  ggsave(plot = x, filename = paste0("Output/QC/10x/Mito/Mito_", all_samples[i], ".png"))
}

# Check the number of useful principal components
for (i in 1:length(temp_seurat)) {
  x <- DimHeatmap(temp_seurat[[i]], reduction = "pca", dims = 1:15, cells = 500,
                  fast = FALSE, combine = TRUE)
  ggsave(plot = x, 
         filename = paste0("Output/QC/10x/PCA/PCA_", all_samples[i], ".png"), scale = 2:1.8)
}

# Plot Elbow plots to determine appropriate dimensionality 
for (i in 1:length(temp_seurat)) {
  x <- ElbowPlot(temp_seurat[[i]])
  ggsave(plot = x, 
         filename = paste0("Output/QC/10x/Elbow_plots/Elbow_plots_", 
                           all_samples[i], ".png"))
}


# 4.3 Determine genes that may add unhelpful variance to data --------------------------------------

# Create empty list to store Integration features
sce_list <- list()

# Convert the temp Seurat object to an sce object and determine which genes may need to be culled
for (i in names(temp_seurat)) {
  temp <- as.SingleCellExperiment(temp_seurat[[i]])
  sce_list[[i]] <- temp
}

# # Plot the 50 highest expressing genes across all samples
for (i in 1:length(sce_list)) {
  x <- plotHighestExprs(sce_list[[i]], exprs_values = "counts", n = 50, colour_cells_by = "Patient")
  ggsave(plot = x, 
         filename = paste0("Output/QC/10x/Expression/Counts_", all_samples[i], ".png"))
}

# Looks like it may be necessary to regress percent. mito but not cell cycle
# Need to remove TCR genes due to biology (cells should express only 1 alpha/beta pair)
# Remove mitochondrial, ribosomal and histone genes as well as some other highly expressed genes


# 5. Run SCTransform and remove genes contributing unhelpful variance from var.features slot -------

# Define function to run SCTransform
sc_trans <- function(input_seurat){
  input_seurat <- SCTransform(input_seurat, vars.to.regress = "percent.mt")
  return(input_seurat)
}

# Apply to filtered Seurat objects
all_seurat_filtered_SCT <- lapply(all_seurat_filtered, sc_trans)

# Create an empty list to store the original var.features
orig_var.features_list <- list()

# Run loop to backup original var.features
for (i in names(all_seurat_filtered_SCT)) {
  temp <- all_seurat_filtered_SCT[[i]]@assays[["SCT"]]@var.features
  name <- i
  orig_var.features_list[[name]] <- temp
}

# Create list of genes to exclude from clustering
# NB: Other genes that may be worthwhile excluding: |ACTB|B2M|EEF1A1|TMSB4X|TMSB10|
# Create a temporary Seurat object containing all 33358 genes in our data set and exclude genes that
# influence clustering in a non helpful manner
temp_seurat <- CreateSeuratObject(all_gex$PB43)
black_list <- 
  grep(pattern = "^MT-|^MTRNR|^HIST|RPS|RPL|TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC|MALAT1", 
       x <- rownames(x = temp_seurat@assays$RNA@counts), value = TRUE)

# Define function for removing TCR genes from variable feature slot
remove_genes <- function(input_seurat){
  input_seurat@assays[["SCT"]]@var.features <- 
    input_seurat@assays[["SCT"]]@var.features[!(input_seurat@assays[["SCT"]]@var.features %in% black_list)]
  return(input_seurat)
}

# Apply function on SCTransformed Seurat object
all_seurat_filtered_SCT <- lapply(all_seurat_filtered_SCT, remove_genes)


# 6. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_seurat_filtered_SCT, "Data/R_out/10x/all_seurat_filtered_SCT.rds")
saveRDS(black_list, "Data/R_out/black_list.rds")
saveRDS(orig_var.features_list, "Data/R_out/10x/orig_var.features_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '4. Further QC and SCTransform' script")