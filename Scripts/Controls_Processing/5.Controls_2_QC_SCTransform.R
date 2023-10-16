# Script information -------------------------------------------------------------------------------

# Title: Control dataset #2 - QC and SCTransform
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will load the cleaned control dataset from Zavidij et al., 2019 and
# create some control Seurat objects for our project. We will process these samples as per our 10x
# data, however, given the limited size of samples (some HD in single digits!) we will merge all
# samples by disease state.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '5. Control dataset #2 - QC and SCTransform' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, SingleCellExperiment, scater)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load filtered CD8+ T-cells
all_seurat_controls_filtered_CD8 <- 
  readRDS("Data/R_out/Controls_2/all_seurat_controls_filtered_CD8.rds")

# Load cycle genes
load("Data/Public_data/cycle.rda")

# Load Blacklist
black_list <- readRDS("Data/R_out/black_list.rds")

# Create a character vector for all samples IDs
all_samples_controls <- apply(expand.grid(c("HD", "MGUS", "MM", "SMM"), 1:11), 
                              1, paste0, collapse=" ") %>% str_replace_all(fixed(" "), "")
discard <- c("HD10", "HD11", 
             "MGUS6", "MGUS7", "MGUS8", "MGUS9", "MGUS10", "MGUS11", 
             "MM8", "MM9", "MM10", "MM11")
all_samples_controls <- setdiff(all_samples_controls, discard)
rm(discard)


# 3. Perform additional QC and merge samples by disease state --------------------------------------

# 3.1 Define function to add the cell cycle stage --------------------------------------------------
# NB: Have checked through each sample to see if cell cycle needs to be regressed - doesn't appear
# to; several genes not found in our data.
add_cycle <- function(input_seurat){
  input_seurat <- CellCycleScoring(input_seurat, g2m.features = g2m_genes, s.features = s_genes)
  return(input_seurat)
}

# Apply to all filtered Seurat objects
all_seurat_controls_filtered_CD8 <- lapply(all_seurat_controls_filtered_CD8, add_cycle)


# 3.2 Merge samples by disease state ---------------------------------------------------------------
# NB: This could likely be cleaned up a bit but it is functional

# Create some merged Seurat objects
MGUS <- merge(all_seurat_controls_filtered_CD8$`MGUS1`, y = 
                c(all_seurat_controls_filtered_CD8$`MGUS2`, 
                  all_seurat_controls_filtered_CD8$`MGUS3`, 
                  all_seurat_controls_filtered_CD8$`MGUS4`, 
                  all_seurat_controls_filtered_CD8$`MGUS5`))

MM <- merge(all_seurat_controls_filtered_CD8$`MM1`, y = 
              c(all_seurat_controls_filtered_CD8$`MM2`, 
                all_seurat_controls_filtered_CD8$`MM3`, 
                all_seurat_controls_filtered_CD8$`MM4`, 
                all_seurat_controls_filtered_CD8$`MM5`, 
                all_seurat_controls_filtered_CD8$`MM6`, 
                all_seurat_controls_filtered_CD8$`MM7`))

HD <- merge(all_seurat_controls_filtered_CD8$`HD1`, y = 
              c(all_seurat_controls_filtered_CD8$`HD2`, 
                all_seurat_controls_filtered_CD8$`HD3`, 
                all_seurat_controls_filtered_CD8$`HD4`, 
                all_seurat_controls_filtered_CD8$`HD5`, 
                all_seurat_controls_filtered_CD8$`HD6`, 
                all_seurat_controls_filtered_CD8$`HD7`, 
                all_seurat_controls_filtered_CD8$`HD8`, 
                all_seurat_controls_filtered_CD8$`HD9`))

SMM <- merge(all_seurat_controls_filtered_CD8$`SMM1`, y = 
               c(all_seurat_controls_filtered_CD8$`SMM2`, 
                 all_seurat_controls_filtered_CD8$`SMM3`, 
                 all_seurat_controls_filtered_CD8$`SMM4`, 
                 all_seurat_controls_filtered_CD8$`SMM5`, 
                 all_seurat_controls_filtered_CD8$`SMM6`, 
                 all_seurat_controls_filtered_CD8$`SMM7`, 
                 all_seurat_controls_filtered_CD8$`SMM8`, 
                 all_seurat_controls_filtered_CD8$`SMM9`, 
                 all_seurat_controls_filtered_CD8$`SMM10`, 
                 all_seurat_controls_filtered_CD8$`SMM11`))

Controls_2 <- merge(HD, c(MM, MGUS, SMM))

# Move merged Seurat objects into a list
Controls_2 <- list("Controls_2" = Controls_2, "HD" = HD, "MGUS" = MGUS, "MM" = MM, "SMM" = SMM)

# Create a character vector of all merged objects for output of figures
controls <- c("Controls_2", "HD", "MGUS", "MM", "SMM")


# 3.3 Determine unhelpful variance in data ---------------------------------------------------------

# Define function to run NormalizeData
# NB: Need to re-run normalizeData, or merge with the argument "merge.data = TRUE"
norm_data <- function(input_seurat){
  input_seurat <- NormalizeData(input_seurat)
  return(input_seurat)
}

# Apply to filtered Seurat objects
Controls_2 <- lapply(Controls_2, norm_data)

# Create a temp list of Seurat objects to QC
temp_seurat <- Controls_2

# Process data the "old way" using FindVariableFeatures/ScaleData/RunPCA
for (i in 1:length(temp_seurat)) {
  x <- FindVariableFeatures(temp_seurat[[i]], nfeatures = 2000) %>% ScaleData() %>% RunPCA()
  temp_seurat[[i]] <- x
}

# Set names of temp list
names(temp_seurat) <- names(Controls_2)


# 3.4 Determine if cell cycle genes contribute variance to data ------------------------------------

# Plot the PCA colored by cell cycle phase across all samples
for (i in 1:length(temp_seurat)) {
  x <- DimPlot(temp_seurat[[i]], reduction = "pca", group.by= "Phase", split.by = "Phase")
  ggsave(plot = x, filename = paste0("Output/QC/Controls_2/Phase/Controls_2_Phase_", 
                                     controls[i], ".png"))
}


# 3.5 Determine if mitochondrial genes contribute variance to data ---------------------------------

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
  ggsave(plot = x, filename = paste0("Output/QC/Controls_2/Mito/Controls_2_Mito_", 
                                     controls[i], ".png"))
}

# Check the number of useful principal components
for (i in 1:length(temp_seurat)) {
  x <- DimHeatmap(temp_seurat[[i]], reduction = "pca", dims = 1:15, cells = 500, 
                  fast = FALSE, combine = TRUE) + theme(text = element_text(size = 2))
  ggsave(plot = x, filename = paste0("Output/QC/Controls_2/PCA/Controls_2_PCA_", 
                                     controls[i], ".png"), scale = 2:1.8)
}


# 3.6 Determine dimensionality of data -------------------------------------------------------------

# Plot Elbow plots to determine appropriate dimensionality 
for (i in 1:length(temp_seurat)) {
  x <- ElbowPlot(temp_seurat[[i]])
  ggsave(plot = x, filename = 
           paste0("Output/QC/Controls_2/Elbow_plots/Controls_2_Elbow_plots_", 
                  controls[i], ".png"))
}


# 3.7 Determine genes that may add unhelpful variance to data --------------------------------------

# Create empty list to store Integration features
sce_list <- list()

# Convert the temp Seurat object to an sce object and determine which genes may need to be culled
for (i in names(temp_seurat)) {
  temp <- as.SingleCellExperiment(temp_seurat[[i]])
  sce_list[[i]] <- temp
}

# # Plot the 50 highest expressing genes across all samples
for (i in 1:length(sce_list)) {
  x <- plotHighestExprs(sce_list[[i]], 
                        exprs_values = "counts", n = 50, colour_cells_by = "orig.ident")
  ggsave(plot = x, filename = 
           paste0("Output/QC/Controls_2/Expression/Controls_2_Counts_", 
                  controls[i], ".png"))
}

# Looks like it may be necessary to regress percent. mito but not cell cycle


# 4. Run SCTransform and remove genes contributing unhelpful variance from var.features slot -------

# Define function to run SCTransform
sc_trans <- function(input_seurat){
  input_seurat <- SCTransform(input_seurat, vars.to.regress = "percent.mt")
  return(input_seurat)
}

# Apply to filtered Seurat objects
Controls_2_SCT <- lapply(Controls_2, sc_trans)

# Create an empty list to store the original var.features
orig_var.features_list <- list()

# Run loop to backup original var.features
for (i in names(Controls_2_SCT)) {
  temp <- Controls_2_SCT[[i]]@assays[["SCT"]]@var.features
  name <- i
  orig_var.features_list[[name]] <- temp
}

# Define function for removing TCR genes from variable feature slot
remove_genes <- function(input_seurat){
  input_seurat@assays[["SCT"]]@var.features <- 
    input_seurat@assays[["SCT"]]@var.features[!(input_seurat@assays[["SCT"]]@var.features %in% black_list)]
  return(input_seurat)
}

# Apply function on SCTransformed Seurat objects
Controls_2_SCT <- lapply(Controls_2_SCT, remove_genes)


# 5. Save out data and print stats -----------------------------------------------------------------

saveRDS(Controls_2_SCT, "Data/R_out/Controls_2/Controls_2_SCT.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '5. Control dataset #2 - QC and SCTransform' script")