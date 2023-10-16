# Script information -------------------------------------------------------------------------------

# Title: Control dataset #1 - QC and SCTransform
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will load the cleaned control dataset from Szabo et al., 2019 and
# associated metadata and create approximate age-matched control Seurat objects for our project.
# We will process these samples as per our 10x data, and create a number of Seurat objects with
# various combinations of tissue and stimulation status.
# NB: We will process these samples in a near identical way to our data, however, as the file 
# formats are different, slight modifications will be required.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '2. Control dataset #1 - QC and SCTransform' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, cowplot, SingleCellExperiment, scater)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Import GEX data and metadata
raw_data_clean <- readRDS("Data/R_out/Controls_1/raw_data_clean.rds")
metadata_list <- readRDS("Data/R_out/Controls_1/metadata_list.rds")

# Load cell cycle genes
load("Data/Public_data/cycle.rda")

# Load blacklist
black_list <- readRDS("Data/R_out/black_list.rds")

# Define a function to pull metadata columns from our list of Seurat objects
get_metadata_column <- function(input_seurat_list, column_name){
  # pull all of the column from the Seurat object
  all_cols <- lapply(input_seurat_list, function(x) x[[]][,column_name])
  return(data.frame(lapply(all_cols, "length<-", max(lengths(all_cols)))) %>% gather) 
  # convert list to data.frame, then from wide to long format
}

# Create a character vector for all controls
controls <- list.files("Data/Public_data/Controls_1", pattern = "*.gz", full.names = TRUE)
controls <- gsub("Data/Public_data/Controls_1/", "", gsub(".txt.gz", "", controls))


# 3. Create Control Seurat objects -----------------------------------------------------------------

# 3.1 Create Seurat objects ------------------------------------------------------------------------

# Create an empty list to store our Seurat objects
all_seurat_controls <- list()

# Create the Seurat objects using the same minimums as our 10x data
for (i in names(raw_data_clean)) {
  x <- as.matrix(raw_data_clean[[i]][,2:ncol(raw_data_clean[[i]])])
  row.names(x) <- raw_data_clean[[i]][["Gene"]]
  all_seurat_controls[[i]] <- 
    CreateSeuratObject(x, project = "i", min.cells = 3, min.features = 200)
}


# 3.2 Use metadata to remove non CD8+ T-cells from Control dataset and add metadata ----------------

# Remove non CD8+ T-cells from Controls
for (i in names(all_seurat_controls)) {
  x <- subset(all_seurat_controls[[i]], cells = metadata_list[[i]][["barcode"]])
  all_seurat_controls[i] <- x           
}

# Set the rownames of the metadata_list to the cells barcode
for (i in names(metadata_list)) {
  rownames(metadata_list[[i]]) <-  metadata_list[[i]][["barcode"]]      
}

# Add remaining metadata
for (i in names(all_seurat_controls)) {
  all_seurat_controls[[i]] <- AddMetaData(all_seurat_controls[[i]], metadata_list[[i]])
}

# Add percent mitochondrial genes for QC
for (i in names(all_seurat_controls)) {
  all_seurat_controls[[i]][["percent.mt"]] <-  
    PercentageFeatureSet(all_seurat_controls[[i]], pattern = "^MT-")
}


# 3.3 Apply some filtering to remove poor quality cells --------------------------------------------

# Get some QC metrics for plotting
nCount_RNA <- get_metadata_column(input_seurat_list <- 
                                    all_seurat_controls, column_name = "nCount_RNA")
nFeature_RNA <- get_metadata_column(input_seurat_list <- 
                                      all_seurat_controls, column_name = "nFeature_RNA")
percent.mt <- get_metadata_column(input_seurat_list <- 
                                    all_seurat_controls, column_name = "percent.mt")

# Define function to produce summary statistics (median and +/- 1, 1.5 & 2*sd) for plotting
for (i in c(1,1.5,2)){
  data_summary <- function(x) {
    m <- median(x)
    ymin <- m-(i*sd(x))
    ymax <- m+i*(sd(x))
    return(c(y=m,ymin=ymin,ymax=ymax))}
  
  p1 <- ggplot(nCount_RNA, aes(key, value)) + geom_violin() + ylim(NA, 1E4) + 
    ggtitle("nCount_RNA") + xlab("sample") + ylab("UMIs detected") + 
    stat_summary(fun.data=data_summary, color = "red", geom = "pointrange", size =1)
  p2 <- ggplot(nFeature_RNA, aes(key, value)) + geom_violin() + ylim(NA, 3E3) + 
    ggtitle("nFeature_RNA") + xlab("sample") + ylab("Genes detected") + 
    stat_summary(fun.data=data_summary, color = "red", geom = "pointrange", size =1)
  p3 <- ggplot(percent.mt, aes(key, value)) + geom_violin() + ylim(NA, 20) + 
    ggtitle("percent.mt") + xlab("sample") + ylab("Mitochondrial gene %") + 
    stat_summary(fun.data=data_summary, color = "red", geom = "pointrange", size =1)
  ggsave(plot = plot_grid(p1,p2,p3, ncol = 3), 
         filename = paste0("Output/QC/Controls_1/Filter/QC_metrics_overview_", 
                           i,"SD", ".png"), width = 30, height = 10)
}


# 4. Filter data -----------------------------------------------------------------------------------

# NB: Analysis suggests +/- 2 SD for all QC metrics looks appropriate
# Define function to filter Seurat objects based on the median +/- 2 SD for all three QC metrics
filter_seurat <- function(input_seurat){
  filtered_seurat <- subset(input_seurat, subset =
                              nFeature_RNA > median(input_seurat$nFeature_RNA) - 
                              (2*sd(input_seurat$nFeature_RNA)) &
                              nFeature_RNA < median(input_seurat$nFeature_RNA) + 
                              (2*sd(input_seurat$nFeature_RNA)) &
                              nCount_RNA > median(input_seurat$nCount_RNA) - 
                              (2*sd(input_seurat$nCount_RNA)) &
                              nCount_RNA < median(input_seurat$nCount_RNA) + 
                              (2*sd(input_seurat$nCount_RNA)) &
                              percent.mt > median(input_seurat$percent.mt) - 
                              (2*sd(input_seurat$percent.mt)) &
                              percent.mt < median(input_seurat$percent.mt) + 
                              (2*sd(input_seurat$percent.mt)))
  return(filtered_seurat)
}

# Apply filter
all_seurat_controls_filtered <- lapply(all_seurat_controls, filter_seurat)


# 5. Normalise RNA data and add cell cycle metrics -------------------------------------------------

# Define function to run NormalizeData
# NB: This is fine to do here before SCT - have compared pre and post with same results
norm_data <- function(input_seurat){
  input_seurat <- NormalizeData(input_seurat)
  return(input_seurat)
}

# Apply to filtered Seurat objects
all_seurat_controls_filtered <- lapply(all_seurat_controls_filtered, norm_data)

# Define function to add the cell cycle stage
# NB: Have checked through each sample to see if cell cycle needs to be regressed - doesn't appear
# to; several genes not found in our data.
add_cycle <- function(input_seurat){
  input_seurat <- CellCycleScoring(input_seurat, g2m.features = g2m_genes, s.features = s_genes)
  return(input_seurat)
}

# Apply to all filtered Seurat objects
all_seurat_controls_filtered <- lapply(all_seurat_controls_filtered, add_cycle)


# 6. Determine unhelpful variance in data ----------------------------------------------------------

# Create a temp list of Seurat objects to QC
temp_seurat <- all_seurat_controls_filtered

# Process data the "old way" using FindVariableFeatures/ScaleData/RunPCA
for (i in 1:length(temp_seurat)) {
  x <- FindVariableFeatures(temp_seurat[[i]], nfeatures = 2000) %>% ScaleData() %>% RunPCA()
  temp_seurat[[i]] <- x
}

# Set names of temp list
names(temp_seurat) <- names(all_seurat_controls_filtered)


# 6.1 Determine if cell cycle genes contribute variance to data ------------------------------------

# Plot the PCA colored by cell cycle phase across all samples
for (i in 1:length(temp_seurat)) {
  x <- DimPlot(temp_seurat[[i]], reduction = "pca", group.by= "Phase", split.by = "Phase")
  ggsave(plot = x, filename = paste0("Output/QC/Controls_1/Phase/Controls_1_Phase_", 
                                     controls[i], ".png"))
}


# 6.2 Determine if mitochondrial genes contribute variance to data ---------------------------------

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
  ggsave(plot = x, filename = paste0("Output/QC/Controls_1/Mito/Controls_1_Mito_", 
                                     controls[i], ".png"))
}

# Check the number of useful principal components
for (i in 1:length(temp_seurat)) {
  x <- DimHeatmap(temp_seurat[[i]], reduction = "pca", dims = 1:15, cells = 500, 
                  fast = FALSE, combine = TRUE)
  ggsave(plot = x, filename = paste0("Output/QC/Controls_1/PCA/Controls_1_PCA_", 
                                     controls[i], ".png"), scale = 2:1.8)
}


# 6.3 Determine dimensionality of data -------------------------------------------------------------

# Plot Elbow plots to determine appropriate dimensionality 
for (i in 1:length(temp_seurat)) {
  x <- ElbowPlot(temp_seurat[[i]])
  ggsave(plot = x, filename = 
           paste0("Output/QC/Controls_1/Elbow_plots/Controls_1_Elbow_plots_", 
                  controls[i], ".png"))
}


# 6.4 Determine genes that may add unhelpful variance to data --------------------------------------

# Create empty list to store Integration features
sce_list <- list()

# Convert the temp Seurat object to an sce object and determine which genes may need to be culled
for (i in names(temp_seurat)) {
  temp <- as.SingleCellExperiment(temp_seurat[[i]])
  sce_list[[i]] <- temp
}

# Plot the 50 highest expressing genes across all samples
for (i in 1:length(sce_list)) {
  x <- plotHighestExprs(sce_list[[i]], exprs_values = "counts", n = 50, 
                        colour_cells_by = "orig.ident")
  ggsave(plot = x, filename = 
           paste0("Output/QC/Controls_1/Expression/Controls_1_Counts_", 
                  controls[i], ".png"))
}

# Looks like it may be necessary to regress percent. mito but not cell cycle


# 5. Run SCTransform and remove genes contributing unhelpful variance from var.features slot -------

# Define function to run SCTransform
sc_trans <- function(input_seurat){
  input_seurat <- SCTransform(input_seurat, vars.to.regress = "percent.mt")
  return(input_seurat)
}

# Apply to filtered Seurat objects
all_seurat_controls_filtered_SCT <- lapply(all_seurat_controls_filtered, sc_trans)

# Create an empty list to store the original var.features
orig_var.features_list <- list()

# Run loop to backup original var.features
for (i in names(all_seurat_controls_filtered_SCT)) {
  temp <- all_seurat_controls_filtered_SCT[[i]]@assays[["SCT"]]@var.features
  name <- i
  orig_var.features_list[[name]] <- temp
}

# Define function for removing TCR genes from variable feature slot
remove_genes <- function(input_seurat){
  input_seurat@assays[["SCT"]]@var.features <- 
    input_seurat@assays[["SCT"]]@var.features[!(input_seurat@assays[["SCT"]]@var.features %in% black_list)]
  return(input_seurat)
}

# Apply function on SCTransformed seurat object
all_seurat_controls_filtered_SCT <- lapply(all_seurat_controls_filtered_SCT, remove_genes)


# 6. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_seurat_controls, "Data/R_out/Controls_1/all_seurat_controls.rds")
saveRDS(all_seurat_controls_filtered, "Data/R_out/Controls_1/all_seurat_controls_filtered.rds")
saveRDS(all_seurat_controls_filtered_SCT, "Data/R_out/Controls_1/all_seurat_controls_filtered_SCT.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '2. Control dataset #1 - QC and SCTransform' script")