# Script information -------------------------------------------------------------------------------

# Title: Control dataset #3 - Complete workflow
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will process data from Mogilenko et al., 2021: 
# https://doi.org/10.1016/j.immuni.2020.11.005. This dataset comprises CD8+ T-cells sourced from 
# the PB of 11 young and 10 ald (age-matched) controls. We will replicate the authors original
# processing (modifying slightly to allow segregation of young and old) then perform an integrated
# analysis to allow downstream use with ProjecTILs to allow comparison with our 10x dataset.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '7. Control dataset #3 - Complete workflow' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, Matrix, readr, openxlsx)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Read in datamatrix from https://www.cell.com/immunity/fulltext/S1074-7613(20)30492-1
load("Data/Public_data/Controls_3/expData.Rda")


# 3. Create Seurat object --------------------------------------------------------------------------

whole <- CreateSeuratObject(expData)

# Add metadata separating old and young donors
# NB: Uncertain why some slots return NA
young <- c("A06", "A08", "A10", "A11", "A12", "A18", "A20", "A21", "A23", "A24", "A26")
old <- c("D03", "D15", "D22", "D24", "E04", "E05", "E08", "E10", "E16", "E17")
whole@meta.data$orig.ident <- as.character(whole@meta.data$orig.ident)
whole@meta.data <- whole@meta.data %>%
  mutate(age = case_when(
    startsWith(orig.ident, "A") ~ "young",
    startsWith(orig.ident, c("D", "E")) ~ "old"
  ))
whole@meta.data$age <- whole@meta.data$age %>% replace_na('old') # Why is this necessary?
whole@meta.data$orig.ident <- as.factor(whole@meta.data$orig.ident)


# 4. Run the authors original script to process their data -----------------------------------------
# NB: Script has been modified to remove pre-processing of the raw data as it is already provided
# removing of saving intermediate outputs and altering cluster names from GZMB-EM to TTE

## NORMALIZATION
mito.genes <- c(grep("^MT-", rownames(x = whole), value = T),
                grep("^mt-", rownames(x = whole), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts')[mito.genes, ]) / 
  Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts'))
whole[['percent.mito']] <- percent.mito

whole <- subset(x = whole, subset = percent.mito <= 0.05)
whole <- NormalizeData(object = whole, normalization.method = "LogNormalize", scale.factor = 10000)
whole <- FindVariableFeatures(object = whole, selection.method = 'mean.var.plot', 
                              mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
whole@assays$RNA@var.features <-  whole@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^CD138-Antibody|^CD20-Antibody|^IgM-Antibody|^CD5-Antibody|^CD24-Antibody|^IgG-Antibody|^IgD-Antibody|^CD38-Antibody|^CD3-Antibody|^CD4-Antibody|^CD8-Antibody|^CD56-Antibody|^CD45RA-Antibody|^CD194-Antibody|^CD25-Antibody|^CD45RO-Antibody|^CD279-Antibody|^TIGIT-Antibody|^CD294-Antibody|^CD183-Antibody|^CD195-Antibody|^CD196-Antibody|^CD185-Antibody|^CD103-Antibody|^CD69-Antibody|^CD62L-Antibody|^CD197-Antibody|^CD152-Antibody|^CD223-Antibody|^KLRG1-Antibody|^CD27-Antibody|^CD314-Antibody|^CD57-Antibody|^CD366-Antibody|^CX3CR1-Antibody|^TCR-Antibody|^CD357-Antibody|^CD184-Antibody|^CD28-Antibody|^CD127-Antibody", whole@assays$RNA@var.features)]
whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), 
                   vars.to.regress = c("nCount_RNA", "percent.mito"))
gc()

## PCA
whole <- RunPCA(object = whole,
                features =  VariableFeatures(object = whole),
                dims = 1:30,
                verbose=FALSE)


## TSNE
whole <- RunTSNE(object = whole, dims = 1:30)

whole <- RunUMAP(object = whole, dims = 1:30)

## CLUSTERING
whole <- FindNeighbors(object = whole, dims = 1:30)
whole <- FindClusters(object = whole, resolution = 0.3, genes.use = whole@assays$RNA@var.features)

dataForPlot <- as.data.frame(whole@reductions$tsne@cell.embeddings)
dataForPlot$Sample <- whole@meta.data$orig.ident
dataForPlot$Cluster <-  Idents(object = whole)
dataForPlot$nUmi <- whole@meta.data$nCount_RNA
dataForPlot$nGene <- whole@meta.data$nFeature_RNA
dataForPlot$nUmiLog2 <- log2(whole@meta.data$nCount_RNA)
dataForPlot$nGeneLog2 <- log2(whole@meta.data$nFeature_RNA)

## FINDING MARKERS
whole.markers <- FindAllMarkers(object = whole,
                                only.pos = TRUE,
                                min.pct = 0.10,
                                thresh.use = 0.10)

# Rename Idents
# NB: These are named based on the output of FindAllMarkers using biological knowledge
whole <- RenameIdents(whole, '0' = "TN", '1' = "TTE", '2' = "TEM", '3' = "TCM", '4' = "MAIT")

# Check clusters look as per paper
DimPlot(whole, reduction = "umap", label = TRUE)


# 5. Perform integrated analysis to allow compatibility with ProjecTILs reference atlas ------------

# Split data into young and old to allow integration based on age
whole.list <- SplitObject(whole, split.by = "age")
whole.list <- lapply(X = whole.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Remove genes from var.features list based on authors selected parameters
whole.list$young@assays$RNA@var.features <-  
  whole.list$young@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^CD138-Antibody|^CD20-Antibody|^IgM-Antibody|^CD5-Antibody|^CD24-Antibody|^IgG-Antibody|^IgD-Antibody|^CD38-Antibody|^CD3-Antibody|^CD4-Antibody|^CD8-Antibody|^CD56-Antibody|^CD45RA-Antibody|^CD194-Antibody|^CD25-Antibody|^CD45RO-Antibody|^CD279-Antibody|^TIGIT-Antibody|^CD294-Antibody|^CD183-Antibody|^CD195-Antibody|^CD196-Antibody|^CD185-Antibody|^CD103-Antibody|^CD69-Antibody|^CD62L-Antibody|^CD197-Antibody|^CD152-Antibody|^CD223-Antibody|^KLRG1-Antibody|^CD27-Antibody|^CD314-Antibody|^CD57-Antibody|^CD366-Antibody|^CX3CR1-Antibody|^TCR-Antibody|^CD357-Antibody|^CD184-Antibody|^CD28-Antibody|^CD127-Antibody", whole.list$young@assays$RNA@var.features)]

whole.list$old@assays$RNA@var.features <-  
  whole.list$old@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^CD138-Antibody|^CD20-Antibody|^IgM-Antibody|^CD5-Antibody|^CD24-Antibody|^IgG-Antibody|^IgD-Antibody|^CD38-Antibody|^CD3-Antibody|^CD4-Antibody|^CD8-Antibody|^CD56-Antibody|^CD45RA-Antibody|^CD194-Antibody|^CD25-Antibody|^CD45RO-Antibody|^CD279-Antibody|^TIGIT-Antibody|^CD294-Antibody|^CD183-Antibody|^CD195-Antibody|^CD196-Antibody|^CD185-Antibody|^CD103-Antibody|^CD69-Antibody|^CD62L-Antibody|^CD197-Antibody|^CD152-Antibody|^CD223-Antibody|^KLRG1-Antibody|^CD27-Antibody|^CD314-Antibody|^CD57-Antibody|^CD366-Antibody|^CX3CR1-Antibody|^TCR-Antibody|^CD357-Antibody|^CD184-Antibody|^CD28-Antibody|^CD127-Antibody", whole.list$old@assays$RNA@var.features)]

# Select features that are repeatedly variable across datasets for integration
# NB: This portion of the script is directly taken from satijalab
features <- SelectIntegrationFeatures(object.list = whole.list)
immune.anchors <- FindIntegrationAnchors(object.list = whole.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering 
# NB: Set cluster resolution lower than default as these are CD8+ T-cells
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.15)

# Rename Seruat object so we can just copy ProjecTILS vignette going forward
data.seurat <- immune.combined

# Run FindAllMarkers to identify clusters - use same parameters as original authors
# NB: We will repeat this later on so the data is properly annotated
all_markers <- FindAllMarkers(object = data.seurat, 
                              assay = "RNA",
                              only.pos = TRUE,
                              min.pct = 0.10,
                              thresh.use = 0.10)

# NB: We will rename Idents later on so the DE script is easier to loop
# data.seurat <- RenameIdents(data.seurat, '0' = "TN", '1' = "TTE", '2' = "TEM", 
#'3' = "MAIT", '4' = "TCM")

# Check clustering
DimPlot(data.seurat, label = TRUE)


# 6. Save out data and print stats -----------------------------------------------------------------

saveRDS(data.seurat, "Data/R_out/Controls_3/integrated_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '7. Control dataset #3 - Complete workflow' script")