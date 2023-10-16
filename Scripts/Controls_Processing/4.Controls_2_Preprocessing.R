# Script information -------------------------------------------------------------------------------

# Title: Control dataset #2 - Pre-processing
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will load a control dataset of BM from 32 individuals at various
# stages of PC dyscrasia, including an age-matched cohort from Zavidij et al., 2019: 
# https://doi.org/10.1016/j.clml.2019.09.040.
# This data is available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124310.
# As the data consists of BM-MNC, we will create some Seurat objects then use SingleR and celldex 
# to identify CD8+ T-cells and select these for further processing using the vignette:
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html.
# NB: Prior to running this script, the reference dataset must be downloaded to a cache folder.
# If this folder does not exist it will need to be created. Simply run the below code and type 'yes'
# library(celldex)
# ref <- BlueprintEncodeData()

# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '4. Control dataset #2 - Pre-processing' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, cowplot, SingleCellExperiment, SingleR, celldex)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Get all combinations of MGUS/SMM/MM/HD and 1:11 and discard unused combinations
all_samples_controls <- apply(expand.grid(c("HD", "MGUS", "MM", "SMM"), 1:11), 
                              1, paste0, collapse=" ") %>% str_replace_all(fixed(" "), "")
discard <- c("HD10", "HD11", 
             "MGUS6", "MGUS7", "MGUS8", "MGUS9", "MGUS10", "MGUS11", 
             "MM8", "MM9", "MM10", "MM11")
all_samples_controls <- setdiff(all_samples_controls, discard)

# Import GEX data and rename
all_gex_controls <- 
  lapply(all_samples_controls, 
         function(ID)  Read10X(data.dir = 
                                 str_glue('Data/Public_data/Controls_2/{ID}/filtered_gene_bc_matrices/GRCh38'))) 
names(all_gex_controls) <- all_samples_controls

# Define a function to pull metadata columns from our list of Seurat objects
get_metadata_column <- function(input_seurat_list, column_name){
  # pull all of the column from the Seurat object
  all_cols <- lapply(input_seurat_list, function(x) x[[]][,column_name])
  return(data.frame(lapply(all_cols, "length<-", max(lengths(all_cols)))) %>% gather) 
  # convert list to data.frame, then from wide to long format
}


# 3. Create Seurat objects -------------------------------------------------------------------------

# 3.1 Create Seurat Objects with default parameters ------------------------------------------------
# min.cells = gene must be observed in at least 3 cells
# min.features = cells must express at least 200 different genes
all_seurat_controls <- 
  lapply(1:length(all_samples_controls), 
         function(x) CreateSeuratObject(counts = all_gex_controls[[x]], 
                                        project = str_to_lower(all_samples_controls[x]), 
                                        min.cells = 3, min.features = 200))
names(all_seurat_controls) <- all_samples_controls


# 3.2 Fix cell barcodes to prevent conflict downstream and add some metadata -----------------------

# Add sample ID before cell barcode in "RNA" assay to match metadata
for (i in all_samples_controls) {
  all_seurat_controls[[i]] = RenameCells(all_seurat_controls[[i]], add.cell.id = i)
}

# Add barcode to metadata
for (i in all_samples_controls) {
  all_seurat_controls[[i]]@meta.data[["barcode"]] = 
    all_seurat_controls[[i]]@assays[["RNA"]]@data@Dimnames[[2]]
}

# Capitilize value in 'orig.ident' for consistency and prettier plotting
for (i in all_samples_controls) {
  all_seurat_controls[[i]]@meta.data[["orig.ident"]] = 
    toupper(all_seurat_controls[[i]]@meta.data[["orig.ident"]])
}

# Add disease state and patient ID to metadata - I'm sure there's a cleaner way to do this...
for (disease in c("HD")){
  for (patient in c(1:9)) {
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Disease"]] = disease
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Patient"]] = patient
  }
}

for (disease in c("MGUS")){
  for (patient in c(1:5)) {
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Disease"]] = disease
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Patient"]] = patient
  }
}

for (disease in c("MM")){
  for (patient in c(1:7)) {
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Disease"]] = disease
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Patient"]] = patient
  }
}

for (disease in c("SMM")){
  for (patient in c(1:11)) {
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Disease"]] = disease
    all_seurat_controls[[paste0(disease, patient)]]@meta.data[["Patient"]] = patient
  }
}

# Define function to add the mitochondrial %
add_mito <- function(input_seurat){
  input_seurat[["percent.mt"]] <- PercentageFeatureSet(input_seurat, pattern = "^MT-")
  return(input_seurat)
}

# Apply to all Seurat objects
all_seurat_controls <- lapply(all_seurat_controls, add_mito)


# 4. Apply some filtering to remove poor quality cells and isolate CD8+ T-cells --------------------

# 4.1 Define function to produce summary statistics (median and +/- 1, 1.5 & 2*sd) for plotting ----
# Get some QC metrics for plotting
nCount_RNA <- 
  get_metadata_column(input_seurat_list = all_seurat_controls, column_name = "nCount_RNA")
nFeature_RNA <- 
  get_metadata_column(input_seurat_list = all_seurat_controls, column_name = "nFeature_RNA")
percent.mt <- 
  get_metadata_column(input_seurat_list = all_seurat_controls, column_name = "percent.mt")

# Define function to produce summary statistics (median and +/- 1, 1.5 & 2*sd) for plotting
for (i in c(1,1.5,2)){
  data_summary <- function(x) {
    m = median(x)
    ymin = m-(i*sd(x))
    ymax = m+i*(sd(x))
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
         filename = paste0("Output/QC/Controls_2/Filter/QC_metrics_overview_", 
                           i,"SD", ".png"), width = 30, height = 10)
}


# 4.2 Filter Seurat objects ------------------------------------------------------------------------

# NB: Analysis suggests +/- 1.5 SD for all QC metrics looks appropriate
# Define function to filter Seurat objects based on the median +/- 1.5 sd for all three QC metrics
filter_seurat <- function(input_seurat){
  filtered_seurat = 
    subset(input_seurat, subset =
             nFeature_RNA > median(input_seurat$nFeature_RNA) - (1.5*sd(input_seurat$nFeature_RNA)) &
             nFeature_RNA < median(input_seurat$nFeature_RNA) + (1.5*sd(input_seurat$nFeature_RNA)) &
             nCount_RNA > median(input_seurat$nCount_RNA) - (1.5*sd(input_seurat$nCount_RNA)) &
             nCount_RNA < median(input_seurat$nCount_RNA) + (1.5*sd(input_seurat$nCount_RNA)) &
             percent.mt > median(input_seurat$percent.mt) - (1.5*sd(input_seurat$percent.mt)) &
             percent.mt < median(input_seurat$percent.mt) + (1.5*sd(input_seurat$percent.mt)))
  return(filtered_seurat)
}

# Apply filter
all_seurat_controls_filtered <- lapply(all_seurat_controls, filter_seurat)


# 4.3 Filter Seurat objects to select for CD8+ T-cells using SingleR -------------------------------
# NB: As SingleR relies on normalised counts, we will normalise data here

# 4.3.1 Normalise data -----------------------------------------------------------------------------
# Define function to run NormalizeData
norm_data <- function(input_seurat){
  input_seurat = NormalizeData(input_seurat)
  return(input_seurat)
}

# Apply to filtered Seurat objects
all_seurat_controls_filtered <- lapply(all_seurat_controls_filtered, norm_data)


# 4.3.2 Identify and select out CD8+ T-cells from our control dataset ------------------------------

# Load a reference dataset to annotate our data
ref <- celldex::BlueprintEncodeData()

# Create empty list to store our data
sce_list <- list()

# Convert our Seurat objects into SingleCellExperiment objects
for (i in names(all_seurat_controls_filtered)) {
  temp <- as.SingleCellExperiment(all_seurat_controls_filtered[[i]])
  sce_list[[i]] <- temp
}

# Create an empty list to store results from SingleR
predictions <- list()

# Loop over our control dataset to identify cell types
for (i in 1:length(all_samples_controls)){
  temp <- SingleR(test = sce_list[[i]], ref = ref, assay.type.test = 1, labels = ref$label.main)
  predictions[[length(predictions) +1]] <- temp
}

# Rename our annotated data
names(predictions) <- all_samples_controls

# Add cell type prediction to metadata
for (i in all_samples_controls) {
  all_seurat_controls_filtered[[i]]@meta.data[["cell_type"]] <- predictions[[i]]@listData[["labels"]]
}

# Create an empty list to store subsetted CD8+ T-cells
all_seurat_controls_filtered_CD8 <- list()

# Subset dataset to select only CD8+ T-cells
for (i in 1:length(all_samples_controls)){
  temp <- subset(all_seurat_controls_filtered[[i]], subset = cell_type == "CD8+ T-cells")
  all_seurat_controls_filtered_CD8[[length(all_seurat_controls_filtered_CD8) +1]] <- temp
}

# Rename our annotated data
names(all_seurat_controls_filtered_CD8) <- all_samples_controls

# Look at the number of cells in the objects before and after filtering

# Before filtering Seurat objects
sapply(all_seurat_controls, ncol)

#HD1 MGUS1   MM1  SMM1   HD2 MGUS2   MM2  SMM2   HD3 MGUS3   MM3  SMM3   HD4 MGUS4   MM4  SMM4   HD5
#130  1880  2187   297  1778   229  2047   280   960   741  1059   474   712   391  1374  2295   150
#MGUS5   MM5  SMM5   HD6   MM6  SMM6   HD7   MM7  SMM7   HD8  SMM8   HD9  SMM9 SMM10 SMM11 
#213   207  1024   164  1007   608  1164  1624   226   743   266   635   960  1130   237

# After filtering seurat
sapply(all_seurat_controls_filtered, ncol)

#HD1 MGUS1   MM1  SMM1   HD2 MGUS2   MM2  SMM2   HD3 MGUS3   MM3  SMM3   HD4 MGUS4   MM4  SMM4   HD5
#111  1586  1622   238  1520   165  1702   227   811   613   921   394   577   320  1104  1874   126
#MGUS5   MM5  SMM5   HD6   MM6  SMM6   HD7   MM7  SMM7   HD8  SMM8   HD9  SMM9 SMM10 SMM11 
#178   171   884   144   848   469   921  1248   202   625   221   529   818   973   204 

# After culling non-CD8+ T-cells
sapply(all_seurat_controls_filtered_CD8, ncol)

#HD1 MGUS1   MM1  SMM1   HD2 MGUS2   MM2  SMM2   HD3 MGUS3   MM3  SMM3   HD4 MGUS4   MM4  SMM4   HD5
#15   275   190    73   339    42   198    31   287    66   198   163   113   124   480   670     4 
#MGUS5   MM5  SMM5   HD6   MM6  SMM6   HD7   MM7  SMM7   HD8  SMM8   HD9  SMM9 SMM10 SMM11 
#32    28   365     7   231    89    73   256    33   102    52    40   182   259    47 


# 5. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_gex_controls, "Data/R_out/Controls_2/all_gex_controls.rds")
saveRDS(all_seurat_controls, "Data/R_out/Controls_2/all_seurat_controls.rds")
saveRDS(all_seurat_controls_filtered, "Data/R_out/Controls_2/all_seurat_controls_filtered.rds")
saveRDS(all_seurat_controls_filtered_CD8, "Data/R_out/Controls_2/all_seurat_controls_filtered_CD8.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '4. Control dataset #2 - Pre-processing' script")