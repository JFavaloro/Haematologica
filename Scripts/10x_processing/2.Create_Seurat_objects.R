# Script information -------------------------------------------------------------------------------

# Title: Create Seurat object
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will create Seurat objects for each sample and then add the V(D)J 
# data so that everything is conveniently stored in one object.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '2. Create Seurat objects' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Get all combinations of PB/BM and 13/31/43/63
all_samples <- apply(expand.grid(c("BM", "PB"), c("13","31","43","63")), 1, paste, collapse="")

# Import GEX data and V(D)J data
all_gex <- readRDS("Data/R_out/10x/all_gex.rds")
all_vdj <- readRDS("Data/R_out/10x/all_vdj.rds")


# 3. Create Seurat objects and add additional data for QC ------------------------------------------

# Create Seurat Objects with default parameters:
# min.cells = gene must be observed in at least 3 cells
# min.features = cells must express at least 200 different genes
# NB: Have run this without additional filtering arguments to allow visualisation of all raw data
all_seurat <- lapply(1:length(all_samples), 
                     function(x) CreateSeuratObject(counts = 
                                                      all_gex[[x]], project = 
                                                      str_to_upper(all_samples[x]), 
                                                    min.cells = 3, min.features = 200))

# Add mitochondrial counts
# Define function to add the mitochondrial %
add_mito <- function(input_seurat){
  input_seurat[["percent.mt"]] <- PercentageFeatureSet(input_seurat, pattern = "^MT-")
  return(input_seurat)
}

# Apply to all Seurat objects
all_seurat <- lapply(all_seurat, add_mito)


# 4. Add V(D)J data --------------------------------------------------------------------------------

# Add patient specific V(D)J data to Seurat objects
# NB: At this time point there are no prefixes to any barcodes
# NB: The values for clonal bins have been determined after QC
# NB: The function combineExpression has changed to now calculate proportion - downsampling is no
# longer required to create clonal bins and this can now be determine by frequency. This will also
# error including the groupBy argument - removed. 2022-02-02
all_seurat <- 
  lapply(1:length(all_samples),
         function(x) combineExpression(all_vdj[[x]],
                                       all_seurat[[x]],
                                       cloneCall = "gene+nt",
                                       cloneTypes = c(None = 0, Small = 0.001, Medium = 0.01, 
                                                      Large = 0.1, Expanded = 1)))

# This can be used to classify cells based on pre-defined bins as originally performed.
# all_seurat <- lapply(1:length(all_samples), function(x) combineExpression(all_vdj[[x]], 
# all_seurat[[x]], cloneCall = "gene+nt", proportion = FALSE, groupBy = "sample", cloneTypes = 
# c(None = 0, Single = 1, Small = 6, Medium = 63, Large = 628, Expanded = 6275)))

# Name Seurat objects 
names(all_seurat) <- all_samples

# Fix names for "cloneType" to allow easier manipulation of data
for (i in 1:length(all_seurat)) {
  all_seurat[[i]]@meta.data <- 
    all_seurat[[i]]@meta.data %>% 
    separate(cloneType, into = "cloneType", sep = "\\(")
}


# 5. Rename cell barcodes and add additional meta data to Seurat object ----------------------------
# NB: Attempting this earlier results in an error when merging the information in step 4

# Add sample ID before cell barcode in metadata to prevent conflict
for (i in all_samples) {
  all_seurat[[i]]@meta.data[["barcode"]] <- 
    paste0(as.character(i), sep = "_", all_seurat[[i]]@meta.data[["barcode"]])
}

# Add sample ID before cell barcode in "RNA" assay to match metadata
for (i in all_samples) {
  all_seurat[[i]] <- RenameCells(all_seurat[[i]], add.cell.id = i)
}

# Add in Tissue and Patient metadata
for (Tissue in c("BM", "PB")){
  for (Patient in c(13, 31, 43, 63)) {
    all_seurat[[paste0(Tissue, Patient)]]@meta.data[["Tissue"]] <- Tissue
    all_seurat[[paste0(Tissue, Patient)]]@meta.data[["Patient"]] <- Patient
  }
}


# 5. Check how many cells have matching V(D)J data -------------------------------------------------

# This is # NA, Small, Medium, Large, Expanded in the order of the all_seurat object
value_filtered = c(2083, 4700, 2645, 572, 0, 1858, 4905, 2245, 977, 0, 
                   234, 1559, 962, 326, 0, 2509, 5258, 1343, 522, 0,
                   829, 5220, 1398, 1215, 0, 1192, 3102, 454, 1120, 2172,
                   1138, 3735, 3158, 1969, 0, 799, 1494, 2133, 2335, 751)

# These are the values when creating the Seurat objects with no filtering
value_unfiltered = c(2083, 4700, 2645, 572, 0, 1873, 4905, 2245, 977, 0, 
                     3839, 3644, 1794, 723, 0, 2877, 5258, 1343, 522, 0, 
                     2165, 5221, 1398, 1216, 0, 3148, 3105, 454, 1120, 2173, 
                     1138, 3735, 3158, 1969, 0, 3285, 1494, 2135, 2335, 751)

# Create some character vectors to plot
Sample = c(rep("BM13", 5), rep("PB13", 5), 
           rep("BM31", 5), rep("PB31", 5), 
           rep("BM43", 5), rep("PB43", 5), 
           rep("BM63", 5), rep("PB63", 5) )
expans_level = rep(c("NA", "Small", "Medium", "Large", "Expanded") , 8)

# Plot - Fix up the y-axis title
p1 <- data.frame(Sample,expans_level,value_filtered) %>%
  mutate(expans_level = fct_relevel(expans_level,"Expanded","Large","Medium","Small","NA"))
p1 = ggplot(p1, aes(fill=expans_level, y=value_filtered, x=Sample)) + 
  geom_bar(position="stack", stat="identity")
p1
ggsave(plot = p1, filename = "Output/QC/10x/Cells_with_matching_VDJ_filtered.png")

p1 <- data.frame(Sample,expans_level,value_unfiltered) %>%
  mutate(expans_level = fct_relevel(expans_level,"Expanded","Large","Medium","Small","NA"))
p1 = ggplot(p1, aes(fill=expans_level, y=value_unfiltered, x=Sample)) + 
  geom_bar(position="stack", stat="identity")
p1
ggsave(plot = p1, filename = "Output/QC/10x/Cells_with_matching_VDJ_unfiltered.png")


# 6. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_seurat, "Data/R_out/10x/all_seurat.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '2. Create Seurat objects' script")