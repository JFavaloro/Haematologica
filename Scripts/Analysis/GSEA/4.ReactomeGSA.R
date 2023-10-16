# Script information -------------------------------------------------------------------------------

# Title: ReactomeGSA
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will use ReactomeGSA to determine enrichment across a number of 
# cellular processes using the vignette:
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html
# in conjunction with reactome.org.

# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '4. ReactomeGSA' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, ReactomeGSA)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load GSEA results
all_10x <- readRDS("Data/R_out/GSEA/all_10x.rds")
merged_list <- readRDS("Data/R_out/GSEA/merged_list.rds")


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 all_10x pre-processing -----------------------------------------------------------------------

# Split the all_10x object into BM and PB then recombine to allow comparisons
all_10x_BM <- subset(all_10x, subset = Tissue == "BM")
all_10x_PB <- subset(all_10x, subset = Tissue == "PB")

# Rename Idents
all_10x_BM <- RenameIdents(all_10x_BM, `TEM` = "T[EM]_BM", `TTE` = "T[TE]_BM", `TN` = "T[N]_BM", 
                          `Cyto_TEM` = "Cyto_T[EM]_BM", `PRE_EX` = "P[RE]_EX_BM", `TCM` = "T[CM]_BM")
all_10x_PB <- RenameIdents(all_10x_PB, `TEM` = "T[EM]_PB", `TTE` = "T[TE]_PB", `TN` = "T[N]_PB", 
                          `Cyto_TEM` = "Cyto_T[EM]_PB", `PRE_EX` = "P[RE]_EX_PB", `TCM` = "T[CM]_PB")

# Merge Seurat objects back together
all_10x <- merge(all_10x_BM, all_10x_PB)

# Normalise the data again after merge
all_10x <- NormalizeData(all_10x)


# 3.2 Control_data pre-processing ------------------------------------------------------------------

# 3.2A Remove data for clusters with insufficient data to compare and fix object for plotting ------
for (i in 1:length(merged_list)) {
  temp <- subset(merged_list[[i]], subset = (functional.cluster == "TN" | functional.cluster == "TCM" | functional.cluster == "TEM" | functional.cluster == "TTE"))
  merged_list[[i]] <- temp
}

#Reset levels of Sample in merged_list objects
merged_list$`Activated BM`@meta.data$Sample <- 
  factor(merged_list$`Activated BM`@meta.data$Sample, 
         levels = as.factor(c("MM Reference (BM)", "Activated BM")))
merged_list$`Resting BM`@meta.data$Sample <- 
  factor(merged_list$`Resting BM`@meta.data$Sample, 
         levels = as.factor(c("MM Reference (BM)", "Resting BM")))
merged_list$`Activated PB`@meta.data$Sample <- 
  factor(merged_list$`Activated PB`@meta.data$Sample, 
         levels = as.factor(c("MM Reference (PB)", "Activated PB")))
merged_list$`Resting PB`@meta.data$Sample <- 
  factor(merged_list$`Resting PB`@meta.data$Sample, 
         levels = as.factor(c("MM Reference (PB)", "Resting PB")))
merged_list$`Age-matched PB`@meta.data$Sample <- 
  factor(merged_list$`Age-matched PB`@meta.data$Sample, 
         levels = as.factor(c("MM Reference (PB)", "Age-matched Controls PB")))
merged_list$`Young PB`@meta.data$Sample <- 
  factor(merged_list$`Young PB`@meta.data$Sample, 
         levels = as.factor(c("MM Reference (PB)", "Young PB")))

# Set Idents to functional.cluster
merged_list <- map(merged_list, `Idents<-`, value = "functional.cluster")

# Reorder Idents for plotting
for (i in names(merged_list)) {
  merged_list[[`i`]]@meta.data$functional.cluster <- 
    factor(merged_list[[`i`]]@meta.data$functional.cluster, levels = c("TN", "TCM", "TEM", "TTE"))
}


# 3.2B Split the merged_list objects into reference and control  -----------------------------------

# Create empty lists to store data
merged_reference_BM <- list()
merged_control_BM <- list()
merged_reference_PB <- list()
merged_control_PB <- list()

# Subset data into reference and controls
for (i in names(merged_list)[1:2]) {
  merged_reference_BM[[`i`]] <- subset(merged_list[[`i`]], subset = Sample == "MM Reference (BM)")
}

for (i in names(merged_list)[1:2]) {
  merged_control_BM[[`i`]] <- subset(merged_list[[`i`]], subset = Sample != "MM Reference (BM)")
}

for (i in names(merged_list)[3:6]) {
  merged_reference_PB[[`i`]] <- subset(merged_list[[`i`]], subset = Sample == "MM Reference (PB)")
}

for (i in names(merged_list)[3:6]) {
  merged_control_PB[[`i`]] <- subset(merged_list[[`i`]], subset = Sample != "MM Reference (PB)")
}

# Rename Idents
for (i in names(merged_reference_BM)) {
  merged_reference_BM[[`i`]] <- RenameIdents(merged_reference_BM[[`i`]], 
                                             `TN` = "T[N]_BM_ref", `TCM` = "T[CM]_BM_ref", 
                                             `TEM` = "T[EM]_BM_ref", `TTE` = "T[TE]_BM_ref")
}

# Rename Idents #NB: No TN cells
for (i in names(merged_control_BM)) {
  merged_control_BM[[`i`]] <- RenameIdents(merged_control_BM[[`i`]], 
                                           `TN` = "T[N]_BM_control", `TCM` = "T[CM]_BM_control", 
                                           `TEM` = "T[EM]_BM_control", `TTE` = "T[TE]_BM_control")
}

# Rename Idents
for (i in names(merged_reference_PB)) {
  merged_reference_PB[[`i`]] <- RenameIdents(merged_reference_PB[[`i`]], 
                                             `TN` = "T[N]_PB_ref", `TCM` = "T[CM]_PB_ref", 
                                             `TEM` = "T[EM]_PB_ref", `TTE` = "T[TE]_PB_ref")
}

# Rename Idents
for (i in names(merged_control_PB)) {
  merged_control_PB[[`i`]] <- RenameIdents(merged_control_PB[[`i`]], 
                                           `TN` = "T[N]_PB_control", `TCM` = "T[CM]_PB_control", 
                                           `TEM` = "T[EM]_PB_control", `TTE` = "T[TE]_PB_control")
}


# 3.2C Merge Seurat objects back together ----------------------------------------------------------

merged_list$`Activated BM` <- merge(merged_reference_BM$`Activated BM`, 
                                    merged_control_BM$`Activated BM`)
merged_list$`Resting BM` <- merge(merged_reference_BM$`Resting BM`, 
                                  merged_control_BM$`Resting BM`)
merged_list$`Activated PB` <- merge(merged_reference_PB$`Activated PB`, 
                                    merged_control_PB$`Activated PB`)
merged_list$`Resting PB` <- merge(merged_reference_PB$`Resting PB`, 
                                  merged_control_PB$`Resting PB`)
merged_list$`Age-matched PB` <- merge(merged_reference_PB$`Age-matched PB`, 
                                      merged_control_PB$`Age-matched PB`)
merged_list$`Young PB` <- merge(merged_reference_PB$`Young PB`, 
                                merged_control_PB$`Young PB`)
# Normalise data after merge
merged_list <- map(merged_list, NormalizeData)


# 4. Run ReactomeGSA  ------------------------------------------------------------------------------

# 4.1 Run GSVA for all_10x object ------------------------------------------------------------------
# NB: 20220502 - fixing names is no longer required with ReactomeGSA v.1.10.0
options(reactome_gsa.url = "http://gsa.reactome.org")

all_10x_gsva_results <- analyse_sc_clusters(all_10x, verbose = TRUE, 
                                            create_reactome_visualization = TRUE)

# Fix names in gsva_result
#names(all_10x_gsva_results@results[["Seurat"]][["pathways"]]) <- 
#  gsub(x = names(all_10x_gsva_results@results[["Seurat"]][["pathways"]]), 
#       pattern = "\\X", replacement = "")  
#names(all_10x_gsva_results@results[["Seurat"]][["fold_changes"]]) <- 
#  gsub(x = names(all_10x_gsva_results@results[["Seurat"]][["fold_changes"]]), 
#       pattern = "\\X", replacement = "")  


# 4.2 Run GSVA for Controls ------------------------------------------------------------------------

# Create empty list to store results
merged_gsva_results <- list()

# Run a loop for running GSVA
for (i in names(merged_list)) {
  temp <- analyse_sc_clusters(merged_list[[`i`]], verbose = TRUE, 
                              create_reactome_visualization = TRUE)
  merged_gsva_results[[length(merged_gsva_results) +1]] <- temp
}

# Rename objects in merged_gsva_results
names(merged_gsva_results) = names(merged_list)


# Fix names in gsva_result
#for (i in names(merged_gsva_results)) {
#  names(merged_gsva_results[[`i`]]@results[["Seurat"]][["pathways"]]) <- 
#    gsub(x = names(merged_gsva_results[[`i`]]@results[["Seurat"]][["pathways"]]), 
#         pattern = "\\X", replacement = "")
#}

#for (i in names(merged_gsva_results)) {
#  names(merged_gsva_results[[`i`]]@results[["Seurat"]][["fold_changes"]]) <- 
#    gsub(x = names(merged_gsva_results[[`i`]]@results[["Seurat"]][["fold_changes"]]), 
#         pattern = "\\X", replacement = "")
#}


# 5. Save out data and print stats -----------------------------------------------------------------

# Save GSVA results out
saveRDS(all_10x_gsva_results, "Data/R_out/GSEA/all_10x_gsva_results.rds")
saveRDS(merged_gsva_results, "Data/R_out/GSEA/merged_gsva_results.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '4. ReactomeGSA' script")