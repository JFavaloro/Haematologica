# Script information -------------------------------------------------------------------------------

# Title: 10x ProjecTILs/SingleR annotations
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will use ProjecTILs (https://github.com/carmonalab/ProjecTILs) to 
# project our clustered data onto reference atlases of both tumour infiltrating lymphocytes (TILs)
# and a model of chronic infection with CMV, provided as part of the ProjecTILs package using the 
# vignette: https://carmonalab.github.io/ProjecTILs.demo/tutorial.html.
# Additionally, we will try reference datasets in conjunction the inference based package, SingleR 
# (https://bioconductor.org/packages/release/bioc/html/SingleR.html) to determine how our clustering
# compares to other established datasets and confirm how well we've annotated our data.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '1. 10x ProjecTILs/SingleR annotations' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, Matrix, ProjecTILs, SingleR, celldex, ensembldb, openxlsx, EnhancedVolcano, EnsDb.Hsapiens.v86, biomaRt)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load 10x data
all_list <- readRDS("Data/R_out/10x/all_list.rds")

# Load ProjecTILs reference atlases
ProjecTILs_ref <- 
  list("ref_TIL" = load.reference.map(ref = "Data/Public_data/ProjecTILs/ref_TIL_Atlas_mouse_v1.rds"), 
       "ref_CMV" = load.reference.map(ref = "Data/Public_data/ProjecTILs/ref_LCMV_Atlas_mouse_v1.rds"))

# Load reference datasets for SingleR - exclude Immgen until we can fix the issue with mouse genes
SingleR_ref <- 
  list("ref_HPCA" = celldex::HumanPrimaryCellAtlasData(ensembl=TRUE),
       "ref_BE" = celldex::BlueprintEncodeData(ensembl = TRUE),
       #"ref_Immgen" = celldex::ImmGenData(ensembl = TRUE),
       "ref_DICE" = celldex::DatabaseImmuneCellExpressionData(ensembl = TRUE),
       "ref_NH" = celldex::NovershternHematopoieticData(ensembl = TRUE),
       "ref_Monaco" = celldex::MonacoImmuneData(ensembl = TRUE))

# Create character vector for naming the SingleR annotations
SingleR_predictions_names <- 
  apply(expand.grid(c("ref_HPCA_","ref_BE_","ref_DICE_","ref_NH_", "ref_Monaco_"), 
                    c("all_10x", "BM_10x", "PB_10x")),
        1, paste, collapse = "")


# 3. Run ProjectTILs with authors reference atlases ------------------------------------------------

# Copy the all_list into a new list and set the Identity class to RNA
#NB: Default integrated assay doesn't allow comparisons as there are too few genes - use "RNA"
ProjecTILs_projections_10x <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Duplicate the objects in the ProjecTILs_projections list
ProjecTILs_projections_10x  <- list("all_10x_TIL" = ProjecTILs_projections_10x$all_10x,
                                    "all_10x_CMV" = ProjecTILs_projections_10x$all_10x,
                                    "BM_10x_TIL" = ProjecTILs_projections_10x$BM_10x,
                                    "BM_10x_CMV" = ProjecTILs_projections_10x$BM_10x,
                                    "PB_10x_TIL" = ProjecTILs_projections_10x$PB_10x,
                                    "PB_10x_CMV" = ProjecTILs_projections_10x$PB_10x)

# Duplicate the objects in the ProjecTILs_ref list to allow us to use map2
ProjecTILs_ref <- rep(ProjecTILs_ref, 3)

# Compare our data to ProjecTILs reference atlases
ProjecTILs_projections_10x <- 
  map2(ProjecTILs_projections_10x, ProjecTILs_ref,
       function(x, y){
         make.projection(query = x, ref = y, 
                         filter.cells = FALSE)
       }
  )

# Embed functional cluster data into projections
ProjecTILs_projections_10x <- 
  map2(ProjecTILs_projections_10x, ProjecTILs_ref,
       function(x, y){
         cellstate.predict(query = x, ref = y)
       }
  )


# 4. Run SingleR with authors reference atlases ----------------------------------------------------

# Copy the all_list into a new list and set the Identity class to RNA
all_list_sce <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Convert Seurat objects to SingleCellExperiment objects to allow SingleR to function
all_list_sce <- lapply(all_list_sce, as.SingleCellExperiment)

# Convert gene names of the SingleCellExperiment objects to Ensemble standards
for (i in 1:length(all_list_sce)){
  ens <- mapIds(EnsDb.Hsapiens.v86,
                keys = rownames(all_list_sce[[i]]),
                column = 'GENEID',
                keytype = 'SYMBOL')
  keep <- !is.na(ens)
  ens <- ens[keep]
  all_list_sce[[i]] <- all_list_sce[[i]][keep,]
  rownames(all_list_sce[[i]]) <- ens
  rm(ens, i, keep)
}

# Create an empty list to store predictions
SingleR_predictions_10x <- list()

# Compare our data to other datasets
for (i in names(all_list_sce)) {
  for (j in 1:length(SingleR_ref)){
    temp <- SingleR(test = all_list_sce[[i]], ref = SingleR_ref[[j]], 
                    labels = SingleR_ref[[j]]@colData@listData[["label.fine"]])
    SingleR_predictions_10x[[length(SingleR_predictions_10x) +1]] <- temp
  }
}

# Rename the SingleR_10x_predictions list
names(SingleR_predictions_10x) <- SingleR_predictions_names

# Copy the SingleR annotations back to our Seurat objects as metadata
for (i in names(SingleR_predictions_10x)[1:5]){
  all_list$all_10x[[i]] <- SingleR_predictions_10x[[i]]$labels
}

for (i in names(SingleR_predictions_10x)[6:10]){
  all_list$BM_10x[[i]] <- SingleR_predictions_10x[[i]]$labels
}

for (i in names(SingleR_predictions_10x)[11:15]){
  all_list$PB_10x[[i]] <- SingleR_predictions_10x[[i]]$labels
}


# 5. Save out data and print stats -----------------------------------------------------------------

# Save the ProjecTILs 10x projections and SingleR predictions
saveRDS(ProjecTILs_projections_10x, "Data/R_out/ProjecTILs/ProjecTILs_projections_10x.rds")
saveRDS(SingleR_predictions_10x, "Data/R_out/ProjecTILs/SingleR_predictions_10x.rds")
saveRDS(SingleR_ref, "Data/R_out/ProjecTILs/SingleR_ref.rds")
saveRDS(all_list, "Data/R_out/10x/all_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '1. 10x ProjecTILs/SingleR annotations' script")