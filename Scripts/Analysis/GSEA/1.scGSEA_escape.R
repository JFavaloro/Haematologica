# Script information -------------------------------------------------------------------------------

# Title: scGSEA
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will use escape to perform single cell gene set enrichment analysis
# (GSEA) across our complete data set to determine enrichment across a number of cellular processes
# using gene sets from the molecular signature database: https://www.gsea-msigdb.org/gsea/msigdb/ 
# and the vignette:
# http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '1. scGSEA with escape' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, escape, ProjecTILs)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed 10x data
all_list <- readRDS("Data/R_out/10x/all_list.rds")

# Load custom reference atlases
all_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/all_10x_atlas.rds")
BM_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/BM_10x_atlas.rds")
PB_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/PB_10x_atlas.rds")

# Load the Controls projections
Controls_query_list <- readRDS("Data/R_out/ProjecTILs/Controls_query_list.rds")

# Define function to run enrichIt across all objects in the merged_list
run_enrichIT <- function(input_object){
  input_object <- enrichIt(input_object, gene.sets = GS, groups = 1000, cores = 8)
  return(input_object)
}


# 3. Pre-processing --------------------------------------------------------------------------------

# Extract the all_10x object and set default assay to "RNA"
all_10x <- all_list$all_10x
DefaultAssay(all_10x) <- "RNA"

# Move atlases into a list
all_atlas_list <- list("all_10x_atlas" = all_10x_atlas,
                       "BM_10x_atlas" = BM_10x_atlas,
                       "PB_10x_atlas" = PB_10x_atlas)

# Create an empty list to store merged objects
merged_list <- list()

# Merge objects to allow comparisons across samples
merged_list$`Activated BM` <- 
  merge(all_atlas_list$BM_10x_atlas, Controls_query_list$Controls_1$`Activated BM`)
merged_list$`Resting BM` <- 
  merge(all_atlas_list$BM_10x_atlas, Controls_query_list$Controls_1$`Resting BM`)
merged_list$`Activated PB` <- 
  merge(all_atlas_list$PB_10x_atlas, Controls_query_list$Controls_1$`Activated PB`)
merged_list$`Resting PB` <- 
  merge(all_atlas_list$PB_10x_atlas, Controls_query_list$Controls_1$`Resting PB`)
merged_list$`Age-matched PB` <- 
  merge(all_atlas_list$PB_10x_atlas, Controls_query_list$Controls_3$`Age-matched Controls PB`)
merged_list$`Young PB` <- 
  merge(all_atlas_list$PB_10x_atlas, Controls_query_list$Controls_3$`Young PB`)

# Fix merged objects to allow comparisons across samples
merged_list <- map(merged_list, `DefaultAssay<-`, value = "RNA")
merged_list <- map(merged_list, `Idents<-`, value = "functional.cluster")
merged_list <- map(merged_list, NormalizeData)


# 4. Run Escape  -----------------------------------------------------------------------------------

# Get Hallmark genesets from MSigDB
GS <- getGeneSets(library = "H")

# Run scGSEA across all cells in our combined seurat object NB: This will take a while!
ES_all_10x <- enrichIt(all_10x, gene.sets = GS, groups = 1000, cores = 8)

# Apply run_enrichIT function to merged_list NB: This will take a couple of hours!
ES_merged_list <- lapply(merged_list, run_enrichIT)

# Tidy up the names to make them more presentable for plotting
names(ES_all_10x) <- gsub("HALLMARK_", "", names(ES_all_10x)) %>% str_to_sentence()

for (i in names(ES_merged_list)) {
  names(ES_merged_list[[i]]) <- 
    gsub("HALLMARK_", "", names(ES_merged_list[[i]])) %>% str_to_sentence()
}

# Append the GSEA results back to the all_10x Seurat object as metadata
all_10x <- AddMetaData(all_10x, ES_all_10x)

# Append the GSEA results back to the merged_list objects as metadata
merged_list <- map2(merged_list, ES_merged_list,
                    function(x, y){
                      AddMetaData(
                        object = x,
                        metadata = y
                      )
                    }
)


# 5. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_10x, "Data/R_out/GSEA/all_10x.rds")
saveRDS(merged_list, "Data/R_out/GSEA/merged_list.rds")
saveRDS(ES_all_10x, "Data/R_out/GSEA/ES_all_10x.rds")
saveRDS(ES_merged_list, "Data/R_out/GSEA/ES_merged_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '1. scGSEA with escape' script")