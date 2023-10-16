# Script information -------------------------------------------------------------------------------

# Title: Projection of public data on 10x Reference atlases
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will use ProjectTILs (https://github.com/carmonalab/ProjecTILs) to
# project public data onto the custom reference atlases we created in the previous scripts. We will 
# then output projections for downstream analysis.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '6. Project public data on 10x Reference atlases' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire, openxlsx, multtest, metap, EnhancedVolcano, ProjecTILs)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed public data
Controls_1_complete <- readRDS("Data/R_out/Controls_1/Controls_1_complete.rds")
Controls_2_complete <- readRDS("Data/R_out/Controls_2/Controls_2_complete.rds")
Controls_3_complete <- readRDS("Data/R_out/Controls_3/Controls_3_complete.rds")

# Load custom reference atlases
all_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/all_10x_atlas.rds")
BM_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/BM_10x_atlas.rds")
PB_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/PB_10x_atlas.rds")


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 Subset Controls objects ----------------------------------------------------------------------

# 3.1A Subset Controls_1 object for testing  and add additional metadata ---------------------------
query_list_Controls_1 <- 
  list("Controls 1 complete" = Controls_1_complete,
       "Activated" = subset(Controls_1_complete, subset = stimulation_status == "act"),
       "Resting" = subset(Controls_1_complete, subset = stimulation_status == "rest"),
       "Activated BM" = subset(Controls_1_complete, subset = Tissue == "BM" & 
                                 stimulation_status == "act"),
       "Resting BM" = subset(Controls_1_complete, subset = Tissue == "BM" & 
                               stimulation_status == "rest"),
       "Activated PB" = subset(Controls_1_complete, subset = Tissue == "PB" & 
                                 stimulation_status == "act"),
       "Resting PB" = subset(Controls_1_complete, subset = Tissue == "PB" & 
                               stimulation_status == "rest"))

query_list_Controls_1[["Controls 1 complete"]][["Sample"]] = "Controls 1 complete"
query_list_Controls_1[["Activated"]][["Sample"]] = "Activated"
query_list_Controls_1[["Resting"]][["Sample"]] = "Resting"
query_list_Controls_1[["Activated BM"]][["Sample"]] = "Activated BM"
query_list_Controls_1[["Resting BM"]][["Sample"]] = "Resting BM"
query_list_Controls_1[["Activated PB"]][["Sample"]] = "Activated PB"
query_list_Controls_1[["Resting PB"]][["Sample"]] = "Resting PB"

# 3.1B Subset Controls_2 object for testing and add additional metadata ----------------------------
query_list_Controls_2 <- 
  list("Controls 2 complete" = Controls_2_complete,
       "Age-matched Controls BM" = subset(Controls_2_complete, subset = Disease == "HD"),
       "MGUS BM" = subset(Controls_2_complete, subset = Disease == "MGUS"),
       "SMM BM" = subset(Controls_2_complete, subset = Disease == "SMM"),
       "MM BM" = subset(Controls_2_complete, subset = Disease == "MM"))

query_list_Controls_2[["Controls 2 complete"]][["Sample"]] = "Controls 2 complete"
query_list_Controls_2[["Age-matched Controls BM"]][["Sample"]] = "Age-matched Controls BM"
query_list_Controls_2[["MGUS BM"]][["Sample"]] = "MGUS BM"
query_list_Controls_2[["SMM BM"]][["Sample"]] = "SMM BM"
query_list_Controls_2[["MM BM"]][["Sample"]] = "MM BM"

# 3.1C Subset Controls_3 object for testing --------------------------------------------------------
query_list_Controls_3 <- 
  list("Controls 3 complete" = Controls_3_complete,
       "Age-matched Controls PB" = subset(Controls_3_complete, subset = age == "old"),
       "Young PB" = subset(Controls_3_complete, subset = age == "young"))

query_list_Controls_3[["Controls 3 complete"]][["Sample"]] = "Controls 3 complete"
query_list_Controls_3[["Age-matched Controls PB"]][["Sample"]] = "Age-matched Controls PB"
query_list_Controls_3[["Young PB"]][["Sample"]] = "Young PB"


# 3.2 Set correct parameters for projection --------------------------------------------------------

# Switch Identity class of query lists to "seurat_clusters"
query_list_Controls_1 <- map(query_list_Controls_1, `Idents<-`, value = "seurat_clusters")
query_list_Controls_2 <- map(query_list_Controls_2, `Idents<-`, value = "seurat_clusters")
query_list_Controls_3 <- map(query_list_Controls_3, `Idents<-`, value = "seurat_clusters")

# Switch default assay of query lists to "RNA" to permit projections
#NB: Default integrated assay doesn't allow comparisons as there are too few genes - use "RNA"
query_list_Controls_1 <- map(query_list_Controls_1, `DefaultAssay<-`, value = "RNA")
query_list_Controls_2 <- map(query_list_Controls_2, `DefaultAssay<-`, value = "RNA")
query_list_Controls_3 <- map(query_list_Controls_3, `DefaultAssay<-`, value = "RNA")

# 3.3 Prepare lists of reference atlases to allow use of map2 --------------------------------------

Controls_1_ref <- list("all_10x_atlas" = all_10x_atlas,
                       "all_10x_atlas" = all_10x_atlas,
                       "all_10x_atlas" = all_10x_atlas,
                       "BM_10x_atlas" = BM_10x_atlas,
                       "BM_10x_atlas" = BM_10x_atlas,
                       "PB_10x_atlas" = PB_10x_atlas,
                       "PB_10x_atlas" = PB_10x_atlas)

Controls_2_ref <- list("BM_10x_atlas" = BM_10x_atlas) %>% rep(5)

Controls_3_ref <- list("PB_10x_atlas" = PB_10x_atlas) %>% rep(3)

# 3.4 Move lists into lists ------------------------------------------------------------------------

# Move all queries into a list
query_list <- list("Controls_1" = query_list_Controls_1,
                   "Controls_2" = query_list_Controls_2,
                   "Controls_3" = query_list_Controls_3)

# Move all reference atlases into a list
ref_list <- list("Controls_1_ref" = Controls_1_ref,
                 "Controls_2_ref" = Controls_2_ref,
                 "Controls_3_ref" = Controls_3_ref)

# 4. Run ProjectTILs -------------------------------------------------------------------------------

# Compare our data to other datasets
Controls_query_list <- map2(query_list, ref_list,
                            function(query, reference) {
                              temp <- map2(query, reference,
                                           function(x, y){
                                             make.projection(
                                               query = x, 
                                               ref = y, 
                                               filter.cells = FALSE
                                             )
                                           })
                              
                              # Embed functional cluster data into projections
                              temp <- map2(temp, reference,
                                           function(x, y){
                                             cellstate.predict(
                                               query = x, 
                                               ref = y)
                                           }
                              )
                              temp
                            })


# 5. Save out data and print stats -----------------------------------------------------------------

# Save embeded data
saveRDS(Controls_query_list, "Data/R_out/ProjecTILs/Controls_query_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '6. Project public data on 10x Reference atlases' script")