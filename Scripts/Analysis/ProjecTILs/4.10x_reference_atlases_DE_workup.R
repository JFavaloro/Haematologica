# Script information -------------------------------------------------------------------------------

# Title: ProjecTILs/10x Reference atlases DE workup
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will largely recreate our initial DE work up comparing between a 
# number of clusters as found by clustering using the recalculated PC. 
# NB: All processed data has been coerced into 7 clusters allowing us to easily loop over this
# workflow with all data. This is then changed to 6 clusters, then 5 clusters after analysis and 
# the procedure repeated. After analysis the "correct" number of clusters will be determined in the 
# next script and the Seurat objects exported for further analysis.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '4. ProjecTILs/10x Reference atlases DE workup' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, Matrix, ProjecTILs, openxlsx, EnhancedVolcano)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load data for analysis
all_refs <- readRDS("Data/R_out/ProjecTILs/all_refs.rds")

# Load additional resources
black_list <- readRDS("Data/R_out/black_list.rds")
annotations <- read.csv("Data/Public_data/annotations.csv")

# Modify black_list to retain TCR genes
black_list_TCR <- 
  black_list[!grepl("TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC", black_list)]

# Define %!in% function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Define function to remove black_listed genes from discriminant genes
remove_genes = function(x){
  x = x[!((x$gene) %in% black_list),]
  return(x)
}

# Define function for annotating markers
# NB: Pass to unique at the end removes redundancy but also seems to remove some genes (?)
annotate_markers <- function(input_markers){
  input_markers <- inner_join(x = input_markers, 
                              y = annotations[, c("gene_name", "description")],
                              by = c("gene" = "gene_name")) %>% unique()
  return(input_markers)
}

# Create character vector of samples names
sample_names <- c("all_10x", "BM_10x", "PB_10x")

# Define cluster names for naming of lists - 7 clusters
cluster_names_7 <- apply(expand.grid(c("cluster_"), c(0:6)), 1, paste, collapse="") %>% 
  rep(3) %>% paste0(c(rep("_all_10x_atlas", 7), rep("_BM_10x_atlas", 7), rep("_PB_10x_atlas", 7)))

# Define cluster names for naming of lists - 6 clusters
cluster_names_6 <- apply(expand.grid(c("cluster_"), c(0:5)), 1, paste, collapse="") %>% 
  rep(3) %>% paste0(c(rep("_all_10x_atlas", 6), rep("_BM_10x_atlas", 6), rep("_PB_10x_atlas", 6)))

# Define cluster names for naming of lists - 5 clusters
cluster_names_5 <- apply(expand.grid(c("cluster_"), c(0:4)), 1, paste, collapse="") %>% 
  rep(3) %>% paste0(c(rep("_all_10x_atlas", 5), rep("_BM_10x_atlas", 5), rep("_PB_10x_atlas", 5)))

# Create character vectors of colours for plotting clusters
scales::show_col(scales::hue_pal()(6))
scales::show_col(scales::hue_pal()(12))
scales::show_col(scales::hue_pal()(25))
ccolours_7 <- c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")
ccolours_6 <- c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF")
ccolours_5 <- c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB")


# 3. Pre-processing --------------------------------------------------------------------------------

# Copy the all_refs object into a new object to allow reuse of existing code
all_list <- all_refs

# Switch default assay to "RNA" for differential expression testing
all_list <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Remove non helpful genes from DE results #NB: This will remove genes from the object
for (i in names (all_list)) {
  counts <- GetAssayData(all_list[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  all_list[[i]] <- subset(all_list[[i]], features = rownames(counts))
}


# 4. Analysis of 7 cluster Seurat objects ----------------------------------------------------------
# 4.A DE testing - FindAllMarkers

# 4.A.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

# Create an empty list to store DE results
FindAllMarkers_All_list_7 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_list)){
  temp <- FindAllMarkers(all_list[[i]])
  FindAllMarkers_All_list_7[[length(FindAllMarkers_All_list_7) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_All_list_7) <- sample_names

# Annotate the results
FindAllMarkers_All_list_7 <- lapply(FindAllMarkers_All_list_7, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_All_list_7,
           "Output/DE/10x_atlases/Workup/7_clusters/FindAllMarkers_All_7_clusters.xlsx")

# 4.A.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 
# Copy the all_10x object and switch the identity class to "Tissue"
all_10x_Tissue <- all_list$all_10x
Idents(all_10x_Tissue) <- "Tissue"

# Run FindAllMarkers
FindAllMarkers_Tissue_list_7 <- FindAllMarkers(all_10x_Tissue)

# Annotate the results
FindAllMarkers_Tissue_list_7 <- inner_join(x = FindAllMarkers_Tissue_list_7, 
                                           y = annotations[, c("gene_name", "description")],
                                           by = c("gene" = "gene_name")) %>% unique()

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = BM, pct.2 = PB
FindAllMarkers_Tissue_list_7 <- 
  FindAllMarkers_Tissue_list_7[!(FindAllMarkers_Tissue_list_7$cluster == "PB"),]

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list_7,
           "Output/DE/10x_atlases/Workup/7_clusters/FindAllMarkers_Tissue_7_clusters.xlsx")

# Output Volcano plots for Tissue restricted differences
EnhancedVolcano(FindAllMarkers_Tissue_list_7, 
                lab = FindAllMarkers_Tissue_list_7$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "BM vs. PB",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 12, pointSize = 2, labSize = 4, legendLabSize = 12,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/7_clusters/FindAllMarkers_Tissue_",
                "BM vs. PB", ".png"), scale = 2:1)

# 4.A.iii DE testing - FindAllMarkers_clone - Assess global differences between clonal bins --------

# Copy the all_10x object and switch the identity class to "cloneType"
all_10x_Clone <- all_list
all_10x_Clone <- map(all_10x_Clone, `Idents<-`, value = "cloneType")

# Create an empty list to store DE results
FindAllMarkers_Clone_list_7 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_10x_Clone)){
  temp <- FindAllMarkers(all_10x_Clone[[i]])
  FindAllMarkers_Clone_list_7[[length(FindAllMarkers_Clone_list_7) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Clone_list_7) <- sample_names

# Annotate the results
FindAllMarkers_Clone_list_7 <- lapply(FindAllMarkers_Clone_list_7, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Clone_list_7, 
           "Output/DE/10x_atlases/Workup/7_clusters/FindAllMarkers_Clone_7_clusters.xlsx")

# 4.A.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Copy the all_10x list and switch the identity class to "Patient"
all_10x_Patient <- all_list
all_10x_Patient <- map(all_10x_Patient, `Idents<-`, value = "Patient")

# Create an empty list to store DE results
FindAllMarkers_Patient_list_7 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_10x_Patient)){
  temp = FindAllMarkers(all_10x_Patient[[i]])
  FindAllMarkers_Patient_list_7[[length(FindAllMarkers_Patient_list_7) +1]] = temp
}

# Rename the objects in the list
names(FindAllMarkers_Patient_list_7) <- sample_names

# Annotate the results 
FindAllMarkers_Patient_list_7 <- lapply(FindAllMarkers_Patient_list_7, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_7, 
           "Output/DE/10x_atlases/Workup/7_clusters/FindAllMarkers_Patient_7_clusters.xlsx")


# 4.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 4.B.i FindMarkers - Itemise differences between Clusters -----------------------------------------

# Create an empty list to store results in
FindMarkers_All_list_7 <- list()

# Run for loop to iterate across all clusters
# NB: This will fail unless each input sample is 7 clusters
for (i in names(all_list)){
  for (j in 1:7) {
    temp <- FindMarkers(all_list[[i]], ident.1 = j-1)
    FindMarkers_All_list_7[[length(FindMarkers_All_list_7) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list. NB: These need to be unique or export to *.xlsx will fail
names(FindMarkers_All_list_7) <- cluster_names_7

# Save out the DE results
write.xlsx(FindMarkers_All_list_7, 
           "Output/DE/10x_atlases/Workup/7_clusters/FindMarkers_All_7_clusters.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_All_list_7)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_All_list_7[[`i`]], 
                     lab = FindMarkers_All_list_7[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_All_list_7)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/7_clusters/",
                                    names(FindMarkers_All_list_7)[i], ".png"), scale = 2:1)
}

# 4.B.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create an empty list to store results in
FindMarkers_Tissue_list_7 <- list()

# Run for loop to iterate across all clusters
for (i in 1:7) {
  temp <- FindMarkers(all_list$all_10x, ident.1 = "BM", group.by = "Tissue", subset.ident = i-1)
  FindMarkers_Tissue_list_7[[i]] <- temp %>% rownames_to_column %>% rename(gene = rowname) %>%
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>%
    unique()
}

# Rename the items in the list
names(FindMarkers_Tissue_list_7) <- cluster_names_7[1:7]

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list_7, 
           "Output/DE/10x_atlases/Workup/7_clusters/FindMarkers_Tissue_7_clusters.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Tissue_list_7)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Tissue_list_7[[`i`]], 
                     lab = FindMarkers_Tissue_list_7[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Tissue_list_7)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/7_clusters/FindMarkers_Tissue_",
                                   names(FindMarkers_Tissue_list_7)[i], ".png"), scale = 2:1)
}


# 4.C Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 4.C.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Tissue_list_7 <- list()

# Run for loop to iterate across all clusters
for (i in 1:7) {
  temp <- FindConservedMarkers(all_list$all_10x, ident.1 = i-1, grouping.var = "Tissue")
  FindConservedMarkers_Tissue_list_7[[i]] <- temp %>% 
    rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>% 
    unique()
}

# Rename the items in the list
names(FindConservedMarkers_Tissue_list_7) <- cluster_names_7[1:7]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_7, 
           "Output/DE/10x_atlases/Workup/7_clusters/FindConservedMarkers_Tissue_7_clusters.xlsx")


# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_7) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/10x_atlases/Workup/7_clusters/",
                  sample_names[i], ".png"))
}

# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Calculate distribution of clusters and move results into a temp_list
temp_list <- 
  list("temp_all_10x" = as.data.frame(table(all_list$all_10x$seurat_clusters, 
                                            all_list$all_10x$orig.ident)),
       "temp_BM_10x" = as.data.frame(table(all_list$BM_10x$seurat_clusters, 
                                           all_list$BM_10x$orig.ident)),
       "temp_PB_10x" = as.data.frame(table(all_list$BM_10x$seurat_clusters, 
                                           all_list$BM_10x$orig.ident)))

# Output barplots of cluster distribution
for (i in 1:length(sample_names)) {
  ggplot(temp_list[[i]], 
         aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = rev(ccolours_7), name = "Cluster",
                      (rev(levels(all_list[[i]]@meta.data$seurat_clusters)))) + 
    labs(x = "Sample", y = "Cells", title = sample_names[i]) + theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Barplots/10x_atlases/Workup/7_clusters/",
                  sample_names[i], ".png"))
}


# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all_10x object and output a Dotplot separated by tissue
top_10_markers_all_10x <- 
  FindAllMarkers_All_list_7$all_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all_10x <- unique(top_10_markers_all_10x$gene)
DotPlot(all_list$all_10x, 
        features = top_10_unique_all_10x, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/7_clusters/",
                "BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the BM_10x object and output a Dotplot separated by patient
top_10_markers_BM_10x <- 
  FindAllMarkers_All_list_7$BM_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM_10x <- unique(top_10_markers_BM_10x$gene)
DotPlot(all_list$BM_10x, 
        features = top_10_unique_BM_10x, split.by = "Patient", cols = c("green", "purple")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "BM_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/7_clusters/",
                "BM pt43 vs. pt 63", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB_10x object and output a Dotplot separated by patient
top_10_markers_PB_10x <- 
  FindAllMarkers_All_list_7$PB_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB_10x <- unique(top_10_markers_PB_10x$gene)
DotPlot(all_list$PB_10x, 
        features = top_10_unique_PB_10x, split.by = "Patient", 
        cols = c("orange", "turquoise", "green", "purple")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "PB_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/7_clusters/",
                "PB all patients", ".png"), scale = 2:1)


# 5. Analysis of 6 cluster Seurat objects ----------------------------------------------------------

# 5.A Set resolutions of Seurat objects to allow clustering into 6 clusters
# Set resolution to 0.32 for all_10x object
all_list$all_10x$seurat_clusters <- 
  all_list$all_10x$integrated_snn_res.0.32
Idents(all_list$all_10x) <- "seurat_clusters"

# Set resolution to 0.36 for BM_10x object
all_list$BM_10x$seurat_clusters <- 
  all_list$BM_10x$integrated_snn_res.0.36
Idents(all_list$BM_10x) <- "seurat_clusters"

# Set resolution to 0.25 for PB_10x object
all_list$PB_10x$seurat_clusters <- 
  all_list$PB_10x$integrated_snn_res.0.25
Idents(all_list$PB_10x) <- "seurat_clusters"


# 5.B DE testing - FindAllMarkers ------------------------------------------------------------------

# 5.B.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

# Create an empty list to store DE results
FindAllMarkers_All_list_6 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_list)){
  temp <- FindAllMarkers(all_list[[i]])
  FindAllMarkers_All_list_6[[length(FindAllMarkers_All_list_6) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_All_list_6) <- sample_names

# Annotate the results 
FindAllMarkers_All_list_6 <- lapply(FindAllMarkers_All_list_6, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_All_list_6,
           "Output/DE/10x_atlases/Workup/6_clusters/FindAllMarkers_All_6_clusters.xlsx")

# 5.B.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 
# Copy the all_10x object and switch the identity class to "Tissue"
all_10x_Tissue <- all_list$all_10x
Idents(all_10x_Tissue) <- "Tissue"

# Run FindAllMarkers
FindAllMarkers_Tissue_list_6 <- FindAllMarkers(all_10x_Tissue)

# Annotate the results
FindAllMarkers_Tissue_list_6 <- inner_join(x = FindAllMarkers_Tissue_list_6, 
                                           y = annotations[, c("gene_name", "description")],
                                           by = c("gene" = "gene_name")) %>% unique()

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = BM, pct.2 = PB
FindAllMarkers_Tissue_list_6 <- 
  FindAllMarkers_Tissue_list_6[!(FindAllMarkers_Tissue_list_6$cluster == "PB"),]

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list_6,
           "Output/DE/10x_atlases/Workup/6_clusters/FindAllMarkers_Tissue_6_clusters.xlsx")

# Output Volcano plots for Tissue restricted differences
EnhancedVolcano(FindAllMarkers_Tissue_list_6, 
                lab = FindAllMarkers_Tissue_list_6$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "BM vs. PB",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 12, pointSize = 2, labSize = 4, legendLabSize = 12,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/6_clusters/FindAllMarkers_Tissue_",
                "BM vs. PB", ".png"), scale = 2:1)

# 5.B.iii DE testing - FindAllMarkers_clone - Assess global differences between clonal bins --------

# Copy the all_10x object and switch the identity class to "cloneType"
all_10x_Clone <- all_list
all_10x_Clone <- map(all_10x_Clone, `Idents<-`, value = "cloneType")

# Create an empty list to store DE results
FindAllMarkers_Clone_list_6 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_10x_Clone)){
  temp <- FindAllMarkers(all_10x_Clone[[i]])
  FindAllMarkers_Clone_list_6[[length(FindAllMarkers_Clone_list_6) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Clone_list_6) <- sample_names

# Annotate the results
FindAllMarkers_Clone_list_6 <- lapply(FindAllMarkers_Clone_list_6, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Clone_list_6, 
           "Output/DE/10x_atlases/Workup/6_clusters/FindAllMarkers_Clone_6_clusters.xlsx")

# 5.B.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Copy the all_10x list and switch the identity class to "Patient"
all_10x_Patient <- all_list
all_10x_Patient <- map(all_10x_Patient, `Idents<-`, value = "Patient")

# Create an empty list to store DE results
FindAllMarkers_Patient_list_6 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_10x_Patient)){
  temp <- FindAllMarkers(all_10x_Patient[[i]])
  FindAllMarkers_Patient_list_6[[length(FindAllMarkers_Patient_list_6) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Patient_list_6) <- sample_names

# Annotate the results
FindAllMarkers_Patient_list_6 <- lapply(FindAllMarkers_Patient_list_6, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_6, 
           "Output/DE/10x_atlases/Workup/6_clusters/FindAllMarkers_Patient_6_clusters.xlsx")


# 5.C FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 5.C.i FindMarkers - Itemise differences between Clusters -----------------------------------------

# Create an empty list to store results in
FindMarkers_All_list_6 <- list()

# Run for loop to iterate across all clusters
# NB: This will fail unless each input sample is 6 clusters
for (i in names(all_list)){
  for (j in 1:6) {
    temp <- FindMarkers(all_list[[i]], ident.1 = j-1)
    FindMarkers_All_list_6[[length(FindMarkers_All_list_6) +1]] <- temp %>% rownames_to_column %>% 
      rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list. NB: These need to be unique or export to *.xlsx will fail
names(FindMarkers_All_list_6) <- cluster_names_6

# Save out the DE results
write.xlsx(FindMarkers_All_list_6, 
           "Output/DE/10x_atlases/Workup/6_clusters/FindMarkers_All_6_clusters.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_All_list_6)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_All_list_6[[`i`]], 
                     lab = FindMarkers_All_list_6[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_All_list_6)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/6_clusters/",
                                   names(FindMarkers_All_list_6)[i], ".png"), scale = 2:1)
}

# 5.C.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create an empty list to store results in
FindMarkers_Tissue_list_6 <- list()

# Run for loop to iterate across all clusters
for (i in 1:6) {
  temp <- FindMarkers(all_list$all_10x, ident.1 = "BM", group.by = "Tissue", subset.ident = i-1)
  FindMarkers_Tissue_list_6[[i]] = temp %>% rownames_to_column %>% rename(gene = rowname) %>%
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>%
    unique()
}

# Rename the items in the list
names(FindMarkers_Tissue_list_6) <- cluster_names_6[1:6]

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list_6, 
           "Output/DE/10x_atlases/Workup/6_clusters/FindMarkers_Tissue_6_clusters.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Tissue_list_6)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Tissue_list_6[[`i`]], 
                     lab = FindMarkers_Tissue_list_6[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Tissue_list_6)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/6_clusters/FindMarkers_Tissue_",
                                   names(FindMarkers_Tissue_list_6)[i], ".png"), scale = 2:1)
}


# 5.D Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 5.D.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Tissue_list_6 <- list()

# Run for loop to iterate across all clusters
for (i in 1:6) {
  temp <- FindConservedMarkers(all_list$all_10x, ident.1 = i-1, grouping.var = "Tissue")
  FindConservedMarkers_Tissue_list_6[[i]] <- temp %>% 
    rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>% 
    unique()
}

# Rename the items in the list
names(FindConservedMarkers_Tissue_list_6) <- cluster_names_6[1:6]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_6, 
           "Output/DE/10x_atlases/Workup/6_clusters/FindConservedMarkers_Tissue_6_clusters.xlsx")


# 5.E Create and export figures to assist in cluster identification and annotation -----------------

# 5.E.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_6) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/10x_atlases/Workup/6_clusters/",
                  sample_names[i], ".png"))
}

# 5.E.ii Stacked barplots of cluster distribution --------------------------------------------------

# Calculate distribution of clusters and move results into a temp_list
temp_list <- 
  list("temp_all_10x" = as.data.frame(table(all_list$all_10x$seurat_clusters, 
                                            all_list$all_10x$orig.ident)),
       "temp_BM_10x" = as.data.frame(table(all_list$BM_10x$seurat_clusters, 
                                           all_list$BM_10x$orig.ident)),
       "temp_PB_10x" = as.data.frame(table(all_list$BM_10x$seurat_clusters, 
                                           all_list$BM_10x$orig.ident)))

# Output barplots of cluster distribution
for (i in 1:length(sample_names)) {
  ggplot(temp_list[[i]], 
         aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = rev(ccolours_6), name = "Cluster",
                      (rev(levels(all_list[[i]]@meta.data$seurat_clusters)))) + 
    labs(x = "Sample", y = "Cells", title = sample_names[i]) + theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Barplots/10x_atlases/Workup/6_clusters/",
                  sample_names[i], ".png"))
}


# 5.E.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all_10x object and output a Dotplot separated by tissue
top_10_markers_all_10x <- 
  FindAllMarkers_All_list_6$all_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all_10x <- unique(top_10_markers_all_10x$gene)
DotPlot(all_list$all_10x, 
        features = top_10_unique_all_10x, split.by = "Tissue", cols = c("blue", "red")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/6_clusters/",
                "BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the BM_10x object and output a Dotplot separated by patient
top_10_markers_BM_10x <- 
  FindAllMarkers_All_list_6$BM_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM_10x <- unique(top_10_markers_BM_10x$gene)
DotPlot(all_list$BM_10x, 
        features = top_10_unique_BM_10x, split.by = "Patient", cols = c("green", "purple")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "BM_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/6_clusters/",
                "BM pt43 vs. pt 63", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB_10x object and output a Dotplot separated by patient
top_10_markers_PB_10x <- 
  FindAllMarkers_All_list_6$PB_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB_10x <- unique(top_10_markers_PB_10x$gene)
DotPlot(all_list$PB_10x, 
        features = top_10_unique_PB_10x, split.by = "Patient", 
        cols = c("orange", "turquoise", "green", "purple")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "PB_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/6_clusters/",
                "PB all patients", ".png"), scale = 2:1)

# 6. Analysis of 5 cluster Seurat objects ----------------------------------------------------------

# 5.A Set resolutions of Seurat objects to allow clustering into 5 clusters
# Set resolution to 0.25 for all_10x object
all_list$all_10x$seurat_clusters <- 
  all_list$all_10x$integrated_snn_res.0.25
Idents(all_list$all_10x) <- "seurat_clusters"

# Set resolution to 0.25 for BM_10x object
all_list$BM_10x$seurat_clusters <- 
  all_list$BM_10x$integrated_snn_res.0.25
Idents(all_list$BM_10x) <- "seurat_clusters"

# Set resolution to 0.2 for PB_10x object
all_list$PB_10x$seurat_clusters <- 
  all_list$PB_10x$integrated_snn_res.0.2
Idents(all_list$PB_10x) <- "seurat_clusters"


# 6.B DE testing - FindAllMarkers ------------------------------------------------------------------

# 6.B.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

# Create an empty list to store DE results
FindAllMarkers_All_list_5 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_list)){
  temp <- FindAllMarkers(all_list[[i]])
  FindAllMarkers_All_list_5[[length(FindAllMarkers_All_list_5) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_All_list_5) <- sample_names

# Annotate the results 
FindAllMarkers_All_list_5 <- lapply(FindAllMarkers_All_list_5, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_All_list_5,
           "Output/DE/10x_atlases/Workup/5_clusters/FindAllMarkers_All_5_clusters.xlsx")

# 6.B.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 
# Copy the all_10x object and switch the identity class to "Tissue"
all_10x_Tissue <- all_list$all_10x
Idents(all_10x_Tissue) <- "Tissue"

# Run FindAllMarkers
FindAllMarkers_Tissue_list_5 <- FindAllMarkers(all_10x_Tissue)

# Annotate the results
FindAllMarkers_Tissue_list_5 <- inner_join(x = FindAllMarkers_Tissue_list_5, 
                                           y = annotations[, c("gene_name", "description")],
                                           by = c("gene" = "gene_name")) %>% unique()

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = BM, pct.2 = PB
FindAllMarkers_Tissue_list_5 <- 
  FindAllMarkers_Tissue_list_5[!(FindAllMarkers_Tissue_list_5$cluster == "PB"),]

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list_5,
           "Output/DE/10x_atlases/Workup/5_clusters/FindAllMarkers_Tissue_5_clusters.xlsx")

# Output Volcano plots for Tissue restricted differences
EnhancedVolcano(FindAllMarkers_Tissue_list_5, 
                lab = FindAllMarkers_Tissue_list_5$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "BM vs. PB",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 12, pointSize = 2, labSize = 4, legendLabSize = 12,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/5_clusters/FindAllMarkers_Tissue_",
                "BM vs. PB", ".png"), scale = 2:1)

# 6.B.iii DE testing - FindAllMarkers_clone - Assess global differences between clonal bins --------

# Copy the all_10x object and switch the identity class to "cloneType"
all_10x_Clone <- all_list
all_10x_Clone <- map(all_10x_Clone, `Idents<-`, value = "cloneType")

# Create an empty list to store DE results
FindAllMarkers_Clone_list_5 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_10x_Clone)){
  temp <- FindAllMarkers(all_10x_Clone[[i]])
  FindAllMarkers_Clone_list_5[[length(FindAllMarkers_Clone_list_5) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Clone_list_5) <- sample_names

# Annotate the results
FindAllMarkers_Clone_list_5 <- lapply(FindAllMarkers_Clone_list_5, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Clone_list_5, 
           "Output/DE/10x_atlases/Workup/5_clusters/FindAllMarkers_Clone_5_clusters.xlsx")

# 6.B.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Copy the all_10x list and switch the identity class to "Patient"
all_10x_Patient <- all_list
all_10x_Patient <- map(all_10x_Patient, `Idents<-`, value = "Patient")

# Create an empty list to store DE results
FindAllMarkers_Patient_list_5 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_10x_Patient)){
  temp <- FindAllMarkers(all_10x_Patient[[i]])
  FindAllMarkers_Patient_list_5[[length(FindAllMarkers_Patient_list_5) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Patient_list_5) <- sample_names

# Annotate the results
FindAllMarkers_Patient_list_5 <- lapply(FindAllMarkers_Patient_list_5, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_5, 
           "Output/DE/10x_atlases/Workup/5_clusters/FindAllMarkers_Patient_5_clusters.xlsx")


# 6.C FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 6.C.i FindMarkers - Itemise differences between Clusters -----------------------------------------

# Create an empty list to store results in
FindMarkers_All_list_5 <- list()

# Run for loop to iterate across all clusters
# NB: This will fail unless each input sample is 5 clusters
for (i in names(all_list)){
  for (j in 1:5) {
    temp <- FindMarkers(all_list[[i]], ident.1 = j-1)
    FindMarkers_All_list_5[[length(FindMarkers_All_list_5) +1]] <- temp %>% rownames_to_column %>% 
      rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list. NB: These need to be unique or export to *.xlsx will fail
names(FindMarkers_All_list_5) <- cluster_names_5

# Save out the DE results
write.xlsx(FindMarkers_All_list_5, 
           "Output/DE/10x_atlases/Workup/5_clusters/FindMarkers_All_5_clusters.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_All_list_5)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_All_list_5[[`i`]], 
                     lab = FindMarkers_All_list_5[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_All_list_5)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/5_clusters/",
                                   names(FindMarkers_All_list_5)[i], ".png"), scale = 2:1)
}

# 6.C.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create an empty list to store results in
FindMarkers_Tissue_list_5 <- list()

# Run for loop to iterate across all clusters
for (i in 1:5) {
  temp <- FindMarkers(all_list$all_10x, ident.1 = "BM", group.by = "Tissue", subset.ident = i-1)
  FindMarkers_Tissue_list_5[[i]] = temp %>% rownames_to_column %>% rename(gene = rowname) %>%
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>%
    unique()
}

# Rename the items in the list
names(FindMarkers_Tissue_list_5) <- cluster_names_5[1:5]

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list_5, 
           "Output/DE/10x_atlases/Workup/5_clusters/FindMarkers_Tissue_5_clusters.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Tissue_list_5)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Tissue_list_5[[`i`]], 
                     lab = FindMarkers_Tissue_list_5[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Tissue_list_5)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x_atlases/Workup/5_clusters/FindMarkers_Tissue_",
                                   names(FindMarkers_Tissue_list_5)[i], ".png"), scale = 2:1)
}


# 6.D Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 6.D.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Tissue_list_5 <- list()

# Run for loop to iterate across all clusters
for (i in 1:5) {
  temp <- FindConservedMarkers(all_list$all_10x, ident.1 = i-1, grouping.var = "Tissue")
  FindConservedMarkers_Tissue_list_5[[i]] <- temp %>% 
    rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>% 
    unique()
}

# Rename the items in the list
names(FindConservedMarkers_Tissue_list_5) <- cluster_names_5[1:5]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_5, 
           "Output/DE/10x_atlases/Workup/5_clusters/FindConservedMarkers_Tissue_5_clusters.xlsx")


# 6.E Create and export figures to assist in cluster identification and annotation -----------------

# 6.E.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_5) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/10x_atlases/Workup/5_clusters/",
                  sample_names[i], ".png"))
}

# 6.E.ii Stacked barplots of cluster distribution --------------------------------------------------

# Calculate distribution of clusters and move results into a temp_list
temp_list <- 
  list("temp_all_10x" = as.data.frame(table(all_list$all_10x$seurat_clusters, 
                                            all_list$all_10x$orig.ident)),
       "temp_BM_10x" = as.data.frame(table(all_list$BM_10x$seurat_clusters, 
                                           all_list$BM_10x$orig.ident)),
       "temp_PB_10x" = as.data.frame(table(all_list$BM_10x$seurat_clusters, 
                                           all_list$BM_10x$orig.ident)))

# Output barplots of cluster distribution
for (i in 1:length(sample_names)) {
  ggplot(temp_list[[i]], 
         aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = rev(ccolours_5), name = "Cluster",
                      (rev(levels(all_list[[i]]@meta.data$seurat_clusters)))) + 
    labs(x = "Sample", y = "Cells", title = sample_names[i]) + theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Barplots/10x_atlases/Workup/5_clusters/",
                  sample_names[i], ".png"))
}


# 6.E.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all_10x object and output a Dotplot separated by tissue
top_10_markers_all_10x <- 
  FindAllMarkers_All_list_5$all_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all_10x <- unique(top_10_markers_all_10x$gene)
DotPlot(all_list$all_10x, 
        features = top_10_unique_all_10x, split.by = "Tissue", cols = c("blue", "red")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/5_clusters/",
                "BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the BM_10x object and output a Dotplot separated by patient
top_10_markers_BM_10x <- 
  FindAllMarkers_All_list_5$BM_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM_10x <- unique(top_10_markers_BM_10x$gene)
DotPlot(all_list$BM_10x, 
        features = top_10_unique_BM_10x, split.by = "Patient", cols = c("green", "purple")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "BM_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/5_clusters/",
                "BM pt43 vs. pt 53", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB_10x object and output a Dotplot separated by patient
top_10_markers_PB_10x <- 
  FindAllMarkers_All_list_5$PB_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB_10x <- unique(top_10_markers_PB_10x$gene)
DotPlot(all_list$PB_10x, 
        features = top_10_unique_PB_10x, split.by = "Patient", 
        cols = c("orange", "turquoise", "green", "purple")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "PB_10x") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/10x_atlases/Workup/5_clusters/",
                "PB all patients", ".png"), scale = 2:1)


# 7. Save out data and print stats -----------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean up the environment
rm(all_10x_Clone, all_10x_Patient, all_10x_Tissue, all_list, all_refs, counts, temp, temp_list, 
   top_10_markers_all_10x, top_10_markers_BM_10x, top_10_markers_PB_10x, black_list, ccolours_5,
   ccolours_6, ccolours_7, cluster_names_5, cluster_names_6, cluster_names_7, filename, i, j, 
   top_10_unique_all_10x, top_10_unique_BM_10x, top_10_unique_PB_10x)

# Save out the DE_results for further analysis
save.image("Data/R_out/ProjecTILs/DE_Workup_10x_atlases.rds")

gc()

print("#># Finished running '4. ProjecTILs/10x Reference atlases DE workup' script")