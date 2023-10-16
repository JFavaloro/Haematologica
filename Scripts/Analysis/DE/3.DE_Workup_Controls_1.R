# Script information -------------------------------------------------------------------------------

# Title: Controls dataset #1 - DE Workup
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will repeat the DE analysis we performed on our 10x dataset on 
# Controls dataset #1 from Szabo et al., 2019: https://doi.org/10.1038/s41467-019-12464-3.
# Where insufficient data exists to allow looping through the analysis without erroring, we will 
# comment these out, or modify the testing and note which are failing or too sparse for analysis.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '3. Controls dataset #1 - DE Workup' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, openxlsx, multtest, metap, EnhancedVolcano)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed Controls #1 dataset
integrated_list <- readRDS("Data/R_out/Controls_1/integrated_list.rds")

# Load additional resources
black_list <- readRDS("Data/R_out/black_list.rds")
annotations <- read.csv("Data/Public_data/annotations.csv") 

# Modify black_list to retain TCR genes
black_list_TCR <- 
  black_list[!grepl("TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC", black_list)]

# Define function to remove black listed genes from discriminant genes (excluding TCR genes)
keep_TCR <- function(x){
  x <- x[!((x$gene) %in% black_list_TCR),]
  return(x)
}

# Define function to remove redundant Tissue entries from discriminant genes
remove_PB <- function(x){
  x <- x[!(x$cluster == "PB"),]
  return(x)
}

# Define function to remove redundant Stimulation status entries from discriminant genes
remove_rest <- function(x){
  x <- x[!(x$cluster == "rest"),]
  return(x)
}

# Define function to remove redundant Patient status entries from discriminant genes
remove_Donor <- function(x){
  x <- x[!(x$cluster == "Donor2"),]
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

# Create a character vector of samples names
sample_names <- c("BM.integrated", "PB.integrated", "act.integrated", 
                  "rest.integrated", "all.integrated")

# Define cluster names for naming of lists - 7 clusters
cluster_names_7 <- apply(expand.grid(c("cluster_"), c(0:6)), 1, paste, collapse="") %>% 
  rep(3) %>% paste0(c(rep("_BM.integrated", 7), rep("_PB.integrated", 7), rep("_act.integrated", 7),
                  rep("_rest.integrated", 7), rep("_all.integrated", 7)))

# Define cluster names for naming of lists - 6 clusters
cluster_names_6 <- apply(expand.grid(c("cluster_"), c(0:5)), 1, paste, collapse="") %>% 
  rep(3) %>% paste0(c(rep("_BM.integrated", 6), rep("_PB.integrated", 6), rep("_act.integrated", 6),
                      rep("_rest.integrated", 6), rep("_all.integrated", 6)))

# Define cluster names for naming of lists - 5 clusters
cluster_names_5 <- apply(expand.grid(c("cluster_"), c(0:4)), 1, paste, collapse="") %>% 
  rep(3) %>% paste0(c(rep("_BM.integrated", 5), rep("_PB.integrated", 5), rep("_act.integrated", 5),
                      rep("_rest.integrated", 5), rep("_all.integrated", 5)))

# Create a character vector of colours for plotting clusters on UMAP
# NB: TEM = red ("#F8766D"), TTE = yellow ("#B79F00"), TN = green ("#00BA38"),
# Cyto_TEM = turquoise ("#00BFC4"), PRE_EX = blue ("#619CFF"), TCM = pink ("#F564E3")
# These need to be finalised
scales::show_col(scales::hue_pal()(6))
scales::show_col(scales::hue_pal()(12))
scales::show_col(scales::hue_pal()(25))
ccolours_7 <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3", "#FF0000")
ccolours_6 <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
ccolours_5 <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF")


# 3. Pre-processing --------------------------------------------------------------------------------

# Copy the integrated_list to the all_list to allow copying of previous code and backup
all_list <- integrated_list

# Switch default assay to "RNA" for differential expression testing
all_list <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Remove non helpful genes from DE results. #NB: This will remove genes from the object
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
           "Output/DE/Controls_1/Workup/7_clusters/FindAllMarkers_All_7_clusters.xlsx")

# 4.A.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 
# Copy the all integrated objects into a new list and switch the identity class to "Tissue"
all_Tissue <- all_list
all_Tissue <- map(all_Tissue, `Idents<-`, value = "Tissue")

# Create an empty list to store DE results
FindAllMarkers_Tissue_list_7 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Tissue)){
  temp <- FindAllMarkers(all_Tissue[[i]])
  FindAllMarkers_Tissue_list_7[[length(FindAllMarkers_Tissue_list_7) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Tissue_list_7) <- sample_names

# Remove BM.integrated and PB.integrated objects
FindAllMarkers_Tissue_list_7 <- FindAllMarkers_Tissue_list_7[3:5]

# Annotate the results 
FindAllMarkers_Tissue_list_7 <- lapply(FindAllMarkers_Tissue_list_7, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = BM, pct.2 = PB
FindAllMarkers_Tissue_list_7 <- lapply(FindAllMarkers_Tissue_list_7, remove_PB)

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list_7,
           "Output/DE/Controls_1/Workup/7_clusters/FindAllMarkers_Tissue_7_clusters.xlsx")

# Output Volcanoplots for Tissue restricted differences
for (i in 1:length(FindAllMarkers_Tissue_list_7)) {
  temp <- 
    (EnhancedVolcano(FindAllMarkers_Tissue_list_7[[`i`]], 
                     lab = FindAllMarkers_Tissue_list_7[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindAllMarkers_Tissue_list_7)[i],
                     subtitle = "BM vs. PB Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/7_clusters/FindAllMarkers_Tissue_",
                                   names(FindAllMarkers_Tissue_list_7)[i], ".png"), scale = 2:1)
}

# 4.A.iii DE testing - FindAllMarkers_Stim - Assess global differences between Stimulation status -

# Copy the all_list into a new list and switch the identity class to "stimulation_status"
all_Stim <- all_list
all_Stim <- map(all_Stim, `Idents<-`, value = "stimulation_status")

# Create an empty list to store DE results
FindAllMarkers_Stim_list_7 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Stim)){
  temp <- FindAllMarkers(all_Stim[[i]])
  FindAllMarkers_Stim_list_7[[length(FindAllMarkers_Stim_list_7) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Stim_list_7) <- sample_names

# Remove act.integrated and rest.integrated objects
FindAllMarkers_Stim_list_7 <- FindAllMarkers_Stim_list_7[-c(3:4)]

# Annotate the results
FindAllMarkers_Stim_list_7 <- lapply(FindAllMarkers_Stim_list_7, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. rest results) pct.1 = act, pct.2 = rest
FindAllMarkers_Stim_list_7 <- lapply(FindAllMarkers_Stim_list_7, remove_rest)

# Save out the DE results
write.xlsx(FindAllMarkers_Stim_list_7, 
           "Output/DE/Controls_1/Workup/7_clusters/FindAllMarkers_Stim_7_clusters.xlsx")

# 4.A.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Copy the all_list into a new list and switch the identity class to "orig.ident"
all_Patient <- all_list
all_Patient <- map(all_Patient, `Idents<-`, value = "orig.ident")

# Create an empty list to store DE results
FindAllMarkers_Patient_list_7 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Patient)){
  temp <- FindAllMarkers(all_Patient[[i]])
  FindAllMarkers_Patient_list_7[[length(FindAllMarkers_Patient_list_7) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Patient_list_7) <- sample_names

# Annotate the results
FindAllMarkers_Patient_list_7 <- lapply(FindAllMarkers_Patient_list_7, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = Donor1, pct.2 = Donor2
FindAllMarkers_Patient_list_7 <- lapply(FindAllMarkers_Patient_list_7, remove_Donor)


# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_7, 
           "Output/DE/Controls_1/Workup/7_clusters/FindAllMarkers_Patient_7_clusters.xlsx")


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
           "Output/DE/Controls_1/Workup/7_clusters/FindMarkers_All_7_clusters.xlsx")

# Output Volcanoplots for all tests
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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/7_clusters/",
                                    names(FindMarkers_All_list_7)[i], ".png"), scale = 2:1)
}

# 4.B.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create an empty list to store results in
FindMarkers_Tissue_list_7 <- list()

# Copy the all_list object and remove the BM.integrated and PB.integrated objects
all_list_temp <- all_list[-c(1,2)]

# Run for loop to iterate across all clusters 
for (i in names(all_list_temp)){
  for (j in 1:7) {
    temp <- FindMarkers(all_list_temp[[`i`]], 
                        ident.1 = "BM", group.by = "Tissue", subset.ident = j-1)
    FindMarkers_Tissue_list_7[[length(FindMarkers_Tissue_list_7) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Tissue_list_7) <- cluster_names_7[15:35]

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list_7, 
           "Output/DE/Controls_1/Workup/7_clusters/FindMarkers_Tissue_7_clusters.xlsx")

# Output Volcanoplots for all tests
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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/7_clusters/FindMarkers_Tissue_",
                                   names(FindMarkers_Tissue_list_7)[i], ".png"), scale = 2:1)
}


# 4.B.iii FindMarkers - Itemise differences between activation and resting -------------------------

# Create an empty list to store results in
FindMarkers_Stim_list_7 <- list()

# Copy the all_list object and remove the act.integrated and rest.integrated objects
all_list_temp <- all_list[-c(3,4)]

# Run for loop to iterate across all clusters
# Needs to be run with argument min.cells.group = 1 - data is too sparse
for (i in names(all_list_temp)){
  for (j in 1:7) {
    temp <- FindMarkers(all_list_temp[[`i`]], 
                        ident.1 = "act", group.by = "stimulation_status", subset.ident = j-1, min.cells.group = 1)
    FindMarkers_Stim_list_7[[length(FindMarkers_Stim_list_7) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Stim_list_7) <- cluster_names_7[c(1:14, 29:35)]

# Save out the DE results
write.xlsx(FindMarkers_Stim_list_7, 
           "Output/DE/Controls_1/Workup/7_clusters/FindMarkers_Stim_7_clusters.xlsx")

# Output Volcanoplots for all tests
for (i in 1:length(FindMarkers_Stim_list_7)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Stim_list_7[[`i`]], 
                     lab = FindMarkers_Stim_list_7[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Stim_list_7)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/7_clusters/Stim_",
                                   names(FindMarkers_Stim_list_7)[i], ".png"), scale = 2:1)
}


# 4.C Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 4.C.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Tissue_list_7 <- list()

# Copy the all_list object and remove the BM.integrated and PB.integrated objects
all_list_temp <- all_list[-c(1,2)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:7) {
    temp <- FindConservedMarkers(all_list_temp[[`i`]], 
                                 ident.1 = j-1, grouping.var = "Tissue")
    FindConservedMarkers_Tissue_list_7[[length(FindConservedMarkers_Tissue_list_7) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindConservedMarkers_Tissue_list_7) <- cluster_names_7[15:35]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_7, 
           "Output/DE/Controls_1/Workup/7_clusters/FindConservedMarkers_Tissue_7_clusters.xlsx")

# 4.C.ii FindConservedMarkers - Averages across Stimulation status ---------------------------------

# Create an empty list to store results in
FindConservedMarkers_Stim_list_7 <- list()

# Copy the all_list object and remove the act.integrated and rest.integrated objects
all_list_temp <- all_list[-c(3,4)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:7) {
    temp <- FindConservedMarkers(all_list_temp[[`i`]], 
                                 ident.1 = j-1, grouping.var = "stimulation_status")
    FindConservedMarkers_Stim_list_7[[length(FindConservedMarkers_Stim_list_7) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindConservedMarkers_Stim_list_7) <- cluster_names_7[c(1:14, 29:35)]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_7, 
           "Output/DE/Controls_1/Workup/7_clusters/FindConservedMarkers_Stim_7_clusters.xlsx")


# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_7) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/Controls_1/Workup/7_clusters/",
                  sample_names[i], ".png"))
}

# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Calculate distribution of clusters and move results into a temp_list
temp_list <- 
  list("temp_BM.integrated" = as.data.frame(table(all_list$BM.integrated$seurat_clusters, 
                                                  all_list$BM.integrated$orig.ident)),
       "temp_PB.integrated" = as.data.frame(table(all_list$PB.integrated$seurat_clusters, 
                                                  all_list$PB.integrated$orig.ident)), 
       "temp_act.integrated" = as.data.frame(table(all_list$act.integrated$seurat_clusters,
                                                   all_list$act.integrated$orig.ident)),
       "temp_rest.integrated" = as.data.frame(table(all_list$rest.integrated$seurat_clusters,
                                                    all_list$rest.integrated$orig.ident)), 
       "temp_all.integrated" = as.data.frame(table(all_list$all.integrated$seurat_clusters,
                                                   all_list$all.integrated$orig.ident)))

# Output barplots of cluster distribution
for (i in 1:length(all_list)) {
  ggplot(temp_list[[i]], 
         aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = rev(ccolours_7), name = "Cluster",
                      (rev(levels(all_list[[i]]@meta.data$seurat_clusters)))) + 
    labs(x = "Sample", y = "Cells", title = sample_names[i]) + theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Barplots/Controls_1/Workup/7_clusters/",
                  sample_names[i], ".png"))
}

# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the act.integrated object and output a Dotplot separated by tissue
top_10_markers_act <- 
  FindAllMarkers_All_list_7$act.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_act <- unique(top_10_markers_act$gene)
DotPlot(all_list$act.integrated, 
        features = top_10_unique_act, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "act") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/7_clusters/",
                "Activated BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the rest.integrated object and output a Dotplot separated by tissue
top_10_markers_rest <- 
  FindAllMarkers_All_list_7$rest.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_rest <- unique(top_10_markers_rest$gene)
DotPlot(all_list$rest.integrated, 
        features = top_10_unique_rest, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "rest") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/7_clusters/",
                "Resting BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by tissue
top_10_markers_all <- 
  FindAllMarkers_All_list_7$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/7_clusters/",
                "All BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the BM.integrated object and output a Dotplot separated by activation status
top_10_markers_BM <- 
  FindAllMarkers_All_list_7$BM.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM <- unique(top_10_markers_BM$gene)
DotPlot(all_list$BM.integrated, 
        features = top_10_unique_BM, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "BM") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/7_clusters/",
                "BM act vs. rest", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB.integrated object and output a Dotplot separated by activation status
top_10_markers_PB <- 
  FindAllMarkers_All_list_7$PB.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB <- unique(top_10_markers_PB$gene)
DotPlot(all_list$PB.integrated, 
        features = top_10_unique_PB, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "PB") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/7_clusters/",
                "PB act vs. rest", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by activation status
top_10_markers_all <- 
  FindAllMarkers_All_list_7$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/7_clusters/",
                "All act vs. rest", ".png"), scale = 2:1)


# 5. Analysis of 6 cluster Seurat objects ----------------------------------------------------------

# 5.A Set resolutions of Seurat objects to allow clustering into 6 clusters
# Set resolution to 0.6 for BM.integrated object
all_list$BM.integrated$seurat_clusters <- 
  all_list$BM.integrated$integrated_snn_res.0.6
Idents(all_list$BM.integrated) <- "seurat_clusters"

# Set resolution to 0.4 for PB.integrated object
all_list$PB.integrated$seurat_clusters <- 
  all_list$PB.integrated$integrated_snn_res.0.4
Idents(all_list$PB.integrated) <- "seurat_clusters"

# Set resolution to 0.4 for act.integrated object
all_list$act.integrated$seurat_clusters <- 
  all_list$act.integrated$integrated_snn_res.0.4
Idents(all_list$act.integrated) <- "seurat_clusters"

# Set resolution to 0.44 for rest.integrated object
all_list$rest.integrated$seurat_clusters <- 
  all_list$rest.integrated$integrated_snn_res.0.44
Idents(all_list$rest.integrated) <- "seurat_clusters"

# Set resolution to 0.4 for all.integrated object
all_list$all.integrated$seurat_clusters <- 
  all_list$all.integrated$integrated_snn_res.0.4
Idents(all_list$all.integrated) <- "seurat_clusters"

# 5.B DE testing - FindAllMarkers

# 5.A.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

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
           "Output/DE/Controls_1/Workup/6_clusters/FindAllMarkers_All_6_clusters.xlsx")

# 5.A.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 
# Copy the all integrated objects into a new list and switch the identity class to "Tissue"
all_Tissue <- all_list
all_Tissue <- map(all_Tissue, `Idents<-`, value = "Tissue")

# Create an empty list to store DE results
FindAllMarkers_Tissue_list_6 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Tissue)){
  temp <- FindAllMarkers(all_Tissue[[i]])
  FindAllMarkers_Tissue_list_6[[length(FindAllMarkers_Tissue_list_6) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Tissue_list_6) <- sample_names

# Remove BM.integrated and PB.integrated objects
FindAllMarkers_Tissue_list_6 <- FindAllMarkers_Tissue_list_6[3:5]

# Annotate the results 
FindAllMarkers_Tissue_list_6 <- lapply(FindAllMarkers_Tissue_list_6, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = BM, pct.2 = PB
FindAllMarkers_Tissue_list_6 <- lapply(FindAllMarkers_Tissue_list_6, remove_PB)

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list_6,
           "Output/DE/Controls_1/Workup/6_clusters/FindAllMarkers_Tissue_6_clusters.xlsx")

# Output Volcanoplots for Tissue restricted differences
for (i in 1:length(FindAllMarkers_Tissue_list_6)) {
  temp <- 
    (EnhancedVolcano(FindAllMarkers_Tissue_list_6[[`i`]], 
                     lab = FindAllMarkers_Tissue_list_6[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindAllMarkers_Tissue_list_6)[i],
                     subtitle = "BM vs. PB Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/6_clusters/FindAllMarkers_Tissue_",
                                   names(FindAllMarkers_Tissue_list_6)[i], ".png"), scale = 2:1)
}

# 5.A.iii DE testing - FindAllMarkers_Stim - Assess global differences between Stimulation status -

# Copy the all_list into a new list and switch the identity class to "stimulation_status"
all_Stim <- all_list
all_Stim <- map(all_Stim, `Idents<-`, value = "stimulation_status")

# Create an empty list to store DE results
FindAllMarkers_Stim_list_6 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Stim)){
  temp <- FindAllMarkers(all_Stim[[i]])
  FindAllMarkers_Stim_list_6[[length(FindAllMarkers_Stim_list_6) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Stim_list_6) <- sample_names

# Remove act.integrated and rest.integrated objects
FindAllMarkers_Stim_list_6 <- FindAllMarkers_Stim_list_6[-c(3:4)]

# Annotate the results
FindAllMarkers_Stim_list_6 <- lapply(FindAllMarkers_Stim_list_6, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. rest results) pct.1 = act, pct.2 = rest
FindAllMarkers_Stim_list_6 <- lapply(FindAllMarkers_Stim_list_6, remove_rest)

# Save out the DE results
write.xlsx(FindAllMarkers_Stim_list_6, 
           "Output/DE/Controls_1/Workup/6_clusters/FindAllMarkers_Stim_6_clusters.xlsx")

# 5.A.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Copy the all_list into a new list and switch the identity class to "orig.ident"
all_Patient <- all_list
all_Patient <- map(all_Patient, `Idents<-`, value = "orig.ident")

# Create an empty list to store DE results
FindAllMarkers_Patient_list_6 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Patient)){
  temp <- FindAllMarkers(all_Patient[[i]])
  FindAllMarkers_Patient_list_6[[length(FindAllMarkers_Patient_list_6) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Patient_list_6) <- sample_names

# Annotate the results
FindAllMarkers_Patient_list_6 <- lapply(FindAllMarkers_Patient_list_6, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = Donor1, pct.2 = Donor2
FindAllMarkers_Patient_list_6 <- lapply(FindAllMarkers_Patient_list_6, remove_Donor)


# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_6, 
           "Output/DE/Controls_1/Workup/6_clusters/FindAllMarkers_Patient_6_clusters.xlsx")


# 5.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 5.B.i FindMarkers - Itemise differences between Clusters -----------------------------------------

# Create an empty list to store results in
FindMarkers_All_list_6 <- list()

# Run for loop to iterate across all clusters
# NB: This will fail unless each input sample is 6 clusters
for (i in names(all_list)){
  for (j in 1:6) {
    temp <- FindMarkers(all_list[[i]], ident.1 = j-1)
    FindMarkers_All_list_6[[length(FindMarkers_All_list_6) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list. NB: These need to be unique or export to *.xlsx will fail
names(FindMarkers_All_list_6) <- cluster_names_6

# Save out the DE results
write.xlsx(FindMarkers_All_list_6, 
           "Output/DE/Controls_1/Workup/6_clusters/FindMarkers_All_6_clusters.xlsx")

# Output Volcanoplots for all tests
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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/6_clusters/",
                                    names(FindMarkers_All_list_6)[i], ".png"), scale = 2:1)
}

# 5.B.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create an empty list to store results in
FindMarkers_Tissue_list_6 <- list()

# Copy the all_list object and remove the BM.integrated and PB.integrated objects
# NB: We will exclude the 'rest' object as there are no BM cells in cluster 6 for resolutions 
# 0.42-0.48. This can be run manually if required, but the loop will fail if included
all_list_temp = all_list[-c(1,2,4)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:6) {
    temp <- FindMarkers(all_list_temp[[`i`]], 
                        ident.1 = "BM", group.by = "Tissue", subset.ident = j-1)
    FindMarkers_Tissue_list_6[[length(FindMarkers_Tissue_list_6) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Tissue_list_6) <- cluster_names_6[c(13:18,25:30)]

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list_6, 
           "Output/DE/Controls_1/Workup/6_clusters/FindMarkers_Tissue_6_clusters.xlsx")

# Output Volcanoplots for all tests
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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/6_clusters/FindMarkers_Tissue_",
                                   names(FindMarkers_Tissue_list_6)[i], ".png"), scale = 2:1)
}


# 5.B.iii FindMarkers - Itemise differences between activation and resting -------------------------

# Create an empty list to store results in
FindMarkers_Stim_list_6 <- list()

# Copy the all_list object and remove the act.integrated and rest.integrated objects
all_list_temp <- all_list[-c(3,4)]

# Run for loop to iterate across all clusters
# Needs to be run with argument min.cells.group = 1 - data is too sparse
for (i in names(all_list_temp)){
  for (j in 1:6) {
    temp <- FindMarkers(all_list_temp[[`i`]], 
                        ident.1 = "act", group.by = "stimulation_status", subset.ident = j-1, min.cells.group = 1)
    FindMarkers_Stim_list_6[[length(FindMarkers_Stim_list_6) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Stim_list_6) <- cluster_names_6[c(1:12, 25:30)]

# Save out the DE results
write.xlsx(FindMarkers_Stim_list_6, 
           "Output/DE/Controls_1/Workup/6_clusters/FindMarkers_Stim_6_clusters.xlsx")

# Output Volcanoplots for all tests
for (i in 1:length(FindMarkers_Stim_list_6)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Stim_list_6[[`i`]], 
                     lab = FindMarkers_Stim_list_6[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Stim_list_6)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/6_clusters/Stim_",
                                   names(FindMarkers_Stim_list_6)[i], ".png"), scale = 2:1)
}


# 5.C Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 5.C.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Tissue_list_6 <- list()

# Copy the all_list object and remove the BM.integrated and PB.integrated objects
all_list_temp <- all_list[-c(1,2)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:6) {
    temp <- FindConservedMarkers(all_list_temp[[`i`]], 
                                 ident.1 = j-1, grouping.var = "Tissue")
    FindConservedMarkers_Tissue_list_6[[length(FindConservedMarkers_Tissue_list_6) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindConservedMarkers_Tissue_list_6) <- cluster_names_6[13:30]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_6, 
           "Output/DE/Controls_1/Workup/6_clusters/FindConservedMarkers_Tissue_6_clusters.xlsx")

# 5.C.ii FindConservedMarkers - Averages across Stimulation status ---------------------------------

# Create an empty list to store results in
FindConservedMarkers_Stim_list_6 = list()

# Copy the all_list object and remove the act.integrated and rest.integrated objects
all_list_temp <- all_list[-c(3,4)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:6) {
    temp <- FindConservedMarkers(all_list_temp[[`i`]], 
                                 ident.1 = j-1, grouping.var = "stimulation_status")
    FindConservedMarkers_Stim_list_6[[length(FindConservedMarkers_Stim_list_6) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindConservedMarkers_Stim_list_6) <- cluster_names_6[c(1:12, 25:30)]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_6, 
           "Output/DE/Controls_1/Workup/6_clusters/FindConservedMarkers_Stim_6_clusters.xlsx")


# 5.D Create and export figures to assist in cluster identification and annotation -----------------

# 5.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_6) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/Controls_1/Workup/6_clusters/",
                  sample_names[i], ".png"))
}

# 5.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Calculate distribution of clusters and move results into a temp_list
temp_list <- 
  list("temp_BM.integrated" = as.data.frame(table(all_list$BM.integrated$seurat_clusters, 
                                                  all_list$BM.integrated$orig.ident)),
       "temp_PB.integrated" = as.data.frame(table(all_list$PB.integrated$seurat_clusters, 
                                                  all_list$PB.integrated$orig.ident)), 
       "temp_act.integrated" = as.data.frame(table(all_list$act.integrated$seurat_clusters,
                                                   all_list$act.integrated$orig.ident)),
       "temp_rest.integrated" = as.data.frame(table(all_list$rest.integrated$seurat_clusters,
                                                    all_list$rest.integrated$orig.ident)), 
       "temp_all.integrated" = as.data.frame(table(all_list$all.integrated$seurat_clusters,
                                                   all_list$all.integrated$orig.ident)))

# Output barplots of cluster distribution
for (i in 1:length(all_list)) {
  ggplot(temp_list[[i]], 
         aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = rev(ccolours_6), name = "Cluster",
                      (rev(levels(all_list[[i]]@meta.data$seurat_clusters)))) + 
    labs(x = "Sample", y = "Cells", title = sample_names[i]) + theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Barplots/Controls_1/Workup/6_clusters/",
                  sample_names[i], ".png"))
}

# 5.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the act.integrated object and output a Dotplot separated by tissue
top_10_markers_act <- 
  FindAllMarkers_All_list_6$act.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_act <- unique(top_10_markers_act$gene)
DotPlot(all_list$act.integrated, 
        features = top_10_unique_act, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "act") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/6_clusters/",
                "Activated BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the rest.integrated object and output a Dotplot separated by tissue
top_10_markers_rest <- 
  FindAllMarkers_All_list_6$rest.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_rest <- unique(top_10_markers_rest$gene)
DotPlot(all_list$rest.integrated, 
        features = top_10_unique_rest, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "rest") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/6_clusters/",
                "Resting BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by tissue
top_10_markers_all <- 
  FindAllMarkers_All_list_6$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/6_clusters/",
                "All BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the BM.integrated object and output a Dotplot separated by activation status
top_10_markers_BM <- 
  FindAllMarkers_All_list_6$BM.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM <- unique(top_10_markers_BM$gene)
DotPlot(all_list$BM.integrated, 
        features = top_10_unique_BM, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "BM") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/6_clusters/",
                "BM act vs. rest", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB.integrated object and output a Dotplot separated by activation status
top_10_markers_PB <- 
  FindAllMarkers_All_list_6$PB.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB <- unique(top_10_markers_PB$gene)
DotPlot(all_list$PB.integrated, 
        features = top_10_unique_PB, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "PB") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/6_clusters/",
                "PB act vs. rest", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by activation status
top_10_markers_all <- 
  FindAllMarkers_All_list_6$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/6_clusters/",
                "All act vs. rest", ".png"), scale = 2:1)


# 6. Analysis of 5 cluster Seurat objects ----------------------------------------------------------

# 6.A Set resolutions of Seurat objects to allow clustering into 5 clusters
# Set resolution to 0.6 for BM.integrated object
all_list$BM.integrated$seurat_clusters <- 
  all_list$BM.integrated$integrated_snn_res.0.6
Idents(all_list$BM.integrated) <- "seurat_clusters"

# Set resolution to 0.4 for PB.integrated object
all_list$PB.integrated$seurat_clusters <- 
  all_list$PB.integrated$integrated_snn_res.0.4
Idents(all_list$PB.integrated) <- "seurat_clusters"

# Set resolution to 0.3 for act.integrated object
all_list$act.integrated$seurat_clusters <- 
  all_list$act.integrated$integrated_snn_res.0.3
Idents(all_list$act.integrated) <- "seurat_clusters"

# Set resolution to 0.3 for rest.integrated object
all_list$rest.integrated$seurat_clusters <- 
  all_list$rest.integrated$integrated_snn_res.0.3
Idents(all_list$rest.integrated) <- "seurat_clusters"

# Set resolution to 0.4 for all.integrated object
all_list$all.integrated$seurat_clusters <- 
  all_list$all.integrated$integrated_snn_res.0.4
Idents(all_list$all.integrated) <- "seurat_clusters"

# 6.B DE testing - FindAllMarkers

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
           "Output/DE/Controls_1/Workup/5_clusters/FindAllMarkers_All_5_clusters.xlsx")

# 6.B.ii DE testing - FindAllMarkers_Tissue - Assess global differences between BM and PB ---------- 
# Copy the all integrated objects into a new list and switch the identity class to "Tissue"
all_Tissue <- all_list
all_Tissue <- map(all_Tissue, `Idents<-`, value = "Tissue")

# Create an empty list to store DE results
FindAllMarkers_Tissue_list_5 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Tissue)){
  temp <- FindAllMarkers(all_Tissue[[i]])
  FindAllMarkers_Tissue_list_5[[length(FindAllMarkers_Tissue_list_5) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Tissue_list_5) <- sample_names

# Remove BM.integrated and PB.integrated objects
FindAllMarkers_Tissue_list_5 <- FindAllMarkers_Tissue_list_5[3:5]

# Annotate the results 
FindAllMarkers_Tissue_list_5 <- lapply(FindAllMarkers_Tissue_list_5, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = BM, pct.2 = PB
FindAllMarkers_Tissue_list_5 <- lapply(FindAllMarkers_Tissue_list_5, remove_PB)

# Save out the DE results
write.xlsx(FindAllMarkers_Tissue_list_5,
           "Output/DE/Controls_1/Workup/5_clusters/FindAllMarkers_Tissue_5_clusters.xlsx")

# Output Volcanoplots for Tissue restricted differences
for (i in 1:length(FindAllMarkers_Tissue_list_5)) {
  temp <- 
    (EnhancedVolcano(FindAllMarkers_Tissue_list_5[[`i`]], 
                     lab = FindAllMarkers_Tissue_list_5[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindAllMarkers_Tissue_list_5)[i],
                     subtitle = "BM vs. PB Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/5_clusters/FindAllMarkers_Tissue_",
                                   names(FindAllMarkers_Tissue_list_5)[i], ".png"), scale = 2:1)
}

# 6.B.iii DE testing - FindAllMarkers_Stim - Assess global differences between Stimulation status --

# Copy the all_list into a new list and switch the identity class to "stimulation_status"
all_Stim <- all_list
all_Stim <- map(all_Stim, `Idents<-`, value = "stimulation_status")

# Create an empty list to store DE results
FindAllMarkers_Stim_list_5 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Stim)){
  temp <- FindAllMarkers(all_Stim[[i]])
  FindAllMarkers_Stim_list_5[[length(FindAllMarkers_Stim_list_5) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Stim_list_5) <- sample_names

# Remove act.integrated and rest.integrated objects
FindAllMarkers_Stim_list_5 <- FindAllMarkers_Stim_list_5[-c(3:4)]

# Annotate the results
FindAllMarkers_Stim_list_5 <- lapply(FindAllMarkers_Stim_list_5, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. rest results) pct.1 = act, pct.2 = rest
FindAllMarkers_Stim_list_5 <- lapply(FindAllMarkers_Stim_list_5, remove_rest)

# Save out the DE results
write.xlsx(FindAllMarkers_Stim_list_5, 
           "Output/DE/Controls_1/Workup/5_clusters/FindAllMarkers_Stim_5_clusters.xlsx")

# 6.B.iv DE testing - FindAllMarkers_Patient - Assess global differences between Patients ----------

# Copy the all_list into a new list and switch the identity class to "orig.ident"
all_Patient <- all_list
all_Patient <- map(all_Patient, `Idents<-`, value = "orig.ident")

# Create an empty list to store DE results
FindAllMarkers_Patient_list_5 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Patient)){
  temp <- FindAllMarkers(all_Patient[[i]])
  FindAllMarkers_Patient_list_5[[length(FindAllMarkers_Patient_list_5) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Patient_list_5) <- sample_names

# Annotate the results 
FindAllMarkers_Patient_list_5 <- lapply(FindAllMarkers_Patient_list_5, annotate_markers)

# Remove redundancy by removing duplicate results (i.e. PB results) pct.1 = Donor1, pct.2 = Donor2
FindAllMarkers_Patient_list_5 <- lapply(FindAllMarkers_Patient_list_5, remove_Donor)


# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_5, 
           "Output/DE/Controls_1/Workup/5_clusters/FindAllMarkers_Patient_5_clusters.xlsx")


# 6.C FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 6.C.i FindMarkers - Itemise differences between Clusters -----------------------------------------

# Create an empty list to store results in
FindMarkers_All_list_5 <- list()

# Run for loop to iterate across all clusters
# NB: This will fail unless each input sample is 5 clusters
for (i in names(all_list)){
  for (j in 1:5) {
    temp <- FindMarkers(all_list[[i]], ident.1 = j-1)
    FindMarkers_All_list_5[[length(FindMarkers_All_list_5) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list. NB: These need to be unique or export to *.xlsx will fail
names(FindMarkers_All_list_5) <- cluster_names_5

# Save out the DE results
write.xlsx(FindMarkers_All_list_5, 
           "Output/DE/Controls_1/Workup/5_clusters/FindMarkers_All_5_clusters.xlsx")

# Output Volcanoplots for all tests
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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/5_clusters/",
                                    names(FindMarkers_All_list_5)[i], ".png"), scale = 2:1)
}

# 6.C.ii FindMarkers - Itemise differences between BM and PB ---------------------------------------

# Create an empty list to store results in
FindMarkers_Tissue_list_5 <- list()

# Copy the all_list object and remove the BM.integrated and PB.integrated objects
all_list_temp <- all_list[-c(1,2)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:5) {
    temp <- FindMarkers(all_list_temp[[`i`]], 
                        ident.1 = "BM", group.by = "Tissue", subset.ident = j-1)
    FindMarkers_Tissue_list_5[[length(FindMarkers_Tissue_list_5) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Tissue_list_5) <- cluster_names_5[11:25]

# Save out the DE results
write.xlsx(FindMarkers_Tissue_list_5, 
           "Output/DE/Controls_1/Workup/5_clusters/FindMarkers_Tissue_5_clusters.xlsx")

# Output Volcanoplots for all tests
for (i in 1:length(FindMarkers_Tissue_list_5)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Tissue_list_5[[`i`]], 
                     lab = FindMarkers_Tissue_list_5[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Tissue_list_5)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 6, 
                     legendPosition = "none", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/5_clusters/FindMarkers_Tissue_",
                                   names(FindMarkers_Tissue_list_5)[i], ".png"), scale = 2:1)
}


# 6.C.iii FindMarkers - Itemise differences between activation and resting -------------------------

# Create an empty list to store results in
FindMarkers_Stim_list_5 <- list()

# Copy the all_list object and remove the act.integrated and rest.integrated objects
all_list_temp <- all_list[-c(3,4)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:5) {
    temp <- FindMarkers(all_list_temp[[`i`]], 
                        ident.1 = "act", group.by = "stimulation_status", subset.ident = j-1)
    FindMarkers_Stim_list_5[[length(FindMarkers_Stim_list_5) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Stim_list_5) <- cluster_names_5[c(1:10, 21:25)]

# Save out the DE results
write.xlsx(FindMarkers_Stim_list_5, 
           "Output/DE/Controls_1/Workup/5_clusters/FindMarkers_Stim_5_clusters.xlsx")

# Output Volcanoplots for all tests
for (i in 1:length(FindMarkers_Stim_list_5)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Stim_list_5[[`i`]], 
                     lab = FindMarkers_Stim_list_5[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Stim_list_5)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_1/Workup/5_clusters/Stim_",
                                   names(FindMarkers_Stim_list_5)[i], ".png"), scale = 2:1)
}


# 6.D Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 6.D.i FindConservedMarkers - Averages across Tissue ----------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Tissue_list_5 <- list()

# Copy the all_list object and remove the BM.integrated and PB.integrated objects
all_list_temp <- all_list[-c(1,2)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:5) {
    temp <- FindConservedMarkers(all_list_temp[[`i`]], 
                                 ident.1 = j-1, grouping.var = "Tissue")
    FindConservedMarkers_Tissue_list_5[[length(FindConservedMarkers_Tissue_list_5) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindConservedMarkers_Tissue_list_5) <- cluster_names_5[11:25]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_5, 
           "Output/DE/Controls_1/Workup/5_clusters/FindConservedMarkers_Tissue_5_clusters.xlsx")

# 6.D.ii FindConservedMarkers - Averages across Stimulation status ---------------------------------

# Create an empty list to store results in
FindConservedMarkers_Stim_list_5 <- list()

# Copy the all_list object and remove the act.integrated and rest.integrated objects
all_list_temp <- all_list[-c(3,4)]

# Run for loop to iterate across all clusters
for (i in names(all_list_temp)){
  for (j in 1:5) {
    temp <- FindConservedMarkers(all_list_temp[[`i`]], 
                                 ident.1 = j-1, grouping.var = "stimulation_status")
    FindConservedMarkers_Stim_list_5[[length(FindConservedMarkers_Stim_list_5) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindConservedMarkers_Stim_list_5) <- cluster_names_5[c(1:10, 21:25)]

# Save out the DE results
write.xlsx(FindConservedMarkers_Tissue_list_5, 
           "Output/DE/Controls_1/Workup/5_clusters/FindConservedMarkers_Stim_5_clusters.xlsx")


# 6.E Create and export figures to assist in cluster identification and annotation -----------------

# 6.E.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_5) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/Controls_1/Workup/5_clusters/",
                  sample_names[i], ".png"))
}

# 6.E.ii Stacked barplots of cluster distribution --------------------------------------------------

# Calculate distribution of clusters and move results into a temp_list
temp_list <- 
  list("temp_BM.integrated" = as.data.frame(table(all_list$BM.integrated$seurat_clusters, 
                                                  all_list$BM.integrated$orig.ident)),
       "temp_PB.integrated" = as.data.frame(table(all_list$PB.integrated$seurat_clusters, 
                                                  all_list$PB.integrated$orig.ident)), 
       "temp_act.integrated" = as.data.frame(table(all_list$act.integrated$seurat_clusters,
                                                   all_list$act.integrated$orig.ident)),
       "temp_rest.integrated" = as.data.frame(table(all_list$rest.integrated$seurat_clusters,
                                                    all_list$rest.integrated$orig.ident)), 
       "temp_all.integrated" = as.data.frame(table(all_list$all.integrated$seurat_clusters,
                                                   all_list$all.integrated$orig.ident)))

# Output barplots of cluster distribution
for (i in 1:length(all_list)) {
  ggplot(temp_list[[i]], 
         aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = rev(ccolours_5), name = "Cluster",
                      (rev(levels(all_list[[i]]@meta.data$seurat_clusters)))) + 
    labs(x = "Sample", y = "Cells", title = sample_names[i]) + theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Barplots/Controls_1/Workup/5_clusters/",
                  sample_names[i], ".png"))
}

# 6.E.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the act.integrated object and output a Dotplot separated by tissue
top_10_markers_act <- 
  FindAllMarkers_All_list_5$act.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_act <- unique(top_10_markers_act$gene)
DotPlot(all_list$act.integrated, 
        features = top_10_unique_act, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "act") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/5_clusters/",
                "Activated BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the rest.integrated object and output a Dotplot separated by tissue
top_10_markers_rest <- 
  FindAllMarkers_All_list_5$rest.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_rest <- unique(top_10_markers_rest$gene)
DotPlot(all_list$rest.integrated, 
        features = top_10_unique_rest, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "rest") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/5_clusters/",
                "Resting BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by tissue
top_10_markers_all <- 
  FindAllMarkers_All_list_5$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "Tissue", cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/5_clusters/",
                "All BM vs. PB", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the BM.integrated object and output a Dotplot separated by activation status
top_10_markers_BM <- 
  FindAllMarkers_All_list_5$BM.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_BM <- unique(top_10_markers_BM$gene)
DotPlot(all_list$BM.integrated, 
        features = top_10_unique_BM, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "BM") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/5_clusters/",
                "BM act vs. rest", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the PB.integrated object and output a Dotplot separated by activation status
top_10_markers_PB <- 
  FindAllMarkers_All_list_5$PB.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_PB <- unique(top_10_markers_PB$gene)
DotPlot(all_list$PB.integrated, 
        features = top_10_unique_PB, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "PB") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/5_clusters/",
                "PB act vs. rest", ".png"), scale = 2:1)

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by activation status
top_10_markers_all <- 
  FindAllMarkers_All_list_5$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "stimulation_status", cols = c("purple", "green")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_1/Workup/5_clusters/",
                "All act vs. rest", ".png"), scale = 2:1)


# 7. Save out data and print stats -----------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean up the environment
rm(list = 
     ls()[!ls() %in% c("annotations", "annotate_markers", "black_list_TCR", 
                       "FindAllMarkers_All_list_5", "FindAllMarkers_All_list_6", 
                       "FindAllMarkers_Patient_list_5", "FindAllMarkers_Patient_list_6", 
                       "FindAllMarkers_Stim_list_5", "FindAllMarkers_Stim_list_6", 
                       "FindAllMarkers_Tissue_list_5", "FindAllMarkers_Tissue_list_6", 
                       "FindConservedMarkers_Stim_list_5", "FindConservedMarkers_Stim_list_6", 
                       "FindConservedMarkers_Tissue_list_5", "FindConservedMarkers_Tissue_list_6", 
                       "FindMarkers_All_list_5", "FindMarkers_All_list_6", 
                       "FindMarkers_Stim_list_5", "FindMarkers_Stim_list_6", 
                       "FindMarkers_Tissue_list_5", "FindMarkers_Tissue_list_6",
                       "keep_TCR", "sample_names")])

# Save out the DE_results for further analysis
save.image("Data/R_out/Controls_1/DE_Workup_Controls_1.rds")

gc()

print("#># Finished running '3. Controls dataset #1 - DE Workup' script")