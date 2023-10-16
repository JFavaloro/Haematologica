# Script information -------------------------------------------------------------------------------

# Title: Controls dataset #2 - DE Workup
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will repeat the DE analysis we performed on our 10x dataset on 
# Controls dataset #2 from Zavidij et al., 2019: https://doi.org/10.1016/j.clml.2019.09.040.
# Where insufficient data exists to allow looping through the analysis without erroring, we will 
# comment these out, or modify the testing and note which are failing or too sparse for analysis.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '5. Control dataset #2 - DE Workup' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, openxlsx, multtest, metap, EnhancedVolcano)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed Controls #2 dataset
integrated_list <- readRDS("Data/R_out/Controls_2/integrated_list.rds")

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

# Define function for annotating markers
# NB: Pass to unique at the end removes redundancy but also seems to remove some genes (?)
annotate_markers <- function(input_markers){
  input_markers <- inner_join(x = input_markers, 
                              y = annotations[, c("gene_name", "description")],
                              by = c("gene" = "gene_name")) %>% unique()
  return(input_markers)
}

# Create character vector of samples names
sample_names <- c("all.integrated")

# Define character vector of plotting order for disease
disease_order = c("HD", "MGUS", "SMM", "MM")

# Change the order of disease for more sensible results
integrated_list$all.integrated$Disease <- 
  factor(integrated_list$all.integrated$Disease, levels = disease_order)

# Change the order of orig.ident for more sensible results
integrated_list$all.integrated$orig.ident <- 
  factor(integrated_list$all.integrated$orig.ident, 
         levels = c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", "HD7", "HD8", "HD9", 
                    "MGUS1", "MGUS2", "MGUS3", "MGUS4", "MGUS5", 
                    "SMM1", "SMM2", "SMM3", "SMM4", "SMM5", "SMM6", "SMM7", "SMM8", "SMM9", "SMM10", "SMM11", 
                    "MM1", "MM2", "MM3", "MM4", "MM5", "MM6", "MM7"))


# Define cluster names for naming of lists - 6 clusters
cluster_names_6 <- apply(expand.grid(c("cluster_"), c(0:5)), 1, paste, collapse="")

# Define cluster names for naming of lists - 6 clusters
cluster_names_5 <- apply(expand.grid(c("cluster_"), c(0:4)), 1, paste, collapse="")

# Create a character vector of colours for plotting clusters on UMAP
# NB: TEM = red ("#F8766D"), TTE = yellow ("#B79F00"), TN = green ("#00BA38"),
# Cyto_TEM = turquoise ("#00BFC4"), PRE_EX = blue ("#619CFF"), TCM = pink ("#F564E3")
# These need to be finalised
scales::show_col(scales::hue_pal()(6))
scales::show_col(scales::hue_pal()(12))
scales::show_col(scales::hue_pal()(25))
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


# 4. Analysis of 6 cluster Seurat objects ----------------------------------------------------------
# 4.A DE testing - FindAllMarkers

# 4.A.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

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
           "Output/DE/Controls_2/Workup/6_clusters/FindAllMarkers_All_6_clusters.xlsx")

# 4.A.ii DE testing - FindAllMarkers_Disease - Assess global differences between disease state -----
# Copy the all integrated objects into a new list and switch the identity class to "Disease"
all_Disease <- all_list
all_Disease <- map(all_Disease, `Idents<-`, value = "Disease")

# Create an empty list to store DE results
FindAllMarkers_Disease_list_6 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Disease)){
  temp <- FindAllMarkers(all_Disease[[i]])
  FindAllMarkers_Disease_list_6[[length(FindAllMarkers_Disease_list_6) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Disease_list_6) <- sample_names

# Annotate the results 
FindAllMarkers_Disease_list_6 <- lapply(FindAllMarkers_Disease_list_6, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Disease_list_6,
           "Output/DE/Controls_2/Workup/6_clusters/FindAllMarkers_Disease_6_clusters.xlsx")

# Output Volcano plots for Disease restricted differences
EnhancedVolcano(FindAllMarkers_Disease_list_6$all.integrated, 
                lab = FindAllMarkers_Disease_list_6$all.integrated$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "Age matched controls vs. Disease",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 12, pointSize = 2, labSize = 4, legendLabSize = 12,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/Controls_2/Workup/6_clusters/",
                "HD vs. Disease", ".png"), scale = 2:1)

# 4.A.iii DE testing - FindAllMarkers_Patient - Assess global differences between Patients ---------

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

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_6, 
           "Output/DE/Controls_2/Workup/6_clusters/FindAllMarkers_Patient_6_clusters.xlsx")


# 4.B FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 4.B.i FindMarkers - Itemise differences between Clusters -----------------------------------------

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
           "Output/DE/Controls_2/Workup/6_clusters/FindMarkers_All_6_clusters.xlsx")

# Output Volcano Plots for all tests
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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_2/Workup/6_clusters/",
                                    names(FindMarkers_All_list_6)[i], ".png"), scale = 2:1)
}

# 4.B.ii FindMarkers - Itemise differences between Disease -----------------------------------------

# Create an empty list to store results in
FindMarkers_Disease_list_6 <- list()

# Run for loop to iterate across all clusters
for (i in names(all_list)){
  for (j in 1:6) {
    temp <- FindMarkers(all_list[[`i`]], 
                        ident.1 = "HD", group.by = "Disease", subset.ident = j-1)
    FindMarkers_Disease_list_6[[length(FindMarkers_Disease_list_6) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Disease_list_6) <- cluster_names_6

# Save out the DE results
write.xlsx(FindMarkers_Disease_list_6, 
           "Output/DE/Controls_2/Workup/6_clusters/FindMarkers_Disease_6_clusters.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(FindMarkers_Disease_list_6)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Disease_list_6[[`i`]], 
                     lab = FindMarkers_Disease_list_6[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Disease_list_6)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_2/Workup/6_clusters/",
                                   names(FindMarkers_Disease_list_6)[i], ".png"), scale = 2:1)
}


# 4.C FindConservedMarkers - Find conserved markers across condition - Itemised --------------------

# 4.C.i FindConservedMarkers - Averages across Disease ---------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Disease_list_6 = list()

for (i in 1:6) {
  temp <- FindConservedMarkers(all_list$all.integrated, 
ident.1 = i-1, grouping.var = "Disease")
  FindConservedMarkers_Disease_list_6[[length(FindConservedMarkers_Disease_list_6) +1]] <- 
    temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("gene" = "gene_name")) %>%
    unique()
}

# Rename the items in the list
names(FindConservedMarkers_Disease_list_6) <- cluster_names_6

# Reorder the results until Sam can explain why this is the case
for (i in names(FindConservedMarkers_Disease_list_6)){
  FindConservedMarkers_Disease_list_6[[`i`]] <- 
    FindConservedMarkers_Disease_list_6[[`i`]][, c(1, 7:11, 2:6, 17:21, 12:16, 22:24)]
}

# Save out the DE results
write.xlsx(FindConservedMarkers_Disease_list_6, 
           "Output/DE/Controls_2/Workup/6_clusters/FindConservedMarkers_Disease_6_clusters.xlsx")


# 4.D Create and export figures to assist in cluster identification and annotation -----------------

# 4.D.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_6) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/Controls_2/Workup/6_clusters/",
                  sample_names[i], ".png"))
}

# 4.D.ii Stacked barplots of cluster distribution --------------------------------------------------

# Create stacked bar plots of cluster distribution on a per disease level
# Gather data to plot
temp_Disease <- 
  as.data.frame(table(all_list$all.integrated$seurat_clusters, all_list$all.integrated$Disease))
# Convert cell numbers to proportions based on disease
temp_Disease = temp_Disease %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells, split by disease status
ggplot(temp_Disease, 
       aes(fill = factor(Var1, levels = 5:0), y = Freq, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_6), name = "Cluster",
                    labels = rev(cluster_names_6)) +
  labs(x = "Patient cohort", y = "Cells", 
       title = "Cluster distribution by disease status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/6_clusters/",
                "Cluster_by_disease", ".png"))

# Plot as proportion of total, split by disease
ggplot(temp_Disease, 
       aes(fill = factor(Var1, levels = 5:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_6), name = "Cluster", 
                    labels = rev(cluster_names_6)) + 
  labs(x = "Patient cohort", y = "Proportion (%)", 
       title = "Cluster distribution by stimulation status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/6_clusters/",
                "Cluster_by_disease_proportions", ".png"))

# Create stacked bar plots of cluster distribution on a per sample level
# Gather data to plot
temp_all.integrated <- 
  as.data.frame(table(all_list$all.integrated$seurat_clusters, all_list$all.integrated$orig.ident))
# Convert cell numbers to proportions based on disease
temp_all.integrated = temp_all.integrated %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells, split by sample
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_6), name = "Cluster", 
                    labels = rev(cluster_names_6),
                    (rev(levels(all_list$all.integrated@meta.data$seurat_clusters)))) + 
  geom_text(aes(label = paste(rev(Freq))), position = position_stack(vjust = 0.5)) +
  labs(x = "Sample", y = "Cells", title = "all.integrated") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/6_clusters/",
                "all.integrated", ".png"))

# Plot as proportion of total, split by sample
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Proportion), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_6), name = "Cluster", 
                    labels = rev(cluster_names_6),
                    (rev(levels(all_list$all.integrated@meta.data$seurat_clusters)))) + 
  labs(x = "Sample", y = "Proportion (%)", title = "all.integrated") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/6_clusters/",
                "all.integrated_proportion", ".png"))

# 4.D.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by Disease
top_10_markers_all <- 
  FindAllMarkers_All_list_6$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "Disease", 
        cols = c("orange", "turquoise", "green", "purple")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_2/Workup/6_clusters/",
                "Disease", ".png"), scale = 2:1)


# 5. Analysis of 5 cluster Seurat objects ----------------------------------------------------------

# 5.A Set resolutions of Seurat objects to allow clustering into 5 clusters

# Set resolution to 0.20 for all.integrated object
all_list$all.integrated$seurat_clusters <- 
  all_list$all.integrated$integrated_snn_res.0.2
Idents(all_list$all.integrated) <- "seurat_clusters"

# 5.B DE testing - FindAllMarkers

# 5.B.i DE testing - FindAllMarkers - Compile all genes that define our clusters -------------------

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
           "Output/DE/Controls_2/Workup/5_clusters/FindAllMarkers_All_5_clusters.xlsx")

# 5.B.ii DE testing - FindAllMarkers_Disease - Assess global differences between disease state ------ 
# Copy the all integrated objects into a new list and switch the identity class to "Disease"
all_Disease <- all_list
all_Disease <- map(all_Disease, `Idents<-`, value = "Disease")

# Create an empty list to store DE results
FindAllMarkers_Disease_list_5 <- list()

# Run FindAllMarkers and store results in a list
for (i in names(all_Disease)){
  temp <- FindAllMarkers(all_Disease[[i]])
  FindAllMarkers_Disease_list_5[[length(FindAllMarkers_Disease_list_5) +1]] <- temp
}

# Rename the objects in the list
names(FindAllMarkers_Disease_list_5) <- sample_names

# Annotate the results 
FindAllMarkers_Disease_list_5 <- lapply(FindAllMarkers_Disease_list_5, annotate_markers)

# Save out the DE results
write.xlsx(FindAllMarkers_Disease_list_5,
           "Output/DE/Controls_2/Workup/5_clusters/FindAllMarkers_Disease_5_clusters.xlsx")

# Output Volcano plots for Disease restricted differences
EnhancedVolcano(FindAllMarkers_Disease_list_5$all.integrated, 
                lab = FindAllMarkers_Disease_list_5$all.integrated$gene, 
                x = 'avg_log2FC', y = 'p_val', title = "Age matched controls vs. Disease",
                subtitle = "Differential Expression", drawConnectors = TRUE,
                axisLabSize = 12, pointSize = 2, labSize = 4, legendLabSize = 12,
                legendPosition = "bottom", max.overlaps = Inf)
ggsave(filename = 
         paste0("Output/Figures/Volcanoplots/Controls_2/Workup/5_clusters/",
                "HD vs. Disease", ".png"), scale = 2:1)

# 5.B.iii DE testing - FindAllMarkers_Patient - Assess global differences between Patients ---------

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

# Save out the DE results
write.xlsx(FindAllMarkers_Patient_list_5, 
           "Output/DE/Controls_2/Workup/5_clusters/FindAllMarkers_Patient_5_clusters.xlsx")


# 5.C FindMarkers - Find markers between cluster "A" and "B" - Itemised ----------------------------

# 5.C.i FindMarkers - Itemise differences between Clusters -----------------------------------------

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
           "Output/DE/Controls_2/Workup/5_clusters/FindMarkers_All_5_clusters.xlsx")

# Output Volcano Plots for all tests
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
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_2/Workup/",
                                    names(FindMarkers_All_list_5)[i], ".png"), scale = 2:1)
}

# 5.C.ii FindMarkers - Itemise differences between Disease -----------------------------------------

# Create an empty list to store results in
FindMarkers_Disease_list_5 <- list()

# Run for loop to iterate across all clusters
for (i in names(all_list)){
  for (j in 1:5) {
    temp <- FindMarkers(all_list[[`i`]], 
                        ident.1 = "HD", group.by = "Disease", subset.ident = j-1)
    FindMarkers_Disease_list_5[[length(FindMarkers_Disease_list_5) +1]] <- 
      temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("gene" = "gene_name")) %>%
      unique()
  }
}

# Rename the items in the list
names(FindMarkers_Disease_list_5) <- cluster_names_5

# Save out the DE results
write.xlsx(FindMarkers_Disease_list_5, 
           "Output/DE/Controls_2/Workup/5_clusters/FindMarkers_Disease_5_clusters.xlsx")

# Output Volcano Plots for all tests
for (i in 1:length(FindMarkers_Disease_list_5)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_Disease_list_5[[`i`]], 
                     lab = FindMarkers_Disease_list_5[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_Disease_list_5)[i],
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/Controls_2/Workup/5_clusters/Disease_",
                                   names(FindMarkers_Disease_list_5)[i], ".png"), scale = 2:1)
}


# 5.D Find ConservedMarkers - Find conserved markers across condition - Itemised -------------------

# 5.D.i FindConservedMarkers - Averages across Disease ---------------------------------------------

# Create an empty list to store results in
FindConservedMarkers_Disease_list_5 = list()

# Run for loop to iterate across all clusters
for (i in 1:5) {
  temp <- FindConservedMarkers(all_list$all.integrated, 
                               ident.1 = j-i, grouping.var = "Disease")
  FindConservedMarkers_Disease_list_5[[length(FindConservedMarkers_Disease_list_5) +1]] <- 
    temp %>% rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("gene" = "gene_name")) %>%
    unique()
}


# Rename the items in the list
names(FindConservedMarkers_Disease_list_5) <- cluster_names_5

# Reorder the results until Sam can explain why this is the case
for (i in names(FindConservedMarkers_Disease_list_5)){
  FindConservedMarkers_Disease_list_5[[`i`]] <- 
    FindConservedMarkers_Disease_list_5[[`i`]][, c(1, 7:11, 2:6, 17:21, 12:16, 22:24)]
}

# Save out the DE results
write.xlsx(FindConservedMarkers_Disease_list_5, 
           "Output/DE/Controls_2/Workup/5_clusters/FindConservedMarkers_Disease_5_clusters.xlsx")


# 5.E Create and export figures to assist in cluster identification and annotation -----------------

# 5.E.i Dimplots of cluster distribution -----------------------------------------------------------

# UMAPs of cluster distribution
for (i in 1:length(sample_names)) {
  DimPlot(all_list[[i]], cols = ccolours_5) + labs(title = sample_names[i], color = "Cluster") + 
    theme_bw()
  ggsave(filename = 
           paste0("Output/Figures/Dimplots/Controls_2/Workup/5_clusters/",
                  sample_names[i], ".png"))
}

# 5.E.ii Stacked barplots of cluster distribution --------------------------------------------------

# Create stacked bar plots of cluster distribution on a per disease level
# Gather data to plot
temp_Disease <- 
  as.data.frame(table(all_list$all.integrated$seurat_clusters, all_list$all.integrated$Disease))
# Convert cell numbers to proportions based on disease
temp_Disease = temp_Disease %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells, split by disease status
ggplot(temp_Disease, 
       aes(fill = factor(Var1, levels = 4:0), y = Freq, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_5), name = "Cluster",
                    labels = rev(cluster_names_5)) +
  labs(x = "Patient cohort", y = "Cells", 
       title = "Cluster distribution by disease status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/5_clusters/",
                "Cluster_by_disease", ".png"))

# Plot as proportion of total, split by disease
ggplot(temp_Disease, 
       aes(fill = factor(Var1, levels = 4:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_5), name = "Cluster", 
                    labels = rev(cluster_names_5)) + 
  labs(x = "Patient cohort", y = "Proportion (%)", 
       title = "Cluster distribution by stimulation status") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/5_clusters/",
                "Cluster_by_disease_proportions", ".png"))

# Create stacked bar plots of cluster distribution on a per sample level
# Gather data to plot
temp_all.integrated <- 
  as.data.frame(table(all_list$all.integrated$seurat_clusters, all_list$all.integrated$orig.ident))
# Convert cell numbers to proportions based on disease
temp_all.integrated = temp_all.integrated %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Plot as total cells, split by sample
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Freq), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_5), name = "Cluster", 
                    labels = rev(cluster_names_5),
                    (rev(levels(all_list$all.integrated@meta.data$seurat_clusters)))) + 
  geom_text(aes(label = paste(rev(Freq))), position = position_stack(vjust = 0.5)) +
  labs(x = "Sample", y = "Cells", 
       title = "all.integrated") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/5_clusters/",
                "all.integrated", ".png"))

# Plot as proportion of total, split by sample
ggplot(temp_all.integrated, 
       aes(fill = Var1, y = rev(Proportion), x = rev(Var2))) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_5), name = "Cluster", 
                    labels = rev(cluster_names_5),
                    (rev(levels(all_list$all.integrated@meta.data$seurat_clusters)))) + 
  labs(x = "Sample", y = "Proportion (%)", 
       title = "all.integrated") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Workup/5_clusters/",
                "all.integrated_proportion", ".png"))

# 5.E.iii Dotplots of genes that define clusters based on level of expression ----------------------

# Extract the top 10 unique markers for the all.integrated object and output a Dotplot separated by Disease
top_10_markers_all <- 
  FindAllMarkers_All_list_5$all.integrated %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all <- unique(top_10_markers_all$gene)
DotPlot(all_list$all.integrated, 
        features = top_10_unique_all, split.by = "Disease", 
        cols = c("orange", "turquoise", "green", "purple")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Gene", y = "Cluster", title = "all") + NoLegend()
ggsave(filename = 
         paste0("Output/Figures/Dotplots/Controls_2/Workup/5_clusters/",
                "Disease", ".png"), scale = 2:1)


# 6. Save out data and print stats -----------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean up the environment
rm(list = 
     ls()[!ls() %in% c("annotations", "annotate_markers", "black_list_TCR", 
                       "FindAllMarkers_All_list_5", "FindAllMarkers_All_list_6", 
                       "FindAllMarkers_Disease_list_6", "FindAllMarkers_Disease_list_6", 
                       "FindAllMarkers_Patient_list_5", "FindAllMarkers_Patient_list_6", 
                       "FindAllMarkers_Disease_list_5", "FindAllMarkers_Disease_list_6", 
                       "FindConservedMarkers_Disease_list_5", "FindConservedMarkers_Disease_list_6", 
                       "FindMarkers_All_list_5", "FindMarkers_All_list_6", 
                       "FindMarkers_Disease_list_5", "FindMarkers_Disease_list_6",
                       "keep_TCR", "sample_names")])

# Save out the DE_results for further analysis
save.image("Data/R_out/Controls_2/DE_Workup_Controls_2.rds")

print("#># Finished running '5. Control dataset #2 - DE Workup' script")