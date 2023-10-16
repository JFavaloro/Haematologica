# Script information -------------------------------------------------------------------------------

# Title: Further analysis of projected Controls data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will continue the analysis of the ProjectTILs projection of public 
# data onto our custom reference maps. The DE results will be utilised to prepare figures comparing
# expression levels of a number of genes of interest illustrated through radar plots. This script is
# reasonably modular - the genes_to_plot list can be amended if required - max 13 genes per column.
# NB: Aesthetics of the radar plots can be edited with trace(plot.states.radar, edit = "true")
# Change 'black' to 'blue' for BM and 'red' for PB atlas. Change lines X:
# annotate(geom = "text", fontface = "italic", x = seq(1, 
# length(genes4radar)), y = ymax - 0.05 * ymax, 
# label = genes4radar, size = 7) + coord_polar()


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '8. Further analysis of projected Controls data' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire, openxlsx, multtest, metap, EnhancedVolcano, ProjecTILs, ggpubr)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed 10x and public data
all_list <- readRDS("Data/R_out/10x/all_list.rds")
Controls_1_complete <- readRDS("Data/R_out/Controls_1/Controls_1_complete.rds")
Controls_2_complete <- readRDS("Data/R_out/Controls_2/Controls_2_complete.rds")
Controls_3_complete <- readRDS("Data/R_out/Controls_3/Controls_3_complete.rds")

# Load custom reference atlases
all_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/all_10x_atlas.rds")
BM_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/BM_10x_atlas.rds")
PB_10x_atlas <- load.reference.map("Data/R_out/ProjecTILs/PB_10x_atlas.rds")

# Load the Controls projections
Controls_query_list <- readRDS("Data/R_out/ProjecTILs/Controls_query_list.rds")

# Load the 10x atlases DE analysis
load("Data/R_out/ProjecTILs/DE_Analysis_10x_atlases.rds")

# Load the Projections DE analysis
load("Data/R_out/ProjecTILs/Projections_DE_Analysis.rds")

# Load additional resources
genes_to_plot <- read.csv("Data/Public_data/Genes_for_plotting.csv", na.strings = "") %>% 
  as.list() %>% lapply(na.omit)

# Create character vectors of cluster names
cluster_names_all_10x_atlas <- c("T[EM]", "T[TE]", "T[N]", "IL7R^`+` ~ T[M]", "Cyto ~ T[EM]", "Activated")
cluster_names_BM_10x_atlas <- c("T[EM]", "T[N]", "IL7R^`+` ~ T[M]", "T[TE]", "Activated", "Cyto ~ T[EM]")
cluster_names_PB_10x_atlas <- c("T[TE]", "T[EM]", "IL7R^`+` ~ T[M]", "T[N]", "Activated", "IFN")

excel_cluster_names_all_10x_atlas <- c("TEM", "TTE", "TN", "IL7R_TM", "Cyto_TEM", "Activated")
excel_cluster_names_BM_10x_atlas <- c("TEM", "TN", "IL7R_TM", "TTE", "Activated", "Cyto_TEM")
excel_cluster_names_PB_10x_atlas <- c("TTE", "TEM", "IL7R_TM", "TN", "Activated", "IFN")

# Create character vectors of title plots
Controls_1_titles <- c("Controls 1", "Activated", "Resting", 
                       "Activated BM", "Resting BM", "Activated PB", "Resting PB")
Controls_2_titles <- c("Controls 2", "Age-matched BM", "MGUS BM", "SMM BM", "MM BM")
Controls_3_titles <- c("Controls 3", "Age-matched PB", "Young PB")

# Create character vectors of colours for plotting clusters
all_10x_atlas_plot_cols <- 
  c("T[EM]" = "#F8766D", "T[TE]" = "#B79F00", "T[N]" = "#00BA38", 
    "IL7R^`+` ~ T[M]" = "#F564E3", "Cyto ~ T[EM]" = "#00BFC4", "Activated" = "#8A2BE2")

BM_10x_atlas_plot_cols <- 
  c("T[EM]" = "#F8766D", "T[N]" = "#00BA38", "IL7R^`+` ~ T[M]" = "#F564E3", 
    "T[TE]" = "#B79F00", "Activated" = "#8A2BE2", "Cyto ~ T[EM]" = "#00BFC4")

PB_10x_atlas_plot_cols <- 
  c("T[TE]" = "#B79F00", "T[EM]" = "#F8766D", "IL7R^`+` ~ T[M]" = "#F564E3", 
    "T[N]" = "#00BA38", "Activated" = "#8A2BE2", "IFN" = "#FF69B4")

ccolours_all_10x_atlas <- c("#F8766D", "#B79F00", "#00BA38", "#F564E3", "#00BFC4", "#8A2BE2")
ccolours_BM_10x_atlas <- c("#F8766D", "#00BA38", "#F564E3", "#B79F00", "#8A2BE2", "#00BFC4")
ccolours_PB_10x_atlas <- c("#B79F00", "#F8766D", "#F564E3", "#00BA38", "#8A2BE2", "#FF69B4")

# Create character vectors of title plots
Controls_1_titles <- c("Controls 1", "Activated", "Resting", 
                       "Activated BM", "Resting BM", "Activated PB", "Resting PB")
Controls_2_titles <- c("Controls 2", "Age-matched BM", "MGUS BM", "SMM BM", "MM BM")
Controls_3_titles <- c("Controls 3", "Age-matched PB", "Young PB")


# 3. Pre-processing --------------------------------------------------------------------------------

# Move atlases into a list
all_atlas_list <- list("all_10x_atlas" = all_10x_atlas,
                       "BM_10x_atlas" = BM_10x_atlas,
                       "PB_10x_atlas" = PB_10x_atlas)


# 3.1 Create lists of the top 10 genes defining functional.clusters to plot on radar plots ---------

# Top 10 genes defining clusters in the all_10x_atlas
top_10_FindAllMarkers_all_10x_atlas <- 
  FindAllMarkers_All_list$all_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_FindAllMarkers_all_10x_atlas <- 
  split(top_10_FindAllMarkers_all_10x_atlas, 
        factor(sort(rank(row.names(top_10_FindAllMarkers_all_10x_atlas))%%6)))

# Top 10 genes defining clusters in the BM_10x_atlas
top_10_FindAllMarkers_BM_10x_atlas <- 
  FindAllMarkers_All_list$BM_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_FindAllMarkers_BM_10x_atlas <- 
  split(top_10_FindAllMarkers_BM_10x_atlas, 
        factor(sort(rank(row.names(top_10_FindAllMarkers_BM_10x_atlas))%%6)))

# Top 10 genes defining clusters in the PB_10x_atlas
top_10_FindAllMarkers_PB_10x_atlas <- 
  FindAllMarkers_All_list$PB_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_FindAllMarkers_PB_10x_atlas <- 
  split(top_10_FindAllMarkers_PB_10x_atlas, 
        factor(sort(rank(row.names(top_10_FindAllMarkers_PB_10x_atlas))%%6)))


# 3.2 Create lists of the top 10 genes defining cluster restricted differences ---------------------

# Create lists of the top 10 genes that distinguish Controls #1 from the all_10x reference atlas
top_10_FindMarkers_Controls_1_all_10x_atlas <- list()

for (i in names(Controls_1_vs_all_10x_atlas_FindMarkers_All)) {
  temp <- Controls_1_vs_all_10x_atlas_FindMarkers_All[[i]] %>% top_n(10, avg_log2FC)
  top_10_FindMarkers_Controls_1_all_10x_atlas[[length(top_10_FindMarkers_Controls_1_all_10x_atlas) +1]] <- 
    temp
}

# Rename the list
names(top_10_FindMarkers_Controls_1_all_10x_atlas) <- 
  names(Controls_1_vs_all_10x_atlas_FindMarkers_All)

# Create lists of the top 10 genes that distinguish Controls #1 from the BM_10x reference atlas
top_10_FindMarkers_Controls_1_BM_10x_atlas <- list()

for (i in names(Controls_1_vs_BM_10x_atlas_FindMarkers_All)) {
  temp <- Controls_1_vs_BM_10x_atlas_FindMarkers_All[[i]] %>% top_n(10, avg_log2FC)
  top_10_FindMarkers_Controls_1_BM_10x_atlas[[length(top_10_FindMarkers_Controls_1_BM_10x_atlas) +1]] <- 
    temp
}

# Rename the list
names(top_10_FindMarkers_Controls_1_BM_10x_atlas) <- 
  names(Controls_1_vs_BM_10x_atlas_FindMarkers_All)

# Create lists of the top 10 genes that distinguish Controls #1 from the PB_10x reference atlas
top_10_FindMarkers_Controls_1_PB_10x_atlas <- list()

for (i in names(Controls_1_vs_PB_10x_atlas_FindMarkers_All)) {
  temp <- Controls_1_vs_PB_10x_atlas_FindMarkers_All[[i]] %>% top_n(10, avg_log2FC)
  top_10_FindMarkers_Controls_1_PB_10x_atlas[[length(top_10_FindMarkers_Controls_1_PB_10x_atlas) +1]] <- 
    temp
}

# Rename the list
names(top_10_FindMarkers_Controls_1_PB_10x_atlas) <- 
  names(Controls_1_vs_PB_10x_atlas_FindMarkers_All)

# Create lists of the top 10 genes that distinguish Controls #2 from the BM_10x reference atlas
top_10_FindMarkers_Controls_2_BM_10x_atlas = list()

for (i in names(Controls_2_vs_BM_10x_atlas_FindMarkers_All)) {
  temp = Controls_2_vs_BM_10x_atlas_FindMarkers_All[[i]] %>% top_n(10, avg_log2FC)
  top_10_FindMarkers_Controls_2_BM_10x_atlas[[length(top_10_FindMarkers_Controls_2_BM_10x_atlas) +1]] <- 
    temp
}

# Rename the list
names(top_10_FindMarkers_Controls_2_BM_10x_atlas) <- 
  names(Controls_2_vs_BM_10x_atlas_FindMarkers_All)

# Create lists of the top 10 genes that distinguish Controls #3 from the PB_10x reference atlas
top_10_FindMarkers_Controls_3_PB_10x_atlas = list()

for (i in names(Controls_3_vs_PB_10x_atlas_FindMarkers_All)) {
  temp = Controls_3_vs_PB_10x_atlas_FindMarkers_All[[i]] %>% top_n(10, avg_log2FC)
  top_10_FindMarkers_Controls_3_PB_10x_atlas[[length(top_10_FindMarkers_Controls_3_PB_10x_atlas) +1]] <- 
    temp
}

# Rename the list
names(top_10_FindMarkers_Controls_3_PB_10x_atlas) <- 
  names(Controls_3_vs_PB_10x_atlas_FindMarkers_All)


# 4.2 Plot radar -----------------------------------------------------------------------------------

# 4.2.A.i Plot Controls 1 - all atlas radar plots --------------------------------------------------

# Create an empty list to store results
Controls_1_all_radarplots <- list()

# Run a loop plotting all combinations of the Controls_1 object against the all_10_x_atlas
for (i in names(Controls_query_list[["Controls_1"]][1])) {
  for (j in names(top_10_FindAllMarkers_all_10x_atlas)) {
    temp <- 
      (plot.states.radar(ref = all_10x_atlas, 
                         query = list("Activated" = Controls_query_list[["Controls_1"]][[2]],
                                      "Resting" = Controls_query_list[["Controls_1"]][[3]]),
                         genes4radar = top_10_FindAllMarkers_all_10x_atlas[[`j`]]$gene,
                         cols = c("purple", "green"),
                         return.as.list = TRUE))
    Controls_1_all_radarplots[[length(Controls_1_all_radarplots) +1]] <- temp
  }
}

# Extract the images we're interested in
Controls_1_all_radarplots <- list("T[EM]" = Controls_1_all_radarplots[[1]][[1]],
                                  "T[TE]" = Controls_1_all_radarplots[[2]][[2]],
                                  "T[N]" = Controls_1_all_radarplots[[3]][[3]],
                                  "IL7R^`+` ~ T[M]" = Controls_1_all_radarplots[[4]][[4]],
                                  "Cyto ~ T[EM]" = Controls_1_all_radarplots[[5]][[5]],
                                  "Activated" = Controls_1_all_radarplots[[6]][[6]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_1_all_radarplots[[i]] + theme(text = element_text(size = 12)) + 
           ggtitle(label = "Combined", subtitle = parse(text = names(Controls_1_all_radarplots[i]))), 
         filename = paste0("Output/Figures/Radarplots/Controls_1/Final/all_10x_atlas/Genes_defining_cluster/",
                           names(Controls_1_all_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_1_all_radarplots <- list()

# Run a loop plotting all combinations of the Controls_1 object against the all_10_x_atlas
for (i in names(top_10_FindMarkers_Controls_1_all_10x_atlas)) {
  temp <- 
    (plot.states.radar(ref = all_10x_atlas, 
                       query = list("Activated" = Controls_query_list[["Controls_1"]][[2]],
                                    "Resting" = Controls_query_list[["Controls_1"]][[3]]),
                       genes4radar = top_10_FindMarkers_Controls_1_all_10x_atlas[[`i`]]$gene,
                       cols = c("purple", "green"),  
                       return.as.list = TRUE))
  Controls_1_all_radarplots[[length(Controls_1_all_radarplots) +1]] <- temp
}

# Extract the images we're interested in - [Sample][State]
Controls_1_all_radarplots <- list("Activated ~ T[EM]" = Controls_1_all_radarplots[[4]][[1]],
                                  "Activated ~ T[TE]" = Controls_1_all_radarplots[[5]][[2]],
                                  "Activated ~ T[N]" = Controls_1_all_radarplots[[6]][[3]],
                                  "Resting ~ T[EM]" = Controls_1_all_radarplots[[7]][[1]],
                                  "Resting ~ T[TE]" = Controls_1_all_radarplots[[8]][[2]],
                                  "Resting ~ T[N]" = Controls_1_all_radarplots[[9]][[3]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_1_all_radarplots[[i]] + theme(text = element_text(size = 12)) + 
           ggtitle(label = "Combined", subtitle = parse(text = names(Controls_1_all_radarplots[i]))),
         filename = paste0("Output/Figures/Radarplots/Controls_1/Final/all_10x_atlas/Genes_defining_state/",
                           names(Controls_1_all_radarplots[i]), ".png"),  
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_1_all_radarplots <- list()

# Run a loop plotting genes of interest in the Controls_1 object against the all_10_x_atlas
for (i in names(Controls_query_list[["Controls_1"]][1])) {
  for (j in names(genes_to_plot)) {
    temp <- 
      (plot.states.radar(ref = all_10x_atlas, 
                         query = list("Activated" = Controls_query_list[["Controls_1"]][[2]],
                                      "Resting" = Controls_query_list[["Controls_1"]][[3]]),
                         genes4radar = genes_to_plot[[`j`]],
                         cols = c("purple", "green"),  
                         return.as.list = TRUE))
    Controls_1_all_radarplots[[length(Controls_1_all_radarplots) +1]] <- temp
  }
}

# Rename the objects in the list - [Sample][State]
for (i in 1:8) {
  names(Controls_1_all_radarplots[[i]]) <- cluster_names_all_10x_atlas
}
names(Controls_1_all_radarplots) <- c("Exhaustion", "Activation", "Effector", "Regulators", 
                                      "Residency", "Homing_1", "Homing_2", "Homing_3")

# Export the images
for(i in names(Controls_1_all_radarplots)){
  for (j in names(Controls_1_all_radarplots[[1]])) {
    ggsave(plot = 
             Controls_1_all_radarplots[[i]][[j]] + theme(text = element_text(size = 12)) +
             ggtitle(label = "Combined", subtitle = parse(text = names(Controls_1_all_radarplots[[1]][j]))), 
           filename = paste0("Output/Figures/Radarplots/Controls_1/Final/all_10x_atlas/Genes_of_interest/",
                             names(Controls_1_all_radarplots[i]), "/all_10x_atlas_", j, ".png"), 
           units = "cm", height = 16, width = 18)
  }
}


# 4.2A ii Plot Controls 1 - BM atlas radar plots ---------------------------------------------------

# Create an empty list to store results
Controls_1_BM_radarplots <- list()

# Run a loop plotting BM combinations of the Controls_1 object against the BM_10_x_atlas
for (i in names(Controls_query_list[["Controls_1"]][1])) {
  for (j in names(top_10_FindAllMarkers_BM_10x_atlas)) {
    temp <- 
      (plot.states.radar(ref = BM_10x_atlas, 
                         query = list("Activated BM" = Controls_query_list[["Controls_1"]][[4]],
                                      "Resting BM" = Controls_query_list[["Controls_1"]][[5]]),
                         genes4radar = top_10_FindAllMarkers_BM_10x_atlas[[`j`]]$gene,
                         cols = c("purple", "green"),  
                         return.as.list = TRUE))
    Controls_1_BM_radarplots[[length(Controls_1_BM_radarplots) +1]] <- temp
  }
}

# Extract the images we're interested in
Controls_1_BM_radarplots <- list("T[EM]" = Controls_1_BM_radarplots[[1]][[1]],
                                 "T[N]" = Controls_1_BM_radarplots[[2]][[2]],
                                 "IL7R^`+` ~ T[M]" = Controls_1_BM_radarplots[[3]][[3]], # Rename this to IL7R
                                 "T[TE]" = Controls_1_BM_radarplots[[4]][[4]],
                                 "Activated" = Controls_1_BM_radarplots[[5]][[5]],
                                 "Cyto ~ T[EM]" = Controls_1_BM_radarplots[[6]][[6]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_1_BM_radarplots[[i]] + theme(text = element_text(size = 12)) +  
           ggtitle(label = "BM", subtitle = parse(text = names(Controls_1_BM_radarplots[i]))),
         filename = paste0("Output/Figures/Radarplots/Controls_1/Final/BM_10x_atlas/Genes_defining_cluster/",
                           names(Controls_1_BM_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_1_BM_radarplots <- list()

# Run a loop plotting all combinations of the Controls_1 object against the BM_10_x_atlas
for (i in names(top_10_FindMarkers_Controls_1_BM_10x_atlas)) {
  temp <- 
    (plot.states.radar(ref = BM_10x_atlas, 
                       query = list("Activated BM" = Controls_query_list[["Controls_1"]][[4]],
                                    "Resting BM" = Controls_query_list[["Controls_1"]][[5]]),
                       genes4radar = top_10_FindMarkers_Controls_1_BM_10x_atlas[[`i`]]$gene,
                       cols = c("purple", "green"),  
                       return.as.list = TRUE))
  Controls_1_BM_radarplots[[length(Controls_1_BM_radarplots) +1]] <- temp
}

# Extract the images we're interested in - [Sample][State]
Controls_1_BM_radarplots <- list("Activated ~ BM ~ T[EM]" = Controls_1_BM_radarplots[[1]][[1]],
                                 "Activated ~ BM ~ IL7R^`+` ~ T[M]" = Controls_1_BM_radarplots[[2]][[3]],
                                 "Activated ~ BM ~ T[TE]" = Controls_1_BM_radarplots[[3]][[4]],
                                 "Resting ~ BM ~ T[EM]" = Controls_1_BM_radarplots[[4]][[1]],
                                 "Resting ~ BM ~ IL7R^`+` ~ T[M]" = Controls_1_BM_radarplots[[5]][[3]],
                                 "Resting ~ BM ~ T[TE]" = Controls_1_BM_radarplots[[6]][[4]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_1_BM_radarplots[[i]]  + theme(text = element_text(size = 12)) +  
           ggtitle(label = "BM", subtitle = parse(text = names(Controls_1_BM_radarplots[i]))),
         filename = paste0("Output/Figures/Radarplots/Controls_1/Final/BM_10x_atlas/Genes_defining_state/",
                           names(Controls_1_BM_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_1_BM_radarplots <- list()

# Run a loop plotting genes of interest in the Controls_1 object against the BM_10_x_atlas
for (i in names(Controls_query_list[["Controls_1"]][1])) {
  for (j in names(genes_to_plot)) {
    temp <- 
      (plot.states.radar(ref = BM_10x_atlas, 
                         query = list("Activated BM" = Controls_query_list[["Controls_1"]][[4]],
                                      "Resting BM" = Controls_query_list[["Controls_1"]][[5]]),
                         genes4radar = genes_to_plot[[`j`]],
                         cols = c("purple", "green"),  
                         return.as.list = TRUE))
    Controls_1_BM_radarplots[[length(Controls_1_BM_radarplots) +1]] <- temp
  }
}

# Rename the objects in the list - [Sample][State]
for (i in 1:8) {
  names(Controls_1_BM_radarplots[[i]]) <- cluster_names_BM_10x_atlas
}
names(Controls_1_BM_radarplots) <- c("Exhaustion", "Activation", "Effector", "Regulators", 
                                     "Residency", "Homing_1", "Homing_2", "Homing_3")

# Export the images
for(i in names(Controls_1_BM_radarplots)){
  for (j in names(Controls_1_BM_radarplots[[1]])) {
    ggsave(plot = 
             Controls_1_BM_radarplots[[i]][[j]] + theme(text = element_text(size = 12)) +  
             ggtitle(label = "BM", subtitle = parse(text = names(Controls_1_BM_radarplots[[1]][j]))),
           filename = paste0("Output/Figures/Radarplots/Controls_1/Final/BM_10x_atlas/Genes_of_interest/",
                             names(Controls_1_BM_radarplots[i]), "/BM_10x_atlas_", j, ".png"), 
           units = "cm", height = 16, width = 18)
  }
}


# 4.2A iii Plot Controls 1 - PB atlas radar plots --------------------------------------------------

# Create an empty list to store results
Controls_1_PB_radarplots <- list()

# Run a loop plotting PB combinations of the Controls_1 object against the PB_10_x_atlas
for (i in names(Controls_query_list[["Controls_1"]][1])) {
  for (j in names(top_10_FindAllMarkers_PB_10x_atlas)) {
    temp <- 
      (plot.states.radar(ref = PB_10x_atlas, 
                         query = list("Activated PB" = Controls_query_list[["Controls_1"]][[6]],
                                      "Resting PB" = Controls_query_list[["Controls_1"]][[7]]),
                         genes4radar = top_10_FindAllMarkers_PB_10x_atlas[[`j`]]$gene,
                         cols = c("purple", "green"),  
                         return.as.list = TRUE))
    Controls_1_PB_radarplots[[length(Controls_1_PB_radarplots) +1]] <- temp
  }
}

# Extract the images we're interested in
Controls_1_PB_radarplots <- list("T[TE]" = Controls_1_PB_radarplots[[1]][[1]],
                                 "T[EM]" = Controls_1_PB_radarplots[[2]][[2]],
                                 "IL7R^`+` ~ T[M]" = Controls_1_PB_radarplots[[3]][[3]],
                                 "T[N]" = Controls_1_PB_radarplots[[4]][[4]],
                                 "Activated" = Controls_1_PB_radarplots[[5]][[5]],
                                 "IFN" = Controls_1_PB_radarplots[[6]][[6]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_1_PB_radarplots[[i]] + theme(text = element_text(size = 12)) +  
           ggtitle(label = "PB", subtitle = parse(text = names(Controls_1_PB_radarplots[i]))),
         filename = paste0("Output/Figures/Radarplots/Controls_1/Final/PB_10x_atlas/Genes_defining_cluster/",
                           names(Controls_1_PB_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_1_PB_radarplots <- list()

# Run a loop plotting all combinations of the Controls_1 object against the PB_10_x_atlas
for (i in names(top_10_FindMarkers_Controls_1_PB_10x_atlas)) {
  temp <- 
    (plot.states.radar(ref = PB_10x_atlas, 
                       query = list("Activated PB" = Controls_query_list[["Controls_1"]][[4]],
                                    "Resting PB" = Controls_query_list[["Controls_1"]][[5]]),
                       genes4radar = top_10_FindMarkers_Controls_1_PB_10x_atlas[[`i`]]$gene,
                       cols = c("purple", "green"),  
                       return.as.list = TRUE))
  Controls_1_PB_radarplots[[length(Controls_1_PB_radarplots) +1]] <- temp
}

# Extract the images we're interested in - [Sample][State]
Controls_1_PB_radarplots <- list("Activated ~ PB ~ T[TE]" = Controls_1_PB_radarplots[[1]][[1]],
                                 "Activated ~ PB ~ T[EM]" = Controls_1_PB_radarplots[[2]][[2]],
                                 "Activated ~ PB ~ IL7R^`+` ~ T[M]" = Controls_1_PB_radarplots[[3]][[3]],
                                 "Resting ~ PB ~ T[TE]" = Controls_1_PB_radarplots[[4]][[1]],
                                 "Resting ~ PB ~ T[EM]" = Controls_1_PB_radarplots[[5]][[2]],
                                 "Resting ~ PB ~ IL7R^`+` ~ T[M]" = Controls_1_PB_radarplots[[6]][[3]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_1_PB_radarplots[[i]] + theme(text = element_text(size = 12)) + 
           ggtitle(label = "PB", subtitle = parse(text = names(Controls_1_PB_radarplots[i]))), 
         filename = paste0("Output/Figures/Radarplots/Controls_1/Final/PB_10x_atlas/Genes_defining_state/",
                           names(Controls_1_PB_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_1_PB_radarplots <- list()

# Run a loop plotting genes of interest in the Controls_1 object against the PB_10_x_atlas
for (i in names(Controls_query_list[["Controls_1"]][1])) {
  for (j in names(genes_to_plot)) {
    temp <- 
      (plot.states.radar(ref = PB_10x_atlas, 
                         query = list("Activated PB" = Controls_query_list[["Controls_1"]][[6]],
                                      "Resting PB" = Controls_query_list[["Controls_1"]][[7]]),
                         genes4radar = genes_to_plot[[`j`]],
                         cols = c("purple", "green"),  
                         return.as.list = TRUE))
    Controls_1_PB_radarplots[[length(Controls_1_PB_radarplots) +1]] <- temp
  }
}

# Rename the objects in the list - [Sample][State]
for (i in 1:8) {
  names(Controls_1_PB_radarplots[[i]]) <- cluster_names_PB_10x_atlas
}
names(Controls_1_PB_radarplots) <- c("Exhaustion", "Activation", "Effector", "Regulators", 
                                     "Residency", "Homing_1", "Homing_2", "Homing_3")

# Export the images
for(i in names(Controls_1_PB_radarplots)){
  for (j in names(Controls_1_PB_radarplots[[1]])) {
    ggsave(plot = 
             Controls_1_PB_radarplots[[i]][[j]] + theme(text = element_text(size = 12)) +   
             ggtitle(label = "PB", subtitle = parse(text = names(Controls_1_PB_radarplots[[1]][j]))),
           filename = paste0("Output/Figures/Radarplots/Controls_1/Final/PB_10x_atlas/Genes_of_interest/",
                             names(Controls_1_PB_radarplots[i]), "/PB_10x_atlas_", j, ".png"), 
           units = "cm", height = 16, width = 18)
  }
}


# 4.2B Plot Controls 2 - BM atlas radar plots ------------------------------------------------------

# Create an empty list to store results
Controls_2_radarplots <- list()

# Run a loop plotting all combinations of the Controls_2 object against the BM_10_x_atlas
for (i in names(Controls_query_list[["Controls_2"]][1])) {
  for (j in 1:length(top_10_FindAllMarkers_BM_10x_atlas)) {
    temp <- 
      (plot.states.radar(ref = BM_10x_atlas, 
                         query = list("HD" = Controls_query_list[["Controls_2"]][[2]],
                                      "MGUS" = Controls_query_list[["Controls_2"]][[3]],
                                      "SMM" = Controls_query_list[["Controls_2"]][[4]],
                                      "MM" = Controls_query_list[["Controls_2"]][[5]]), 
                         genes4radar = top_10_FindAllMarkers_BM_10x_atlas[[j]]$gene, 
                         cols = c("green", "pink", "purple", "orange"), 
                         return.as.list = TRUE))
    Controls_2_radarplots[[length(Controls_2_radarplots) +1]] <- temp
  }
}

# Extract the images we're interested in
Controls_2_radarplots <- list("T[EM]" = Controls_2_radarplots[[1]][[1]],
                              "T[N]" = Controls_2_radarplots[[2]][[2]],
                              "IL7R^`+` ~ T[M]" = Controls_2_radarplots[[3]][[3]],
                              "T[TE]" = Controls_2_radarplots[[4]][[4]],
                              "Activated" = Controls_2_radarplots[[5]][[5]],
                              "Cyto ~ T[EM]" = Controls_2_radarplots[[6]][[6]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_2_radarplots[[i]] + theme(text = element_text(size = 12)) +  
           ggtitle(label = "BM", subtitle = parse(text = names(Controls_2_radarplots[i]))),
         filename = paste0("Output/Figures/Radarplots/Controls_2/Final/BM_10x_atlas/Genes_defining_cluster/",
                           names(Controls_2_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_2_radarplots <- list()

# Run a loop plotting all combinations of the Controls_2 object against the BM_10_x_atlas
for (i in names(top_10_FindMarkers_Controls_2_BM_10x_atlas)) {
  temp <- 
    (plot.states.radar(ref = BM_10x_atlas, 
                       query = list("HD" = Controls_query_list[["Controls_2"]][[2]],
                                    "MGUS" = Controls_query_list[["Controls_2"]][[3]],
                                    "SMM" = Controls_query_list[["Controls_2"]][[4]],
                                    "MM" = Controls_query_list[["Controls_2"]][[5]]),
                       genes4radar = top_10_FindMarkers_Controls_2_BM_10x_atlas[[`i`]]$gene,
                       cols = c("green", "pink", "purple", "orange"),  
                       return.as.list = TRUE))
  Controls_2_radarplots[[length(Controls_2_radarplots) +1]] <- temp
}

# Extract the images we're interested in - [Sample][State]
Controls_2_radarplots <- list("Age-matched ~ Controls ~ T[EM]" = Controls_2_radarplots[[3]][[1]],
                              "Age-matched ~ Controls ~ IL7R^`+` ~ T[M]" = Controls_2_radarplots[[4]][[3]],
                              "MGUS ~ T[EM]" = Controls_2_radarplots[[5]][[1]],
                              "MGUS ~ IL7R^`+` ~ T[M]" = Controls_2_radarplots[[6]][[3]],
                              "SMM ~ T[EM]" = Controls_2_radarplots[[7]][[1]],
                              "SMM ~ IL7R^`+` ~ T[M]" = Controls_2_radarplots[[8]][[3]],
                              "MM ~ T[EM]" = Controls_2_radarplots[[9]][[1]],
                              "MM ~ IL7R^`+` ~ T[M]" = Controls_2_radarplots[[10]][[3]])

# Export the images
for(i in 1:8){
  ggsave(plot = 
           Controls_2_radarplots[[i]] + theme(text = element_text(size = 12)) + 
           ggtitle(label = "BM", subtitle = parse(text = names(Controls_2_radarplots[i]))), 
         filename = paste0("Output/Figures/Radarplots/Controls_2/Final/BM_10x_atlas/Genes_defining_state/",
                           names(Controls_2_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_2_radarplots <- list()

# Run a loop plotting genes of interest in the Controls_2 object against the BM_10_x_atlas
for (i in names(Controls_query_list[["Controls_2"]][1])) {
  for (j in names(genes_to_plot)) {
    temp <- 
      (plot.states.radar(ref = BM_10x_atlas, 
                         query = list("Age-matched BM" = Controls_query_list[["Controls_2"]][[2]],
                                      "MGUS BM" = Controls_query_list[["Controls_2"]][[3]],
                                      "SMM BM" = Controls_query_list[["Controls_2"]][[4]]),
                         genes4radar = genes_to_plot[[`j`]],
                         cols = c("green", "pink", "orange"),  
                         return.as.list = TRUE, min.cells = 2))
    Controls_2_radarplots[[length(Controls_2_radarplots) +1]] <- temp
  }
}

# Rename the objects in the list - [Sample][State]
for (i in 1:8) {
  names(Controls_2_radarplots[[i]]) <- cluster_names_BM_10x_atlas
}
names(Controls_2_radarplots) <- c("Exhaustion", "Activation", "Effector", "Regulators", 
                                  "Residency", "Homing_1", "Homing_2", "Homing_3")

# Export the images
for(i in names(Controls_2_radarplots)){
  for (j in names(Controls_2_radarplots[[1]])) {
    ggsave(plot = 
             Controls_2_radarplots[[i]][[j]] + theme(text = element_text(size = 12)) + 
             ggtitle(label = "BM", subtitle = parse(text = names(Controls_2_radarplots[[1]][j]))), 
           filename = paste0("Output/Figures/Radarplots/Controls_2/Final/BM_10x_atlas/Genes_of_interest/",
                             names(Controls_2_radarplots[i]), "/BM_10x_atlas_", j, ".png"), 
           units = "cm", height = 16, width = 18)
  }
}


# 4.2C Plot Controls 3 - PB atlas radar plots ------------------------------------------------------

# Create an empty list to store results
Controls_3_radarplots <- list()

# Run a loop plotting all combinations of the Controls_3 object against the PB_10_x_atlas
for (i in names(Controls_query_list[["Controls_3"]][1])) {
  for (j in 1:length(top_10_FindAllMarkers_PB_10x_atlas)) {
    temp <- 
      (plot.states.radar(ref = PB_10x_atlas, 
                         query = list("Old PB" = Controls_query_list[["Controls_3"]][[2]],
                                      "Young PB" = Controls_query_list[["Controls_3"]][[3]]), 
                         genes4radar = top_10_FindAllMarkers_PB_10x_atlas[[j]]$gene, 
                         cols = c("purple", "green"), 
                         return.as.list = TRUE))
    Controls_3_radarplots[[length(Controls_3_radarplots) +1]] <- temp
  }
}

# Extract the images we're interested in
Controls_3_radarplots <- list("T[TE]" = Controls_3_radarplots[[1]][[1]],
                              "T[EM]" = Controls_3_radarplots[[2]][[2]],
                              "IL7R^`+` ~ T[M]" = Controls_3_radarplots[[3]][[3]],
                              "T[N]" = Controls_3_radarplots[[4]][[4]],
                              "Activated" = Controls_3_radarplots[[5]][[5]],
                              "IFN" = Controls_3_radarplots[[6]][[6]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_3_radarplots[[i]] + theme(text = element_text(size = 12)) + 
           ggtitle(label = "PB", subtitle = parse(text = names(Controls_3_radarplots[i]))), 
         filename = paste0("Output/Figures/Radarplots/Controls_3/Final/PB_10x_atlas/Genes_defining_cluster/",
                           names(Controls_3_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_3_radarplots <- list()

# Run a loop plotting all combinations of the Controls_3 object against the PB_10_x_atlas
for (i in names(top_10_FindMarkers_Controls_3_PB_10x_atlas)) {
  temp <- 
    (plot.states.radar(ref = PB_10x_atlas, 
                       query = list("Old PB" = Controls_query_list[["Controls_3"]][[2]],
                                    "Young PB" = Controls_query_list[["Controls_3"]][[3]]),
                       genes4radar = top_10_FindMarkers_Controls_3_PB_10x_atlas[[`i`]]$gene,
                       cols = c("purple", "green"),  
                       return.as.list = TRUE))
  Controls_3_radarplots[[length(Controls_3_radarplots) +1]] <- temp
}

# Extract the images we're interested in - [Sample][State]
Controls_3_radarplots <- list("Age-matched ~ Controls ~ T[TE]" = Controls_3_radarplots[[4]][[1]],
                              "Age-matched ~ Controls ~ T[EM]" = Controls_3_radarplots[[5]][[2]],
                              "Age-matched ~ Controls ~ T[N]" = Controls_3_radarplots[[6]][[4]],
                              "Young ~ T[TE]" = Controls_3_radarplots[[7]][[1]],
                              "Young ~ T[EM]" = Controls_3_radarplots[[8]][[2]],
                              "Young ~ T[N]" = Controls_3_radarplots[[9]][[4]])

# Export the images
for(i in 1:6){
  ggsave(plot = 
           Controls_3_radarplots[[i]] + theme(text = element_text(size = 12)) + 
           ggtitle(label = "PB", subtitle = parse(text = names(Controls_3_radarplots[i]))), 
         filename = paste0("Output/Figures/Radarplots/Controls_3/Final/PB_10x_atlas/Genes_defining_state/",
                           names(Controls_3_radarplots[i]), ".png"), 
         units = "cm", height = 16, width = 18)
}

# Empty the list
Controls_3_radarplots <- list()

# Run a loop plotting genes of interest in the Controls_3 object against the PB_10_x_atlas
for (i in names(Controls_query_list[["Controls_3"]][1])) {
  for (j in names(genes_to_plot)) {
    temp <- 
      (plot.states.radar(ref = PB_10x_atlas, 
                         query = list("Old PB" = Controls_query_list[["Controls_3"]][[2]],
                                      "Young PB" = Controls_query_list[["Controls_3"]][[3]]),
                         genes4radar = genes_to_plot[[`j`]],
                         cols = c("purple", "green"),  
                         return.as.list = TRUE))
    Controls_3_radarplots[[length(Controls_3_radarplots) +1]] <- temp
  }
}

# Rename the objects in the list - [Sample][State]
for (i in 1:8) {
  names(Controls_3_radarplots[[i]]) <- cluster_names_PB_10x_atlas
}
names(Controls_3_radarplots) <- c("Exhaustion", "Activation", "Effector", "Regulators", 
                                  "Residency", "Homing_1", "Homing_2", "Homing_3")

# Export the images
for(i in names(Controls_3_radarplots)){
  for (j in names(Controls_3_radarplots[[1]])) {
    ggsave(plot = 
             Controls_3_radarplots[[i]][[j]] + theme(text = element_text(size = 12)) +  
             ggtitle(label = "PB", subtitle = parse(text = names(Controls_3_radarplots[[1]][j]))), 
           filename = paste0("Output/Figures/Radarplots/Controls_3/Final/PB_10x_atlas/Genes_of_interest/",
                             names(Controls_3_radarplots[i]), "/PB_10x_atlas_", j, ".png"), 
           units = "cm", height = 16, width = 18)
  }
}


# 5. Save out data and print stats -----------------------------------------------------------------

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '8. Further analysis of projected Controls data' script")