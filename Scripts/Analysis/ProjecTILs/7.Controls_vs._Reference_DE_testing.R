# Script information -------------------------------------------------------------------------------

# Title: Analysis of projected Controls data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will analyse the ProjectTILs projection of public data onto our
# custom reference maps. We will first generate dimplots illustrating overlay on our reference UMAP 
# and bar plots of cluster distribution. We will then perform DE testing which will feed the next 
# script to create additional figures. 
# NB: find.discriminant.genes will fail if clusters are not observable in both reference and query. 
# We will therefore only be able to test cluster vs. cluster with sufficient overlap.
# NB: find.discriminant.genes uses assay <- "RNA" by default, no need to set it manually.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '7. DE testing of projected data' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire, openxlsx, multtest, metap, EnhancedVolcano, ProjecTILs, ggpubr, scales)
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

# Load additional resources
black_list <- readRDS("Data/R_out/black_list.rds")
annotations <- read.csv("Data/Public_data/annotations.csv")
load("Data/R_out/10x/DE_Analysis_10x.rds")

# Modify black_list to retain TCR genes
black_list_TCR <- 
  black_list[!grepl("TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC", black_list)]

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

# Create a character vector of genes to plot
top_10_markers_all_10x <- 
  FindAllMarkers_All_list$all_10x %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique_all_10x <- unique(top_10_markers_all_10x$gene)


# 3. Pre-processing --------------------------------------------------------------------------------

# Move atlases into a list
all_atlas_list <- list("all_10x_atlas" = all_10x_atlas,
                       "BM_10x_atlas" = BM_10x_atlas,
                       "PB_10x_atlas" = PB_10x_atlas)

# Remove non helpful genes for DE analysis #NB: This will remove genes from the object
for (i in names (all_atlas_list)) {
  counts <- GetAssayData(all_atlas_list[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  all_atlas_list[[i]] <- subset(all_atlas_list[[i]], features = rownames(counts))
}

for (i in names (Controls_query_list[["Controls_1"]])) {
  counts <- GetAssayData(Controls_query_list[["Controls_1"]][[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  Controls_query_list[["Controls_1"]][[i]] <- subset(Controls_query_list[["Controls_1"]][[i]], 
                                                     features = rownames(counts))
}

for (i in names (Controls_query_list[["Controls_2"]])) {
  counts <- GetAssayData(Controls_query_list[["Controls_2"]][[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  Controls_query_list[["Controls_2"]][[i]] <- subset(Controls_query_list[["Controls_2"]][[i]], 
                                                     features = rownames(counts))
}

for (i in names (Controls_query_list[["Controls_3"]])) {
  counts <- GetAssayData(Controls_query_list[["Controls_3"]][[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  Controls_query_list[["Controls_3"]][[i]] <- subset(Controls_query_list[["Controls_3"]][[i]], 
                                                     features = rownames(counts))
}


# 4. Output ProjectTILs figures and data -----------------------------------------------------------

# 4.1 Plot projections -----------------------------------------------------------------------------

# 4.1.A.i Plot Controls 1 - all atlas projections --------------------------------------------------

# Plot projections - Controls 1 - All atlas
for (i in 1:length(Controls_query_list$Controls_1[1:3])) {
  temp <- 
    (plot.projection(ref = all_10x_atlas, 
                     query = Controls_query_list$Controls_1[[i]],
                     cols = all_10x_atlas_plot_cols)) 
  temp <- temp[[1]]
  temp <- temp + theme_bw() + scale_color_manual(name = "Cluster", 
                                    values = all_10x_atlas_plot_cols, labels = label_parse()) + 
    ggtitle(Controls_1_titles[i]) +
    theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
          legend.title = element_text(size = 18), legend.text.align = 0)
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Dimplots/Controls_1/Final/Projection/",
                                   names(Controls_query_list$Controls_1)[i], ".png"), 
         units = "cm", height = 20, width = 25, bg = "transparent")
}


# 4.1.A.ii Plot Controls 1 - BM atlas projections --------------------------------------------------

# Plot projections - Controls 1 - BM atlas
for (i in 4:(length(Controls_query_list$Controls_1)-2)) {
  temp <- 
    (plot.projection(ref = BM_10x_atlas, 
                     query = Controls_query_list$Controls_1[[`i`]],
                     cols = BM_10x_atlas_plot_cols))
  temp <- temp[[1]]
  temp <- temp + theme_bw() + scale_color_manual(name = "Cluster", 
                                    values = BM_10x_atlas_plot_cols, labels = label_parse()) + 
    ggtitle(Controls_1_titles[i]) +
    theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
          legend.title = element_text(size = 18), legend.text.align = 0)
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Dimplots/Controls_1/Final/Projection/",
                                   names(Controls_query_list$Controls_1)[i], ".png"), 
         units = "cm", height = 20, width = 25, bg = "transparent")
}

# 4.1.A.iii Plot Controls 1 - PB atlas projections -------------------------------------------------

# Plot projections - Controls 1 - PB atlas
for (i in 6:length(Controls_query_list$Controls_1)) {
  temp <- 
    (plot.projection(ref = PB_10x_atlas, 
                     query = Controls_query_list$Controls_1[[`i`]],
                     cols = PB_10x_atlas_plot_cols))
  temp <- temp[[1]]
  temp <- temp + theme_bw() + scale_color_manual(name = "Cluster", 
                                    values = PB_10x_atlas_plot_cols, labels = label_parse()) + 
    ggtitle(Controls_1_titles[i]) +
    theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
          legend.title = element_text(size = 18), legend.text.align = 0)
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Dimplots/Controls_1/Final/Projection/",
                                   names(Controls_query_list$Controls_1)[i], ".png"), 
         units = "cm", height = 20, width = 25, bg = "transparent")
}

# 4.1.B.i Plot Controls 2 - BM atlas projections ---------------------------------------------------

for (i in 1:length(Controls_query_list$Controls_2)) {
  temp <- 
    (plot.projection(ref = BM_10x_atlas, 
                     query = Controls_query_list$Controls_2[[`i`]],
                     cols = ccolours_BM_10x_atlas))
  temp <- temp[[1]]
  temp <- temp + theme_bw() + scale_color_manual(name = "Cluster", 
                                    values = BM_10x_atlas_plot_cols, labels = label_parse()) + 
    ggtitle(Controls_2_titles[i]) +
    theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
          legend.title = element_text(size = 18), legend.text.align = 0)
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Dimplots/Controls_2/Final/Projection/",
                                   names(Controls_query_list$Controls_2)[i], ".png"), 
         units = "cm", height = 20, width = 25, bg = "transparent")
}

# 4.1.C.i Plot Controls 3 - PB atlas projections ---------------------------------------------------

for (i in 1:length(Controls_query_list$Controls_3)) {
  temp <- 
    (plot.projection(ref = PB_10x_atlas, 
                     query = Controls_query_list$Controls_3[[`i`]],
                     cols = ccolours_PB_10x_atlas))
  temp <- temp[[1]]
  temp <- temp + theme_bw() + scale_color_manual(name = "Cluster", 
                                    values = PB_10x_atlas_plot_cols, labels = label_parse()) + 
    ggtitle(Controls_3_titles[i]) +
    theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
          legend.title = element_text(size = 18), legend.text.align = 0)
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Dimplots/Controls_3/Final/Projection/",
                                   names(Controls_query_list$Controls_3)[i], ".png"), 
         units = "cm", height = 20, width = 25, bg = "transparent")
}


# 4.2 Create bar plots of cluster distribution -----------------------------------------------------

# 4.2.A.i Stacked barplots of cluster distribution - Controls #1 vs. all_10x_atlas -----------------

# Gather data to plot
temp_1 <- as.data.frame(table(all_10x_atlas$functional.cluster, 
                              all_10x_atlas$Sample))
temp_2 <- as.data.frame(table(Controls_query_list$Controls_1$Activated$functional.cluster, 
                              Controls_query_list$Controls_1$Activated$Sample))
temp_3 <- as.data.frame(table(Controls_query_list$Controls_1$Resting$functional.cluster, 
                              Controls_query_list$Controls_1$Resting$Sample))

# rbind all data to allow plotting together
temp <- rbind(temp_1, temp_2, temp_3)

# Convert cell numbers to proportions based on disease
temp <- temp %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Fix the order of plotting from largest to smallest
temp <- temp %>% arrange(desc(Proportion)) %>% arrange(Var2)

# Plot cluster distribution of the Controls 1 object projected on the all_10x atlas ----------------
ggplot(temp, 
       aes(fill = fct_rev(Var1), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_all_10x_atlas), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_all_10x_atlas))) +
  labs(x = "Sample", y = "Proportion (%)", 
       title = "Age-matched controls projected on Combined reference atlas") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_projected_on_all_10x_atlas", ".png"), 
       units = "cm", height = 20, width = 18, bg = "transparent")

# 4.2.A.ii Stacked barplots of cluster distribution - Controls #1 vs. BM_10x_atlas -----------------

# Gather data to plot
temp_1 <- as.data.frame(table(BM_10x_atlas$functional.cluster, 
                              BM_10x_atlas$Sample))
temp_3 <- as.data.frame(table(Controls_query_list$Controls_1$`Activated BM`$functional.cluster, 
                              Controls_query_list$Controls_1$`Activated BM`$Sample))
temp_2 <- as.data.frame(table(Controls_query_list$Controls_1$`Resting BM`$functional.cluster, 
                              Controls_query_list$Controls_1$`Resting BM`$Sample))

# rbind all data to allow plotting together
temp <- rbind(temp_1, temp_2, temp_3)

# Convert cell numbers to proportions based on disease
temp <- temp %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Fix the order of plotting from largest to smallest
temp <- temp %>% arrange(desc(Proportion)) %>% arrange(Var2)

# Plot cluster distribution of the Controls 1 object projected on the BM_10x atlas
ggplot(temp, 
       aes(fill = fct_rev(Var1), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = rev(ccolours_BM_10x_atlas), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_BM_10x_atlas))) +
  labs(x = "Sample", y = "Proportion (%)", 
       title = "Age-matched controls projected on MM BM reference atlas") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_projected_on_BM_10x_atlas", ".png"), 
       units = "cm", height = 20, width = 18, bg = "transparent")

# 4.2.A.iii Stacked barplots of cluster distribution - Controls #1 vs. PB_10x_atlas ----------------

# Gather data to plot
temp_1 <- as.data.frame(table(PB_10x_atlas$functional.cluster, 
                              PB_10x_atlas$Sample))
temp_3 <- as.data.frame(table(Controls_query_list$Controls_1$`Activated PB`$functional.cluster, 
                              Controls_query_list$Controls_1$`Activated PB`$Sample))
temp_2 <- as.data.frame(table(Controls_query_list$Controls_1$`Resting PB`$functional.cluster, 
                              Controls_query_list$Controls_1$`Resting PB`$Sample))

# rbind all data to allow plotting together
temp <- rbind(temp_1, temp_2, temp_3)

# Convert cell numbers to proportions based on disease
temp <- temp %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Fix the order of plotting from largest to smallest
temp <- temp %>% arrange(desc(Proportion)) %>% arrange(Var2)

# Plot cluster distribution of the Controls 1 object projected on the PB_10x atlas
ggplot(temp, 
       aes(fill = fct_rev(Var1), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_PB_10x_atlas), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_PB_10x_atlas))) +
  labs(x = "Sample", y = "Proportion (%)", 
       title = "Age-matched controls projected on MM PB reference atlas") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_1/Final/Barplot_",
                "Controls_1_projected_on_PB_10x_atlas", ".png"), 
       units = "cm", height = 20, width = 18, bg = "transparent")

# 4.2.B.i Stacked barplots of cluster distribution - Controls #2 vs. BM_10x_atlas ------------------

# Gather data to plot
temp_1 <- as.data.frame(table(BM_10x_atlas$functional.cluster, 
                              BM_10x_atlas$Sample))
temp_2 <- as.data.frame(table(Controls_query_list$Controls_2$`Age-matched Controls BM`$functional.cluster, 
                              Controls_query_list$Controls_2$`Age-matched Controls BM`$Sample))
temp_3 <- as.data.frame(table(Controls_query_list$Controls_2$`MGUS BM`$functional.cluster, 
                              Controls_query_list$Controls_2$`MGUS BM`$Sample))
temp_4 <- as.data.frame(table(Controls_query_list$Controls_2$`SMM BM`$functional.cluster, 
                              Controls_query_list$Controls_2$`SMM BM`$Sample))
temp_5 <- as.data.frame(table(Controls_query_list$Controls_2$`MM BM`$functional.cluster, 
                              Controls_query_list$Controls_2$`MM BM`$Sample))

# rbind all data to allow plotting together
temp <- rbind(temp_1, temp_2, temp_3, temp_4, temp_5)

# Convert cell numbers to proportions based on disease
temp <- temp %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Fix the order of plotting from largest to smallest
temp <- temp %>% arrange(desc(Proportion)) %>% arrange(Var2)

# Plot cluster distribution of the Controls 2 object projected on the BM_10x atlas
ggplot(temp, 
       aes(fill = fct_rev(Var1), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_BM_10x_atlas), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_BM_10x_atlas))) +
  labs(x = "Sample", y = "Proportion (%)", 
       title = "Controls_2 projected on MM BM reference atlas") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_2/Final/Barplot_",
                "Controls_2_projected_on_BM_10x_atlas", ".png"), 
       units = "cm", height = 20, width = 18, bg = "transparent")

# 4.2.C.i Stacked barplots of cluster distribution - Controls #3 vs. PB_10x_atlas ------------------

# Gather data to plot
temp_1 <- as.data.frame(table(PB_10x_atlas$functional.cluster, 
                              PB_10x_atlas$Sample))
temp_2 <- as.data.frame(table(Controls_query_list$Controls_3$`Age-matched Controls PB`$functional.cluster, 
                              Controls_query_list$Controls_3$`Age-matched Controls PB`$Sample))
temp_3 <- as.data.frame(table(Controls_query_list$Controls_3$`Young PB`$functional.cluster, 
                              Controls_query_list$Controls_3$`Young PB`$Sample))

# rbind all data to allow plotting together
temp <- rbind(temp_1, temp_2, temp_3)

# Convert cell numbers to proportions based on disease
temp <- temp %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Fix the order of plotting from largest to smallest
temp <- temp %>% arrange(desc(Proportion)) %>% arrange(Var2)

# Plot cluster distribution of the Controls 3 object projected on the PB_10x atlas
ggplot(temp, 
       aes(fill = fct_rev(Var1), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(ccolours_PB_10x_atlas), name = "Cluster", 
                    labels = parse(text = rev(cluster_names_PB_10x_atlas))) +
  labs(x = "Sample", y = "Proportion (%)", 
       title = "Controls_3 projected on MM PB reference atlas") + theme_bw() + 
  theme(legend.text.align = 0) + theme(text=element_text(size = 18))
ggsave(filename = 
         paste0("Output/Figures/Barplots/Controls_3/Final/Barplot_",
                "Controls_3_projected_on_PB_10x_atlas", ".png"), 
       units = "cm", height = 20, width = 18, bg = "transparent")


# 5. Perform DE and output Volcano plots -----------------------------------------------------------

# NB: Can ammend labels with argument, selectLab = top_10_unique_all_10x if required

# 5.1.A DE testing - FindAllMarkers - Compile all genes that broadly segregate data ----------------

# 5.1.A.i DE testing - FindAllMarkers - Controls #1 vs. all_10x_atlas ------------------------------

# Create a list to store DE results
Controls_1_vs_all_10x_atlas_FindAllMarkers_All <- list()

# Iterate over list - FindAllMarkers
for (i in 1:length(Controls_query_list[["Controls_1"]][1:3])) {
  temp <- find.discriminant.genes(all_atlas_list$all_10x_atlas, 
                                  Controls_query_list[["Controls_1"]][[i]], state= "all")
  Controls_1_vs_all_10x_atlas_FindAllMarkers_All[[i]] <- temp %>% 
    rownames_to_column() %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
}

# Rename items in the list
names(Controls_1_vs_all_10x_atlas_FindAllMarkers_All) <- 
  names((Controls_query_list[["Controls_1"]][1:3]))

# Save out the DE results
write.xlsx(Controls_1_vs_all_10x_atlas_FindAllMarkers_All,
           "Output/DE/Controls_1/Final/all_10x_atlas/Controls_1_vs_all_10x_atlas_FindAllMarkers_All.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(Controls_1_vs_all_10x_atlas_FindAllMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_1_vs_all_10x_atlas_FindAllMarkers_All[[`i`]], 
                     lab = Controls_1_vs_all_10x_atlas_FindAllMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(Controls_1_vs_all_10x_atlas_FindAllMarkers_All)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_1/Final/all_10x_atlas/",
                  names(Controls_1_vs_all_10x_atlas_FindAllMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}

# 5.1.A.ii DE testing - FindAllMarkers - Controls #1 vs. BM_10x_atlas ------------------------------

# Create a list to store DE results
Controls_1_vs_BM_10x_atlas_FindAllMarkers_All <- list()

# Iterate over list - FindAllMarkers
for (i in (names(Controls_query_list[["Controls_1"]][4:5]))) {
  temp <- find.discriminant.genes(all_atlas_list$BM_10x_atlas, 
                                  Controls_query_list[["Controls_1"]][[i]], state= "all")
  Controls_1_vs_BM_10x_atlas_FindAllMarkers_All[[i]] <- temp %>% 
    rownames_to_column() %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
}

# Rename items in the list
names(Controls_1_vs_BM_10x_atlas_FindAllMarkers_All) <- 
  names((Controls_query_list[["Controls_1"]][4:5]))

# Save out the DE results
write.xlsx(Controls_1_vs_BM_10x_atlas_FindAllMarkers_All,
           "Output/DE/Controls_1/Final/BM_10x_atlas/Controls_1_vs_BM_10x_atlas_FindAllMarkers_All.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(Controls_1_vs_BM_10x_atlas_FindAllMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_1_vs_BM_10x_atlas_FindAllMarkers_All[[`i`]], 
                     lab = Controls_1_vs_BM_10x_atlas_FindAllMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(Controls_1_vs_BM_10x_atlas_FindAllMarkers_All)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_1/Final/BM_10x_atlas/",
                  names(Controls_1_vs_BM_10x_atlas_FindAllMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}

# 5.1.A.iii DE testing - FindAllMarkers - Controls #1 vs. PB_10x_atlas -----------------------------

# Create a list to store DE results
Controls_1_vs_PB_10x_atlas_FindAllMarkers_All <- list()

# Iterate over list - FindAllMarkers
for (i in (names(Controls_query_list[["Controls_1"]][6:7]))) {
  temp <- find.discriminant.genes(all_atlas_list$PB_10x_atlas, 
                                  Controls_query_list[["Controls_1"]][[i]], state= "all")
  Controls_1_vs_PB_10x_atlas_FindAllMarkers_All[[i]] <- temp %>% 
    rownames_to_column() %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
}

# Rename items in the list
names(Controls_1_vs_PB_10x_atlas_FindAllMarkers_All) <- 
  names((Controls_query_list[["Controls_1"]][6:7]))

# Save out the DE results
write.xlsx(Controls_1_vs_PB_10x_atlas_FindAllMarkers_All,
           "Output/DE/Controls_1/Final/PB_10x_atlas/Controls_1_vs_PB_10x_atlas_FindAllMarkers_All.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(Controls_1_vs_PB_10x_atlas_FindAllMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_1_vs_PB_10x_atlas_FindAllMarkers_All[[`i`]], 
                     lab = Controls_1_vs_PB_10x_atlas_FindAllMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(Controls_1_vs_PB_10x_atlas_FindAllMarkers_All)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_1/Final/PB_10x_atlas/",
                  names(Controls_1_vs_PB_10x_atlas_FindAllMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}

# 5.1.B.i DE testing - FindAllMarkers - Controls #2 vs. BM_10x_atlas -------------------------------

# Create a list to store DE results
Controls_2_vs_BM_10x_atlas_FindAllMarkers_All <- list()

# Iterate over list - FindAllMarkers
for (i in 1:length(Controls_query_list[["Controls_2"]])) {
  temp <- find.discriminant.genes(all_atlas_list$BM_10x_atlas, 
                                  Controls_query_list[["Controls_2"]][[i]], state= "all")
  Controls_2_vs_BM_10x_atlas_FindAllMarkers_All[[i]] <- temp %>% 
    rownames_to_column() %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
}

# Rename items in the list
names(Controls_2_vs_BM_10x_atlas_FindAllMarkers_All) <- 
  names((Controls_query_list[["Controls_2"]]))

# Save out the DE results
write.xlsx(Controls_2_vs_BM_10x_atlas_FindAllMarkers_All,
           "Output/DE/Controls_2/Final/BM_10x_atlas/Controls_2_vs_BM_10x_atlas_FindAllMarkers_All.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(Controls_2_vs_BM_10x_atlas_FindAllMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_2_vs_BM_10x_atlas_FindAllMarkers_All[[`i`]], 
                     lab = Controls_2_vs_BM_10x_atlas_FindAllMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(Controls_2_vs_BM_10x_atlas_FindAllMarkers_All)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_2/Final/BM_10x_atlas/",
                  names(Controls_2_vs_BM_10x_atlas_FindAllMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}

# 5.1.C.i DE testing - FindAllMarkers - Controls #3 vs. PB_10x_atlas -------------------------------

# Create a list to store DE results
Controls_3_vs_PB_10x_atlas_FindAllMarkers_All <- list()

# Iterate over list - FindAllMarkers
for (i in 1:length(Controls_query_list[["Controls_3"]])) {
  temp <- find.discriminant.genes(all_atlas_list$PB_10x_atlas, 
                                  Controls_query_list[["Controls_3"]][[i]], state= "all")
  Controls_3_vs_PB_10x_atlas_FindAllMarkers_All[[i]] <- temp %>% 
    rownames_to_column() %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
}

# Rename items in the list
names(Controls_3_vs_PB_10x_atlas_FindAllMarkers_All) <- 
  names((Controls_query_list[["Controls_3"]]))

# Save out the DE results
write.xlsx(Controls_3_vs_PB_10x_atlas_FindAllMarkers_All,
           "Output/DE/Controls_3/Final/PB_10x_atlas/Controls_3_vs_PB_10x_atlas_FindAllMarkers_All.xlsx")

# Output Volcano plots for all tests
for (i in 1:length(Controls_3_vs_PB_10x_atlas_FindAllMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_3_vs_PB_10x_atlas_FindAllMarkers_All[[`i`]], 
                     lab = Controls_3_vs_PB_10x_atlas_FindAllMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(Controls_3_vs_PB_10x_atlas_FindAllMarkers_All)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_3/Final/PB_10x_atlas/",
                  names(Controls_3_vs_PB_10x_atlas_FindAllMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}


# 5.2 DE testing - FindMarkers - Compile all genes that segregate clusters -------------------------

# 5.2.A.i DE testing - FindMarkers - Controls #1 vs. all_10x_atlas ---------------------------------

# Create a list to store DE results
Controls_1_vs_all_10x_atlas_FindMarkers_All <- list()

# Iterate over list to test TEM, TTE and TN
for (i in (names(Controls_query_list[["Controls_1"]][1:3]))) {
  for (j in cluster_names_all_10x_atlas[1:3]) {
    temp <- find.discriminant.genes(all_atlas_list$all_10x_atlas, 
                                    Controls_query_list$Controls_1[[i]], state = j)
    Controls_1_vs_all_10x_atlas_FindMarkers_All[[i]][[j]] <- temp %>% 
      rownames_to_column() %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
  }
}

# Flatten the list to allow export
Controls_1_vs_all_10x_atlas_FindMarkers_All <- 
  unlist(Controls_1_vs_all_10x_atlas_FindMarkers_All, recursive = FALSE)

# Copy list and fix names to allow export to excel without error
temp <- Controls_1_vs_all_10x_atlas_FindMarkers_All
names(temp) <- gsub("([.])", " ", names(temp))
names(temp) <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "", names(temp))

# Save out the DE results
write.xlsx(temp,
           "Output/DE/Controls_1/Final/all_10x_atlas/Controls_1_vs_all_10x_atlas_FindMarkers_All.xlsx")

# Fix names for items in list for pretty Volcano plots
names(Controls_1_vs_all_10x_atlas_FindMarkers_All) <- 
  gsub("Controls 1 complete", "Complete", names(Controls_1_vs_all_10x_atlas_FindMarkers_All))
names(Controls_1_vs_all_10x_atlas_FindMarkers_All) <- 
  gsub("[//.]", " ~ ", names(Controls_1_vs_all_10x_atlas_FindMarkers_All))

# Output Volcano plots for all tests
for (i in 1:length(Controls_1_vs_all_10x_atlas_FindMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_1_vs_all_10x_atlas_FindMarkers_All[[`i`]], 
                     lab = Controls_1_vs_all_10x_atlas_FindMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = parse(text = names(Controls_1_vs_all_10x_atlas_FindMarkers_All)[i]), 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_1/Final/all_10x_atlas/",
                  names(Controls_1_vs_all_10x_atlas_FindMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}

# 5.2.A.ii DE testing - FindMarkers - Controls #1 vs. BM_10x_atlas ---------------------------------

# Create a list to store DE results
Controls_1_vs_BM_10x_atlas_FindMarkers_All <- list()

# Iterate over list to test TEM, IL7R_TM and TTE
for (i in (names(Controls_query_list[["Controls_1"]][4:5]))) {
  for (j in cluster_names_BM_10x_atlas[c(1,3,4)]) {
    temp <- find.discriminant.genes(all_atlas_list$BM_10x_atlas, 
                                    Controls_query_list$Controls_1[[i]], state = j)
    Controls_1_vs_BM_10x_atlas_FindMarkers_All[[i]][[j]] <- temp %>% 
      rownames_to_column() %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
  }
}

# Flatten the list to allow export and fix names
Controls_1_vs_BM_10x_atlas_FindMarkers_All <- 
  unlist(Controls_1_vs_BM_10x_atlas_FindMarkers_All, recursive = FALSE)

# Copy list and fix names to allow export to excel without error
temp <- Controls_1_vs_BM_10x_atlas_FindMarkers_All
names(temp) <- gsub("([.])", " ", names(temp))
names(temp) <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "", names(temp))

# Save out the DE results
write.xlsx(temp,
           "Output/DE/Controls_1/Final/BM_10x_atlas/Controls_1_vs_BM_10x_atlas_FindMarkers_All.xlsx")

# Fix names for items in list for pretty Volcano plots
names(Controls_1_vs_BM_10x_atlas_FindMarkers_All) <- 
  gsub(" ", " ~ ", names(Controls_1_vs_BM_10x_atlas_FindMarkers_All))
names(Controls_1_vs_BM_10x_atlas_FindMarkers_All) <- 
  gsub("[//.]", " ~ ", names(Controls_1_vs_BM_10x_atlas_FindMarkers_All))

# Output Volcano plots for all tests
for (i in 1:length(Controls_1_vs_BM_10x_atlas_FindMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_1_vs_BM_10x_atlas_FindMarkers_All[[`i`]], 
                     lab = Controls_1_vs_BM_10x_atlas_FindMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = parse(text = names(Controls_1_vs_BM_10x_atlas_FindMarkers_All)[i]), 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_1/Final/BM_10x_atlas/",
                  names(Controls_1_vs_BM_10x_atlas_FindMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
  
}

# 5.2.A.iii DE testing - FindMarkers - Controls #1 vs. PB_10x_atlas --------------------------------

# Create a list to store DE results
Controls_1_vs_PB_10x_atlas_FindMarkers_All <- list()

# Iterate over list to test TTE, TEM and IL7R_TM
for (i in (names(Controls_query_list[["Controls_1"]][6:7]))) {
  for (j in cluster_names_PB_10x_atlas[c(1:3)]) {
    temp <- find.discriminant.genes(all_atlas_list$PB_10x_atlas, 
                                    Controls_query_list$Controls_1[[i]], state = j)
    Controls_1_vs_PB_10x_atlas_FindMarkers_All[[i]][[j]] <- temp %>% 
      rownames_to_column() %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
  }
}

# Flatten the list to allow export and fix names
Controls_1_vs_PB_10x_atlas_FindMarkers_All <- 
  unlist(Controls_1_vs_PB_10x_atlas_FindMarkers_All, recursive = FALSE)

# Copy list and fix names to allow export to excel without error
temp <- Controls_1_vs_PB_10x_atlas_FindMarkers_All
names(temp) <- gsub("([.])", " ", names(temp))
names(temp) <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "", names(temp))

# Save out the DE results
write.xlsx(temp,
           "Output/DE/Controls_1/Final/PB_10x_atlas/Controls_1_vs_PB_10x_atlas_FindMarkers_All.xlsx")

# Fix names for items in list for pretty Volcano plots
names(Controls_1_vs_PB_10x_atlas_FindMarkers_All) = 
  gsub(" ", " ~ ", names(Controls_1_vs_PB_10x_atlas_FindMarkers_All))
names(Controls_1_vs_PB_10x_atlas_FindMarkers_All) = 
  gsub("[//.]", " ~ ", names(Controls_1_vs_PB_10x_atlas_FindMarkers_All))

# Output Volcano plots for all tests
for (i in 1:length(Controls_1_vs_PB_10x_atlas_FindMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_1_vs_PB_10x_atlas_FindMarkers_All[[`i`]], 
                     lab = Controls_1_vs_PB_10x_atlas_FindMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = parse(text = names(Controls_1_vs_PB_10x_atlas_FindMarkers_All)[i]), # Is B a regular expression?
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_1/Final/PB_10x_atlas/",
                  names(Controls_1_vs_PB_10x_atlas_FindMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}

# 5.2.B.i DE testing - FindMarkers - Controls #2 vs. BM_10x_atlas ----------------------------------

# Create a list to store DE results
Controls_2_vs_BM_10x_atlas_FindMarkers_All <- list()

# Iterate over list to test TEM and IL7R_TM
for (i in names(Controls_query_list[["Controls_2"]])) {
  for (j in cluster_names_BM_10x_atlas[c(1,3)]) {
    temp <- find.discriminant.genes(all_atlas_list$BM_10x_atlas, 
                                    Controls_query_list$Controls_2[[i]], state = j)
    Controls_2_vs_BM_10x_atlas_FindMarkers_All[[i]][[j]] <- temp %>% 
      rownames_to_column() %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
  }
}

# Flatten the list to allow export and fix names
Controls_2_vs_BM_10x_atlas_FindMarkers_All <- 
  unlist(Controls_2_vs_BM_10x_atlas_FindMarkers_All, recursive = FALSE)

# Copy list and fix names to allow export to excel without error
temp <- Controls_2_vs_BM_10x_atlas_FindMarkers_All
names(temp) <- gsub("([.])", " ", names(temp))
names(temp) <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "", names(temp))

# Save out the DE results
write.xlsx(temp,
           "Output/DE/Controls_2/Final/BM_10x_atlas/Controls_2_vs_BM_10x_atlas_FindMarkers_All.xlsx")

# Fix names for items in list for pretty Volcano plots
names(Controls_2_vs_BM_10x_atlas_FindMarkers_All) <- 
  gsub("Controls 2 complete", "Complete", names(Controls_2_vs_BM_10x_atlas_FindMarkers_All))
names(Controls_2_vs_BM_10x_atlas_FindMarkers_All) = 
  gsub(" ", " ~ ", names(Controls_2_vs_BM_10x_atlas_FindMarkers_All))
names(Controls_2_vs_BM_10x_atlas_FindMarkers_All) <- 
  gsub("[//.]", " ~ ", names(Controls_2_vs_BM_10x_atlas_FindMarkers_All))

# Output Volcano plots for all tests
for (i in 1:length(Controls_2_vs_BM_10x_atlas_FindMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_2_vs_BM_10x_atlas_FindMarkers_All[[`i`]], 
                     lab = Controls_2_vs_BM_10x_atlas_FindMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = parse(text = names(Controls_2_vs_BM_10x_atlas_FindMarkers_All)[i]), 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_2/Final/BM_10x_atlas/",
                  names(Controls_2_vs_BM_10x_atlas_FindMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}


# 5.2.C.i DE testing - FindMarkers - Controls #3 vs. PB_10x_atlas ----------------------------------

# Create a list to store DE results
Controls_3_vs_PB_10x_atlas_FindMarkers_All <- list()

# Iterate over list to test TTE, TEM & TN
for (i in names(Controls_query_list[["Controls_3"]])) {
  for (j in cluster_names_PB_10x_atlas[c(1:2,4)]) {
    temp <- find.discriminant.genes(all_atlas_list$PB_10x_atlas, 
                                    Controls_query_list$Controls_3[[i]], state = j)
    Controls_3_vs_PB_10x_atlas_FindMarkers_All[[i]][[j]] <- temp %>% 
      rownames_to_column() %>% 
      inner_join(y = annotations[, c("gene_name", "description")], 
                 by = c("rowname" = "gene_name")) %>% unique() %>% rename_at(1,~"gene")
  }
}

# Flatten the list to allow export and fix names
Controls_3_vs_PB_10x_atlas_FindMarkers_All <- 
  unlist(Controls_3_vs_PB_10x_atlas_FindMarkers_All, recursive = FALSE)

# Copy list and fix names to allow export to excel without error
temp <- Controls_3_vs_PB_10x_atlas_FindMarkers_All
names(temp) <- gsub("([.])", " ", names(temp))
names(temp) <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "", names(temp))

# Save out the DE results
write.xlsx(Controls_3_vs_PB_10x_atlas_FindMarkers_All,
           "Output/DE/Controls_3/Final/PB_10x_atlas/Controls_3_vs_PB_10x_atlas_FindMarkers_All.xlsx")

# Fix names for items in list for pretty Volcano plots
names(Controls_3_vs_PB_10x_atlas_FindMarkers_All) <- 
  gsub("Controls 3 complete", "Complete", names(Controls_3_vs_PB_10x_atlas_FindMarkers_All))
names(Controls_3_vs_PB_10x_atlas_FindMarkers_All) = 
  gsub(" ", " ~ ", names(Controls_3_vs_PB_10x_atlas_FindMarkers_All))
names(Controls_3_vs_PB_10x_atlas_FindMarkers_All) <- 
  gsub("[//.]", " ~ ", names(Controls_3_vs_PB_10x_atlas_FindMarkers_All))

# Output Volcano plots for all tests
for (i in 1:length(Controls_3_vs_PB_10x_atlas_FindMarkers_All)) {
  temp <- 
    (EnhancedVolcano(Controls_3_vs_PB_10x_atlas_FindMarkers_All[[`i`]], 
                     lab = Controls_3_vs_PB_10x_atlas_FindMarkers_All[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = parse(text = names(Controls_3_vs_PB_10x_atlas_FindMarkers_All)[i]), 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf,
                     xlim = c(-5,5)))
  ggsave(plot = 
           temp, filename = 
           paste0("Output/Figures/Volcanoplots/Controls_3/Final/PB_10x_atlas/",
                  names(Controls_3_vs_PB_10x_atlas_FindMarkers_All)[i], ".png")
         , width = 66, height = 33, units = "cm", bg = "transparent")
}


# 6. Save out data and print stats -----------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean the environment and save out the DE analysis
rm(list = ls()[!ls() %in% c("Controls_1_vs_all_10x_atlas_FindAllMarkers_All",
                            "Controls_1_vs_all_10x_atlas_FindMarkers_All",
                            "Controls_1_vs_BM_10x_atlas_FindAllMarkers_All",
                            "Controls_1_vs_BM_10x_atlas_FindMarkers_All",
                            "Controls_1_vs_PB_10x_atlas_FindAllMarkers_All",
                            "Controls_1_vs_PB_10x_atlas_FindMarkers_All",
                            "Controls_2_vs_BM_10x_atlas_FindAllMarkers_All",
                            "Controls_2_vs_BM_10x_atlas_FindMarkers_All",
                            "Controls_3_vs_PB_10x_atlas_FindAllMarkers_All",
                            "Controls_3_vs_PB_10x_atlas_FindMarkers_All")])

# Save out the DE_results for further analysis
save.image("Data/R_out/ProjecTILs/Projections_DE_Analysis.rds")

gc()

print("#># Finished running '7. DE testing of projected data' script")