# Script information -------------------------------------------------------------------------------

# Title: scGSEA - Analysis of Controls data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will analyse the results of the GSEA using the vignette:
# http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '3. scGSEA Analysis of Controls data' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, SingleCellExperiment, pheatmap, escape, dittoSeq, openxlsx, multtest, metap, EnhancedVolcano, ProjecTILs, ggpubr)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load Controls GSEA results
merged_list <- readRDS("Data/R_out/GSEA/merged_list.rds")
ES_merged_list <- readRDS("Data/R_out/GSEA/ES_merged_list.rds")

# Create character vectors of cluster names
cluster_names_BM_10x_atlas <- c("TEM", "TN", "IL7R_TM", "TTE", "Activated", "Cyto_TEM")
cluster_names_BM_10x_atlas <- c("TEM", "TN", "IL7R_TM", "TTE")
cluster_names_PB_10x_atlas <- c("TTE", "TEM", "IL7R_TM", "TN", "Activated", "IFN")
cluster_names_PB_10x_atlas <- c("TTE", "TEM", "IL7R_TM", "TN")
trimmed_names <- c("TN", "IL7R_TM", "TEM", "TTE")
null_list = c("Activated", "Cyto_TEM", "IFN")

# Create character vectors of colours for plotting clusters # FIX THIS
# TN = lime green ("32CD32"), IL7R_TM = corn flower blue ("6495ED"), TEM = crimson ("DC143C"),
# Cyto_TEM = dark orange ("FF8C00"), TTE = gold ("FFD700"), Activated = blue violet ("8A2BE2"),
# IFN = hot pink ("FF69B4"), Spare = medium turquoise ("48D1CC")

ccolours_BM_10x_atlas <- c("#DC143C", "#32CD32", "#6495ED", "#FFD700", "#8A2BE2", "#FF8C00")
ccolours_BM_10x_atlas <- c("#DC143C", "#32CD32", "#6495ED", "#FFD700")
ccolours_PB_10x_atlas <- c("#FFD700", "#DC143C", "#6495ED", "#32CD32", "#8A2BE2", "#FF69B4")
ccolours_PB_10x_atlas <- c("#FFD700", "#DC143C", "#6495ED", "#32CD32")

# Get Hallmark genesets from MSigDB
GS <- getGeneSets(library = "H")

# Write genesets to character vector for later plotting
GS_vector <- names(GS)
GS_vector <- gsub("HALLMARK_", "", GS_vector) %>% str_to_sentence()


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 Remove data for clusters with insufficient data to compare and fix object for plotting -------
for (i in 1:length(merged_list)) {
  temp <- subset(merged_list[[i]], 
                 subset = (functional.cluster == "IL7R^`+` ~ T[M]" | functional.cluster == "T[EM]" | functional.cluster == "T[TE]"))
  merged_list[[i]] <- temp
}

# Reorder Idents for plotting
for (i in names(merged_list)) {
  merged_list[[`i`]]@meta.data$functional.cluster <- 
    factor(merged_list[[`i`]]@meta.data$functional.cluster, 
           levels = c("IL7R^`+` ~ T[M]", "T[EM]", "T[TE]"))
}


# 3.2 Subset out metadata to perform statistics ----------------------------------------------------

# Create an empty list to store metadata
meta_list <- list()

for (i in names(merged_list)[1:4]) {
  temp <- merged_list[[`i`]]@meta.data
  colnames(temp)[ncol(temp)] <- "functional.cluster"
  colnames(temp) <- make.unique(names(temp))
  temp <- temp[, -(1:26)]
  temp <- temp[, -(3:13)]
  temp <- temp[, -(53:54)]
  meta_list[[i]] <- temp
}  

for (i in names(merged_list)[5:6]) {
  temp <- merged_list[[`i`]]@meta.data 
  colnames(temp)[ncol(temp)] <- "functional.cluster"
  colnames(temp) <- make.unique(names(temp))
  temp <- temp[, -(1:26)]
  temp <- temp[, -(3:7)]
  temp <- temp[, -(53:54)]
  meta_list[[i]] <- temp
}

# Convert Sample to factor
for (i in 1:length(meta_list)) {
  meta_list[[i]]$Sample <- as.factor(meta_list[[i]]$Sample)
}

# Reset levels of Sample in the meta_list object to graph the reference first
meta_list$`Activated BM`$Sample <- 
  factor(meta_list$`Activated BM`$Sample, 
         levels = as.factor(c("MM Reference (BM)", "Activated BM")))
meta_list$`Resting BM`$Sample <- 
  factor(meta_list$`Resting BM`$Sample, 
         levels = as.factor(c( "MM Reference (BM)", "Resting BM")))
meta_list$`Activated PB`$Sample <- 
  factor(meta_list$`Activated PB`$Sample, 
         levels = as.factor(c("MM Reference (PB)", "Activated PB")))
meta_list$`Resting PB`$Sample <- 
  factor(meta_list$`Resting PB`$Sample, 
         levels = as.factor(c("MM Reference (PB)", "Resting PB")))
meta_list$`Age-matched PB`$Sample <- 
  factor(meta_list$`Age-matched PB`$Sample, 
         levels = as.factor(c("MM Reference (PB)", "Age-matched Controls PB")))
meta_list$`Young PB`$Sample <- 
  factor(meta_list$`Young PB`$Sample, 
         levels = as.factor(c( "MM Reference (PB)", "Young PB")))


# 3.3 Subset metadata to allow comparison between Sample -------------------------------------------

# Subset TN cells
meta_TN <- meta_list

for (i in 1:length(meta_TN)) {
  temp <- subset(meta_TN[[i]], subset = (functional.cluster == "T[N]"))
  meta_TN[[i]] <- temp
}

# Subset IL7R_TM cells
meta_IL7R_TM <- meta_list

for (i in 1:length(meta_IL7R_TM)) {
  temp <- subset(meta_IL7R_TM[[i]], subset = (functional.cluster == "IL7R^`+` ~ T[M]"))
  meta_IL7R_TM[[i]] <- temp
}

# Subset TEM cells
meta_TEM <- meta_list

for (i in 1:length(meta_TEM)) {
  temp <- subset(meta_TEM[[i]], subset = (functional.cluster == "T[EM]"))
  meta_TEM[[i]] <- temp
}

# Subset TTE cells
meta_TTE <- meta_list

for (i in 1:length(meta_TTE)) {
  temp <- subset(meta_TTE[[i]], subset = (functional.cluster == "T[TE]"))
  meta_TTE[[i]] <- temp
}

# Move results into a list
meta_list_subsets <- 
  list("TN" = meta_TN, "IL7R_TM" = meta_IL7R_TM, "TEM" = meta_TEM, "TTE" = meta_TTE)

# Remove objects with insufficient data to test
meta_list_subsets$TN$`Activated BM` <- NULL
meta_list_subsets$TN$`Resting BM` <- NULL
meta_list_subsets$TN$`Resting PB` <- NULL
meta_list_subsets$IL7R_TM$`Age-matched PB` <- NULL

# Flatten the list
meta_list_subsets <- meta_list_subsets %>% unlist(recursive = FALSE)


# 4. Analyse data  ---------------------------------------------------------------------------------

# 4.1 Perform statistical analysis of results ------------------------------------------------------

# Create empty lists to store results of statistics
stats_T_test_10x_atlas <- list()
stats_KW_10x_atlas <- list()
stats_ANOVA_10x_atlas <- list()

# Run statistical analyses
for (i in names(meta_list_subsets)) {
  temp <- getSignificance(meta_list_subsets[[`i`]], 
                          group = "Sample",
                          fit = "T.test") %>% rownames_to_column()
  stats_T_test_10x_atlas[[i]] <- temp
}

# Save out the stats results
write.xlsx(stats_T_test_10x_atlas,
           "Output/Stats/10x_atlas/Final/stats_T_test_10x_atlas.xlsx")

for (i in names(meta_list)) {
  temp <- getSignificance(meta_list[[`i`]], 
                          group = "functional.cluster",
                          fit = "KW") %>% rownames_to_column()
  stats_KW_10x_atlas[[i]] <- temp
}

# Save out the stats results
write.xlsx(stats_KW_10x_atlas,
           "Output/Stats/10x_atlas/Final/stats_KW_10x_atlas.xlsx")

for (i in names(meta_list)) {
  temp <- getSignificance(meta_list[[`i`]], 
                          group = "functional.cluster",
                          fit = "ANOVA") %>% rownames_to_column()
  stats_ANOVA_10x_atlas[[i]] <- temp
}

# Save out the stats results
write.xlsx(stats_ANOVA_10x_atlas,
           "Output/Stats/10x_atlas/Final/stats_ANOVA_10x_atlas.xlsx")


# 4.2 Output figures -------------------------------------------------------------------------------

# 4.2Ai Output Split enrichment Violin plots for Controls_1 ----------------------------------------
for (i in names(meta_list)[1]) {
  for (j in 1:length(GS_vector)) {
    temp <- 
      (splitEnrichment(meta_list[[`i`]],
                       x.axis = "functional.cluster", split = "Sample", gene.set = GS_vector[`j`], 
                       colors = c("blue", "purple"))) +
      scale_x_discrete("Cluster", 
                       labels = parse(text = levels(meta_list[[`i`]]$functional.cluster))) + 
      theme(axis.text.x = element_text(size = 16)) + 
      scale_y_continuous(limits = c(1500, 3500), breaks = seq(1500, 3500, by = 1000)) + 
      theme(axis.text.y = element_text(size = 16)) + theme(legend.position="none")
    ggsave(plot = temp, 
           filename = paste0("Output/Figures/Violinplots/Controls_1/Final/GSEA/Split/", 
                             names(merged_list[`i`]), "/", GS_vector[`j`], ".png"), 
           units = "cm", height = 15, width = 10, bg = "transparent")
  }
}

for (i in names(meta_list)[2]) {
  for (j in 1:length(GS_vector)) {
    temp <- 
      (splitEnrichment(meta_list[[`i`]],
                       x.axis = "functional.cluster", split = "Sample", gene.set = GS_vector[`j`], 
                       colors = c("blue", "green"))) +
      scale_x_discrete("Cluster", 
                       labels = parse(text = levels(meta_list[[`i`]]$functional.cluster))) + 
      theme(axis.text.x = element_text(size = 16)) + 
      scale_y_continuous(limits = c(1500, 3500), breaks = seq(1500, 3500, by = 1000)) + 
      theme(axis.text.y = element_text(size = 16)) + theme(legend.position="none")
    ggsave(plot = temp, 
           filename = paste0("Output/Figures/Violinplots/Controls_1/Final/GSEA/Split/", 
                             names(merged_list[`i`]), "/", GS_vector[`j`], ".png"), 
           units = "cm", height = 15, width = 10, bg = "transparent")
  }
}

for (i in names(meta_list)[3]) {
  for (j in 1:length(GS_vector)) {
    temp <- 
      (splitEnrichment(meta_list[[`i`]],
                       x.axis = "functional.cluster", split = "Sample", gene.set = GS_vector[`j`], 
                       colors = c("red", "purple"))) +
      scale_x_discrete("Cluster", 
                       labels = parse(text = levels(meta_list[[`i`]]$functional.cluster))) + 
      theme(axis.text.x = element_text(size = 16)) + 
      scale_y_continuous(limits = c(1500, 3500), breaks = seq(1500, 3500, by = 1000)) + 
      theme(axis.text.y = element_text(size = 16)) + theme(legend.position="none")
    ggsave(plot = temp, 
           filename = paste0("Output/Figures/Violinplots/Controls_1/Final/GSEA/Split/", 
                             names(merged_list[`i`]), "/", GS_vector[`j`], ".png"), 
           units = "cm", height = 15, width = 10, bg = "transparent")
  }
}

for (i in names(meta_list)[4]) {
  for (j in 1:length(GS_vector)) {
    temp <- 
      (splitEnrichment(meta_list[[`i`]],
                       x.axis = "functional.cluster", split = "Sample", gene.set = GS_vector[`j`], 
                       colors = c("red", "green"))) +
      scale_x_discrete("Cluster", 
                       labels = parse(text = levels(meta_list[[`i`]]$functional.cluster))) + 
      theme(axis.text.x = element_text(size = 16)) + 
      scale_y_continuous(limits = c(1500, 3500), breaks = seq(1500, 3500, by = 1000)) + 
      theme(axis.text.y = element_text(size = 16)) + theme(legend.position="none")
    ggsave(plot = temp, 
           filename = paste0("Output/Figures/Violinplots/Controls_1/Final/GSEA/Split/", 
                             names(merged_list[`i`]), "/", GS_vector[`j`], ".png"), 
           units = "cm", height = 15, width = 10, bg = "transparent")
  }
}

# 4.2Aii Output Split enrichment Violin plots for Controls_3 ---------------------------------------
for (i in names(meta_list)[5:6]) {
  for (j in 1:length(GS_vector)) {
    temp <- 
      (splitEnrichment(meta_list[[`i`]],
                       x.axis = "functional.cluster", split = "Sample", gene.set = GS_vector[`j`])) +
      scale_x_discrete("functional.cluster", 
                       labels = parse(text = levels(meta_list[[`i`]]$functional.cluster)))
    ggsave(plot = temp, 
           filename = paste0("Output/Figures/Violinplots/Controls_3/Final/GSEA/Split/", 
                             names(merged_list[`i`]), "/", GS_vector[`j`], ".png"))
  }
}

# 4.2Bi Output Violin plots of all BM comparisons, arranged by cluster, split by Sample ------------
for (i in names(merged_list)[1:2]) {
  for (j in 1:length(GS_vector)) {
    temp <- 
      (dittoPlot(merged_list[[`i`]], var = GS_vector[j], 
                 group.by = "functional.cluster", split.by = "Sample", 
                 plots = c("jitter", "vlnplot", "boxplot"), ylab = "Enrichment Scores", 
                 color.panel = ccolours_BM_10x_atlas)) + 
      scale_x_discrete("functional.cluster", 
                       labels = parse(text = levels(merged_list[[`i`]]@meta.data$functional.cluster)))
    ggsave(plot = temp, 
           filename = paste0("temp/", 
                             names(merged_list[`i`]), "/", GS_vector[`j`], ".png"))
  }
}

# Move folders to appropriate locations NB: This requires creation of folders PRIOR to moving
file.copy("temp/Activated BM/", "Output/Figures/Violinplots/Controls_1/Final/GSEA/Side_by_side/Activated BM", overwrite = TRUE)
file.copy("temp/Resting BM/", "Output/Figures/Violinplots/Controls_1/Final/GSEA/Side_by_side/Resting BM", overwrite = TRUE)


# 4.2Bii Output Violin plots of all PB comparisons, arranged by cluster, split by Sample -----------
for (i in names(merged_list)[3:6]) {
  for (j in 1:length(GS_vector)) {
    temp <- 
      (dittoPlot(merged_list[[`i`]], var = GS_vector[j], 
                 group.by = "functional.cluster", split.by = "Sample", 
                 plots = c("jitter", "vlnplot", "boxplot"), ylab = "Enrichment Scores", 
                 color.panel = ccolours_PB_10x_atlas)) + 
      scale_x_discrete("functional.cluster", 
                       labels = parse(text = levels(merged_list[[`i`]]@meta.data$functional.cluster)))
    ggsave(plot = temp, 
           filename = paste0("temp/", 
                             names(merged_list[`i`]), "/", GS_vector[`j`], ".png"))
  }
}

# Move folders to appropriate locations NB: This requires creation of folders PRIOR to moving
file.copy("temp/Activated PB/", 
          "Output/Figures/Violinplots/Controls_1/Final/GSEA/Side_by_side/Activated PB", overwrite = TRUE)
file.copy("temp/Resting PB/", 
          "Output/Figures/Violinplots/Controls_1/Final/GSEA/Side_by_side/Resting PB", overwrite = TRUE)
file.copy("temp/Age-matched PB/", 
          "Output/Figures/Violinplots/Controls_3/Final/GSEA/Side_by_side/Age-matched PB", overwrite = TRUE)
file.copy("temp/Young PB/", 
          "Output/Figures/Violinplots/Controls_3/Final/GSEA/Side_by_side/Young PB", overwrite = TRUE)


# 4.2Ci Output heatmaps of all BM comparisons, arranged by cluster, split by Sample ----------------
for (i in names(merged_list)[1:2]) {
  temp <- 
    (dittoHeatmap(merged_list[[`i`]], genes = NULL, metas = names(ES_merged_list[[`i`]]), 
                  annot.by = c("functional.cluster", "Sample"),
                  order.by = c("Sample", "functional.cluster"),
                  fontsize = 7, 
                  cluster_cols = FALSE,
                  annot.colors = c(ccolours_BM_10x_atlas, "green", "purple"),
                  heatmap.colors = colorRampPalette(c("purple", "black", "yellow"))(50)))
  temp[["gtable"]][["grobs"]][[6]][["children"]][["GRID.text.412585"]][["label"]] = parse(text = levels(merged_list[[`i`]]@meta.data$functional.cluster))
  ggsave(plot = temp, 
         filename = paste0("temp/", 
                           names(merged_list[`i`]), "/", "Heatmap", ".png"))
}

# Move folders to appropriate locations NB: This requires creation of folders PRIOR to moving
file.copy("temp/Activated BM/Heatmap.png", 
          "Output/Figures/Heatmaps/Controls_1/Final/GSEA/Activated_BM_Heatmap.png", overwrite = TRUE)
file.copy("temp/Resting BM/Heatmap.png", 
          "Output/Figures/Heatmaps/Controls_1/Final/GSEA/Resting_BM_Heatmap.png", overwrite = TRUE)


# 4.2Cii Output heatmaps of all PB comparisons, arranged by cluster, split by Sample ---------------
for (i in names(merged_list)[3:6]) {
  temp <- 
    (dittoHeatmap(merged_list[[`i`]], genes = NULL, metas = names(ES_merged_list[[`i`]]), 
                  annot.by = c("functional.cluster", "Sample"),
                  order.by = c("Sample", "functional.cluster"),
                  fontsize = 7, 
                  cluster_cols = FALSE,
                  annot.colors = c(ccolours_PB_10x_atlas, "green", "purple"),
                  heatmap.colors = colorRampPalette(c("purple", "black", "yellow"))(50)))
  temp[["gtable"]][["grobs"]][[6]][["children"]][["GRID.text.412585"]][["label"]] = parse(text = levels(merged_list[[`i`]]@meta.data$functional.cluster))
  ggsave(plot = temp, 
         filename = paste0("temp/", 
                           names(merged_list[`i`]), "/", "Heatmap", ".png"))
}

# Move folders to appropriate locations NB: This requires creation of folders PRIOR to moving
file.copy("temp/Activated PB/Heatmap.png", 
          "Output/Figures/Heatmaps/Controls_1/Final/GSEA/Activated_PB_Heatmap.png", overwrite = TRUE)
file.copy("temp/Resting PB/Heatmap.png", 
          "Output/Figures/Heatmaps/Controls_1/Final/GSEA/Resting_PB_Heatmap.png", overwrite = TRUE)
file.copy("temp/Age-matched PB/Heatmap.png", 
          "Output/Figures/Heatmaps/Controls_3/Final/GSEA/Age-matched_PB_Heatmap.png", overwrite = TRUE)
file.copy("temp/Young PB/Heatmap.png", 
          "Output/Figures/Heatmaps/Controls_3/Final/GSEA/Young_PB_Heatmap.png", overwrite = TRUE)


# 5. Save out data and print stats -----------------------------------------------------------------

# Clean the temp directory
unlink("temp", recursive = TRUE)

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '3. scGSEA Analysis of Controls data' script")