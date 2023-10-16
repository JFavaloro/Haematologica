# Script information -------------------------------------------------------------------------------

# Title: scGSEA - Analysis of 10x data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will analyse the results of the GSEA using the vignette:
# http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '2. scGSEA Analysis of 10x data' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, SingleCellExperiment, pheatmap, escape, dittoSeq, openxlsx)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load 10x GSEA results
all_10x <- readRDS("Data/R_out/GSEA/all_10x.rds")
ES_all_10x <- readRDS("Data/R_out/GSEA/ES_all_10x.rds")

# Create character vectors of cluster names
cluster_names_all_10x <- c("TEM", "TTE", "TN", "Cyto_TEM", "PRE_EX", "TCM")

# Create character vectors of colours for plotting clusters # FIX THIS
# TN = lime green ("32CD32"), TCM = corn flower blue ("6495ED"), TEM = crimson ("DC143C"),
# Cyto_TEM = dark orange ("FF8C00"), TTE = gold ("FFD700"), Activated = blue violet ("8A2BE2"),
# IFN = hot pink ("FF69B4"), Spare = medium turquoise ("48D1CC")

ccolours_all_10x <- c("#DC143C", "#FFD700", "#32CD32", "#6495ED", "#FF8C00", "#8A2BE2") #FIX THIS

# Get Hallmark genesets from MSigDB
GS <- getGeneSets(library = "H")

# Write genesets to character vector for later plotting
GS_vector <- names(GS)
GS_vector <- gsub("HALLMARK_", "", GS_vector) %>% str_to_sentence()


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 Subset out metadata to perform statistics ----------------------------------------------------

meta_data <- all_10x@meta.data
meta_data <- meta_data[, -(1:11)]
meta_data <- meta_data[, -(2:8)]
meta_data <- meta_data[, -(3)]

# Reorder clusters as per ontogeny
meta_data <- meta_data %>%
  mutate(Cluster = fct_relevel(Cluster, "T[N]","T[CM]","T[EM]","T[TE]","Cyto ~ T[EM]","P[RE] ~ Ex"))

# 3.2 Subset metadata to allow comparison between Tissue -------------------------------------------

# Create an empty list to store results
meta_list_subsets <- list()

# Subset Clusters
meta_list_subsets$TEM <- subset(meta_data, subset = (Cluster == "T[EM]"))
meta_list_subsets$TTE <- subset(meta_data, subset = (Cluster == "T[TE]"))
meta_list_subsets$TN <- subset(meta_data, subset = (Cluster == "T[N]"))
meta_list_subsets$Cyto_TEM <- subset(meta_data, subset = (Cluster == "Cyto ~ T[EM]"))
meta_list_subsets$PRE_EX <- subset(meta_data, subset = (Cluster == "P[RE] ~ Ex"))
meta_list_subsets$TCM <- subset(meta_data, subset = (Cluster == "T[CM]"))


# 4. Analyse data  ---------------------------------------------------------------------------------

# 4.1 Perform statistical analysis of results ------------------------------------------------------

# Create empty lists to store results of statistics
stats_T_test_10x <- list()

# Run t test comparing tissue restricted differences across clusters analyses
for (i in names(meta_list_subsets)) {
  temp <- getSignificance(meta_list_subsets[[`i`]], 
                          group = "Tissue",
                          fit = "T.test") %>% rownames_to_column()
  stats_T_test_10x[[i]] <- temp
}

# Save out the stats results
write.xlsx(stats_T_test_10x,
           "Output/Stats/10x/Final/stats_T_test_10x.xlsx")

# Run Kruskal Wallis test
stats_KW_10x <- getSignificance(meta_data, 
                                group = "Cluster", fit = "KW") %>% rownames_to_column()

# Save out the stats results
write.xlsx(stats_KW_10x,
           "Output/Stats/10x/Final/stats_KW_10x.xlsx")

# Run ANOVA
stats_ANOVA_10x <- getSignificance(meta_data, 
                                   group = "Cluster", fit = "ANOVA") %>% rownames_to_column()

# Save out the stats results
write.xlsx(stats_ANOVA_10x,
           "Output/Stats/10x/Final/stats_ANOVA_10x.xlsx")


# 4.2 Output figures -------------------------------------------------------------------------------

# 4.2A Output Split enrichment Violin plots --------------------------------------------------------
for (i in 1:length(GS_vector)) {
  temp <- 
    (splitEnrichment(meta_data,
                     x.axis = "Cluster", split = "Tissue", 
                     gene.set = GS_vector[`i`], colors = c("blue", "red"))) +
    scale_x_discrete("Cluster", labels = parse(text = levels(meta_data$Cluster))) + 
    theme(axis.text.x = element_text(size = 16))
  ggsave(plot = temp, 
         filename = paste0("Output/Figures/Violinplots/10x/Final/GSEA/Split/",  
                           GS_vector[`i`], ".png"), 
         units = "cm", height = 15, width = 25, bg = "transparent")
}


# 4.2B Output side-by side Violin plots ------------------------------------------------------------
for (i in 1:length(GS_vector)) {
  temp <- 
    (dittoPlot(all_10x, var = GS_vector[i], 
               group.by = "Cluster", split.by = "Tissue", 
               plots = c("jitter", "vlnplot", "boxplot"), 
               color.panel = ccolours_all_10x, ylab = "Enrichment Scores" + 
                 scale_x_discrete("Cluster", 
                                  labels = parse(text = levels(all_10x@meta.data$Cluster)))))
  ggsave(plot = temp, 
         filename = paste0("Output/Figures/Violinplots/10x/Final/GSEA/Side_by_side/", 
                           GS_vector[`i`], ".png"))
}


# 4.2C Output Heatmap of GSEA results --------------------------------------------------------------
temp <- dittoHeatmap(all_10x, genes = NULL, metas = names(ES_all_10x), 
             annot.by = c("Cluster", "Tissue"),
             order.by = c("Tissue", "Cluster"),
             fontsize = 7, 
             cluster_cols = FALSE,
             annot.colors = c(ccolours_all_10x, "blue", "red"),
             heatmap.colors = colorRampPalette(c("purple", "black", "yellow"))(50))
temp[["gtable"]][["grobs"]][[6]][["children"]][["GRID.text.24644"]][["label"]] = parse(text = levels(all_10x@meta.data$Cluster)) # This appears to get a different name each time the figure is created. I don't know how to automate this...
ggsave(temp, filename = paste0("Output/Figures/Heatmaps/10x/Final/all_10x_GSEA_heatmap.png"))


# 5. Save out data and print stats -----------------------------------------------------------------

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '2. scGSEA Analysis of 10x data' script")