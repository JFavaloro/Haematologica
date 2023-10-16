# Script information -------------------------------------------------------------------------------

# Title: 10x ProjecTILs/SingleR annotations analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will analysis the output of the ProjecTILs projection using the
# vignette: https://carmonalab.github.io/ProjecTILs.demo/tutorial.html. We will additionally assess
# The SingleR annotations to determine how they compare to our preliminary annotations based on
# biological knowledge of the expression of canonical genes using the vignette: 
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '2. 10x annotations analysis' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, Matrix, ProjecTILs, SingleR, celldex, ensembldb, openxlsx, EnhancedVolcano, Scillus, bluster, EnsDb.Hsapiens.v86)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load data
all_list <- readRDS("Data/R_out/10x/all_list.rds")
ProjecTILs_projections_10x <- readRDS("Data/R_out/ProjecTILs/ProjecTILs_projections_10x.rds")
SingleR_predictions_10x <- readRDS("Data/R_out/ProjecTILs/SingleR_predictions_10x.rds")

# Load additional resources
black_list <- readRDS("Data/R_out/black_list.rds")

# Modify black_list to retain TCR genes
black_list_TCR <- 
  black_list[!grepl("TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC", black_list)]

# Load ProjecTILs reference atlases
ProjecTILs_ref <- 
  list("ref_TIL" = load.reference.map(ref = 
                                        "Data/Public_data/ProjecTILs/ref_TIL_Atlas_mouse_v1.rds"), 
       "ref_CMV" = load.reference.map(ref = 
                                        "Data/Public_data/ProjecTILs/ref_LCMV_Atlas_mouse_v1.rds"))

# Load reference datasets for SingleR
SingleR_ref <- readRDS("Data/R_out/ProjecTILs/SingleR_ref.rds")

# Create character vectors of title plots
plot_titles <- c("All data projected on tumour infiltrating lymphocyte dataset",
                 "All data projected on chronic infection dataset",
                 "BM data projected on tumour infiltrating lymphocyte dataset",
                 "BM data projected on chronic infection dataset",
                 "PB data projected on tumour infiltrating lymphocyte dataset",
                 "PB data projected on chronic infection dataset")

# Create character vectors of sample_names
sample_names <- c("All data - TIL", "All data - CMV", 
                  "BM data - TIL", "BM data - CMV", 
                  "PB data - TIL", "PB data - CMV")

# Create character vectors of cluster_names
cluster_names_TILs <- c("CD8_Tex", "CD8_Tpex", "CD8_EffectorMemory",
                        "CD8_EarlyActiv", "CD8_NaiveLike", "CD4_NaiveLike",
                        "Tfh", "Th1", "Treg")
cluster_names_CMV <- c("CD8_Tex", "CD8_Tpex", "SLEC", "Eff Interm",
                       "Eff Cycling", "Eff Early", "Memory Prec")

# Create colours
ccolours_TILs <- c("#EDBE2A", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", 
                   "#D1CfCC", "#FF0000", "#87F6A5", "#E812DD")
ccolours_CMV <- c("#EDBE2A", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#D1CfCC", "#FF0000")


# 3. Pre-processing --------------------------------------------------------------------------------

# Duplicate the ProjecTILs refs to allow plotting
ProjecTILs_ref <- rep(ProjecTILs_ref, 3)

# Switch default assay to "RNA" for differential expression testing
all_list_RNA <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Remove non helpful genes from DE results #NB: This will remove genes from the object
for (i in names (all_list_RNA)) {
  counts <- GetAssayData(all_list_RNA[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  all_list_RNA[[i]] <- subset(all_list_RNA[[i]], features = rownames(counts))
}

# Add additional metadata to 10x projections to allow easier plotting of data
for (i in 1:length(ProjecTILs_projections_10x)){
  ProjecTILs_projections_10x[[i]][["Sample"]] <- sample_names[i]
}

# Set the Identity of the ProjecTILs_projections_10x list to "functional.cluster"
ProjecTILs_projections_10x <- 
  map(ProjecTILs_projections_10x, `Idents<-`, value = "functional.cluster")


# 4. Prepare figures and DE output -----------------------------------------------------------------

# 4.A Plot projections
for (i in 1:length(ProjecTILs_projections_10x)) {
  temp = 
    (plot.projection(ref = ProjecTILs_ref[[`i`]], 
                     query = ProjecTILs_projections_10x[[`i`]]) + ggtitle(plot_titles[i]) + 
       theme(plot.title = element_text(size = 10)))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/ProjecTILs/10x/Projections/Dimplots/Projection_",
                                   names(ProjecTILs_projections_10x)[i], ".png"))
}


# 4.B Plot barplots

# Plot the TIL overlay
# Merge Seurat objects to allow them to be plotted together
temp <- merge(ProjecTILs_projections_10x$all_10x_TIL, 
              c(ProjecTILs_projections_10x$BM_10x_TIL, ProjecTILs_projections_10x$PB_10x_TIL))

# Extract metrics to plot
temp <- table(temp$functional.cluster, temp$Sample) %>% as.data.frame()
temp <- temp %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Barplots of cluster distribution of all TIL projections - ask Sam why this won't plot - levels should be set
ggplot(temp, 
       aes(fill = factor(Var1, levels = 8:0), y = Proportion, x = Var2)) + 
  geom_bar(position = "stack", stat = 'identity') + 
  scale_fill_manual(values = ccolours_TILs, name = "Cluster", 
                    labels = cluster_names_TILs) + 
  labs(x = "Sample", y = "Proportion (%)", title = "Data projected on TIL atlas") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/General/10x/Final/Barplot/Barplot_",
                "Data clustered by TIL atlas", ".png"))

# Plot the CMV overlay
# Merge Seurat objects to allow them to be plotted together
temp <- merge(ProjecTILs_projections_10x$all_10x_CMV, 
              c(ProjecTILs_projections_10x$BM_10x_CMV, ProjecTILs_projections_10x$PB_10x_CMV))

# Extract metrics to plot
temp <- table(temp$functional.cluster, temp$Sample) %>% as.data.frame()
temp <- temp %>% group_by(Var2) %>% mutate(Proportion = Freq/sum(Freq) *100)

# Barplots of cluster distribution of all CMV projections - ask Sam why this won't plot
ggplot(temp, 
       aes(fill = factor(Var1, levels = 6:0), y = Proportion, x = Var2)) + 
  geom_col(position = "stack") + 
  scale_fill_manual(values = ccolours_CMV, name = "Cluster", 
                    labels = cluster_names_CMV) + 
  labs(x = "Sample", y = "Proportion (%)", title = "Data projected on CMV atlas") + theme_bw()
ggsave(filename = 
         paste0("Output/Figures/General/10x/Final/Barplot/Barplot_",
                "Data clustered by CMV atlas", ".png"))


# 4.C Plot projections - Working OK
for (i in 1:length(SingleR_predictions_10x)) {
  temp = 
    (plotScoreHeatmap(SingleR_predictions_10x[[i]]))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/SingleR/10x/Annotations/Heatmaps/Annotation_",
                                   names(SingleR_predictions_10x)[i], ".png"))
}


# 5. Save out data and print stats -----------------------------------------------------------------

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '2. 10x Annotation analysis' script")