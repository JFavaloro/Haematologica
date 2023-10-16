# Script information -----------------------------------------------------------

# title: scGSEA
# author: James Favaloro
# date: 2021-09-02
# description: In this script we will use escape to determine enrichment
# across a number of cellular processess using the vignette:
# http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html


# 1. Import libraries ----------------------------------------------------------

print("#># Start running 'scGSEA' script")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(scRepertoire)
  library(SingleCellExperiment)
  library(pheatmap)
  library(escape)
  library(dittoSeq)
  library(ggpubr)
})

# 2. Load data and define functions---------------------------------------------

integrated.list = readRDS("Data/integrated.list.rds")
black_list = readRDS("Data/black_list.rds")
annotations = read.csv("Data/annotations.csv")

# Define %!in% function
'%!in%' <- function(x,y)!('%in%'(x,y))



# 3. Fix Seurat object ---------------------------------------------------------

# Pull out combined object for analysis and correct clusters and assay
x = integrated.list$all_ds.integrated
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"
DefaultAssay(x) = "RNA"

# Check to ensure settings are correct
DimPlot(x)
VlnPlot(x, "CD69", split.by = "tissue")

# Backup x and remove the integrated list to free up memory
x_backup = x
rm(integrated.list)

# 4. Run Escape

# Get Gene Sets
GS <- getGeneSets(library = "H")

# Write Gene Sets to character vector for latter plotting
GS_vector = names(GS)

# Run scGSEA across all cells in our combined seurat object NB: This will take a while!
ES <- enrichIt(x, gene.sets = GS, groups = 1000, cores = 8)

# Save our the enrichment data
saveRDS(ES, "Data/ES.R")

# Append the results back to the Seurat object as metadata
x <- AddMetaData(x, ES)

# Plot results as heatmap
dittoHeatmap(x, genes = NULL, metas = names(ES), 
             annot.by = c("seurat_clusters", "tissue"),
             order.by = c("tissue", "seurat_clusters"),
             fontsize = 7, 
             cluster_cols = FALSE,
             annot.colors = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#619CFF", "#FF61C3", "red", "blue"),
             heatmap.colors = colorRampPalette(c("purple", "black", "yellow"))(50)) %>% ggsave() # Ask Sam how to fix this to save correctly - low priority.


#Convert Seurat object to SingleCellExperiment and plot gene sets as Violin plots
x.sce = as.SingleCellExperiment(x)

# Copy meta data to SingleCellExperiment
### This is really broken! Run lines 1-4 then run 1,2,4 SKIP 3! ###
met.data <- merge(colData(x.sce), ES, by = "row.names", all=TRUE)
row.names(met.data) <- met.data$Row.names
met.data$Row.names <- NULL # This is breaking things - not sure why this is needed in the first place (?)
colData(x.sce) <- met.data

# Define function to plot all 50 Genesets as Violin plots, arranged by cluster, split by tissue
GS_plot = function(x) {
  multi_dittoPlot(x.sce, vars = x, 
                  group.by = "seurat_clusters", split.by = "tissue", plots = c("jitter", "vlnplot", "boxplot"), color.panel = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#619CFF", "#FF61C3"),
                  ylab = "Enrichment Scores",
                  theme = theme_classic() + theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)))
  paste("Output/Figures/GSEA/GS", x, ".png") %>% ggsave(width = 41.3, height = 23.7, units = "cm")
}

# Run function to plot data
sapply(GS_vector, GS_plot)

# Ask Sam For some reason Androgen response is bugged - seems to be replaced arbitrarily by last figure(?)
# Ask Sam if he can help integrate stat_compare_means to add p values - not easily. Give this a shot myself

# Define function to plot all 50 Genesets as Violin plots, arranged by cluster, split by tissue, including stats
GS_plot_stats = function(x) {
  multi_dittoPlot(x.sce, vars = x, 
                  group.by = "seurat_clusters", split.by = "tissue", plots = c("jitter", "vlnplot", "boxplot"), color.panel = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#619CFF", "#FF61C3"),
                  ylab = "Enrichment Scores",
                  theme = theme_classic() + theme(plot.title = element_text(size = 20)))
  paste("Output/Figures/GSEA/GS", x, ".png") %>% ggsave(width = 41.3, height = 23.7, units = "cm")
}

# Define groups to test
my_compar = list( c("0", "0"), c("0", "1"), c("0", "2"), c("0", "3"), c("0", "4"), c("0", "5"), c("0", "6"), c("1", "1"), c("1", "2"), c("1", "3"), c("1", "4"), c("1", "5"), c("1", "6"), c("2", "2"), c("2", "3"), c("2", "4"), c("2", "5"), c("2", "6"), c("3", "3"), c("3", "4"), c("3", "5"), c("3", "6"), c("4", "4"), c("4", "5"), c("4", "6"), c("5", "5"), c("5", "6"), c("6", "6"))

# stats = list(c("bm", "pb"))
# stats = list(c("0", "1"))

# Not sure how to define y, as it varies across the 50 tests - need a new function?
# compare_means(? ~ c("0", "1", "2", "3", "4", "5", "6"),  data = ES, ref.group = ".all.",      method = "t.test")

# compare_means(Enrichment score ~ dose,  data = ES, method = "t.test")
# Add this to the below function
#+   stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")

# Create ES2 object to work from for split VlnPlot and PCA
ES2 <- data.frame(x[[]], Idents(x))
colnames(ES2)[ncol(ES2)] <- "cluster"

# Write function to output all 50 Gene Sets as split Violin Plots
GS_splitplot = function(x) {
  splitEnrichment(ES2, x.axis = "cluster", split = "tissue", colors = c("red", "blue"), gene.set = x) 
  paste("Output/Figures/GSEA/GS_split", x, ".png") %>% ggsave()
}

# Run function to plot data
sapply(GS_vector, GS_splitplot)

# Ask Sam for help with this - clusters cols as per seurat, tissue red and blue
GS_splitplot = function(x) {
  splitEnrichment(ES2, x.axis = "cluster", split = "tissue", gene.set = x) + 
    theme(scale_fill_manual(values = c("0" = "#F8766D", "1" = "#D39200", "2" = "#93AA00", "3" = "#00BA38", "4" = "#00C19F", "5" = "#619CFF", "6" = "#FF61C3"))) +
    theme(scale_color_manual(values = c("bm" = "red", "pb" = "blue")))
  paste("Output/Figures/GSEA/GS_split", x, ".png") %>% ggsave()
}


# Run PCA to see if we can get a more formal measure of inter-tissue differences
# Cull ES2 to remove NA (missing paired clonotype) and other metadata as PCA cannot run with NA
ES2 = ES2[-c(1:11, 14:25)]

#Run PCA - need to fix this up to save out; piping to ggsave doesn't appear to work(?) but +ggsave does (even if it gives an error) (? ASK SAM)
#  %>% ggsave("Output/Figures/GSEA/PCA.png")
PCA <- performPCA(enriched = ES2, groups = c("tissue", "cluster"))
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE) + ggsave("Output/Figures/GSEA/PCA.png")
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "cluster") + ggsave("Output/Figures/GSEA/PCA_cluster.png")
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = FALSE, facet = "tissue") + ggsave("Output/Figures/GSEA/PCA_tissue.png")
masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)


# Try narrow down our data by subsetting our clutsers and looking at inter-patient differences

# Subset out ES2 by cluster into a list
ES_cluster <- split(ES2, as.factor(ES2$cluster))

# Write function to perfrom PCA over the ES_Cluster list
PCA_cluster = function(x) {
  performPCA(x, groups = c("tissue", "patient"))
}

# Run function
PCA = lapply(ES_cluster, PCA_cluster)


# Write function to perform PCA plotting tissue/cluster #   paste("Output/Figures/GSEA/PCA_plot_tissue", x, ".png") %>% ggsave() ASK SAM - this errors with device null
PCA_plot_tissue = function(x) {
  pcaEnrichment(x, PCx = "PC1", PCy = "PC2", contours = TRUE, facet = "tissue") 
}

# Run function
tissue_plots = lapply(PCA, PCA_plot_tissue)

# Save out plots
for(i in 1:7){
  ggsave(tissue_plots[[i]], file = paste("Output/Figures/GSEA/PCA_plot_tissue",i-1,".png",sep=""))
}


# Write function to perform PCA plotting patient/cluster #   paste("Output/Figures/GSEA/PCA_plot_tissue", x, ".png") %>% ggsave() ASK SAM - this errors with device null
PCA_plot_patient = function(x) {
  pcaEnrichment(x, PCx = "PC1", PCy = "PC2", contours = TRUE, facet = "patient") 
}

# Run function
patient_plots = lapply(PCA, PCA_plot_patient)

# Save out plots
for(i in 1:7){
  ggsave(patient_plots[[i]], file = paste("Output/Figures/GSEA/PCA_plot_patient",i-1,".png",sep=""))
}

# Write function to perform PCA plotting determining which gene sets account for most variance #   paste("Output/Figures/GSEA/PCA_plot_tissue", x, ".png") %>% ggsave() ASK SAM - this errors with device null
PCA_plot_master = function(x) {
  masterPCAPlot(x, PCx = "PC1", PCy = "PC2", top.contribution = 10) 
}

# Run function
master_plots = lapply(ES_cluster, PCA_plot_master)

# Save out plots
for(i in 1:7){
  ggsave(master_plots[[i]], file = paste("Output/Figures/GSEA/PCA_plot_master",i-1,".png",sep=""))
}

# Do some stats - this doesn't help - everything is significant
output <- getSignificance(ES2, group = "cluster", fit = "ANOVA")
output <- getSignificance(ES2, group = "tissue", fit = "linear.model")

