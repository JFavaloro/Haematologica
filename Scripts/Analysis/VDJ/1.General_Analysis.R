# Script information -------------------------------------------------------------------------------

# Title: V(D)J - General analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will use Immunarch (https://immunarch.com/) and scRepertoire 
# (https://github.com/ncborcherding/scRepertoire) to assess clonality and diversity of recovered 
# clonotypes using the vignette: https://ncborcherding.github.io/vignettes/vignette.html. For the
# general analysis, figures will be output on the basis of raw data. Data can be viewed with
# downsampling considered by find-replace "all_vdj" for "all_vdj_ds" for section 4.1.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '1. V(D)J - General analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, immunarch, VennDiagram, scales, circlize, Cairo, ggcorrplot, ggraph)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed 10x data and extract all_10x object for further analysis
all_list <- readRDS("Data/R_out/10x/all_list.rds")
all_10x <- all_list$all_10x

# Load V(D)J data in scRepertoire format
all_vdj <- readRDS("Data/R_out/10x/all_vdj.rds")

# Load V(D)J data in immunarch format #NB These are renamed filtered_contig_annotation *.CSV files
immdata <- repLoad("Data/10x/TCR", .mode = "paired")

# These are the colours to use for 4x clonal bins
plot_cols <- c("Expanded" = "#FF4B20", "Large" = "#FFB433", 
               "Medium" = "#7AC5FF", "Small" = "#0348A6")


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 Downsample scRepertoire data to lowest number of observable cells (BM031; 6262) --------------

# Create empty list for downsampled data
all_vdj_ds <- list()

# Downsample all cells, excluding BM31 (3rd item of list) (lowest common denominator)
for (i in 1:length(all_vdj)) {
  temp <- head(all_vdj[[i]], 6262)
  all_vdj_ds[[i]] <- temp
}

# Rename the items in the list
names(all_vdj_ds) <- names(all_vdj)


# 3.2 Split all_vdj_ds into lists of BM and PB
all_vdj_ds_BM <- list(all_vdj_ds$BM13, all_vdj_ds$BM31, all_vdj_ds$BM43, all_vdj_ds$BM63)
names(all_vdj_ds_BM) <- names(all_vdj_ds)[c(1,3,5,7)]
all_vdj_ds_PB <- list(all_vdj_ds$PB13, all_vdj_ds$PB31, all_vdj_ds$PB43, all_vdj_ds$PB63)
names(all_vdj_ds_PB) <- names(all_vdj_ds)[c(2,4,6,8)]


# 3.3 Downsample scRepertoire data to lowest number of observable clones (PB063; 1,378) ------------
immdata_ds <- repSample(immdata$data, .method = "sample", .prob = TRUE, .n = 1378)


# 3.4 Split immdata_ds into lists of BM and PB -----------------------------------------------------
immdata_ds_BM <- list(immdata_ds$BM13, immdata_ds$BM31, immdata_ds$BM43, immdata_ds$BM63)
names(immdata_ds_BM) <- names(immdata_ds)[1:4]
immdata_ds_PB <- list(immdata_ds$PB13, immdata_ds$PB31, immdata_ds$PB43, immdata_ds$PB63)
names(immdata_ds_PB) <- names(immdata_ds)[5:8]


# 3.5 Split the all_10x object into BM and PB, rename clusters and recombine -----------------------
all_10x_BM <- subset(all_10x, subset = Tissue == "BM")
all_10x_PB <- subset(all_10x, subset = Tissue == "PB")

# Rename Idents - used for Moritisa temporarily
all_10x_BM <- 
  RenameIdents(all_10x_BM, `T[EM]` = "3 ~ T[EM] ~ BM", `T[TE]` = "4 ~ T[TE] ~ BM", `T[N]` = "1 ~ T[N] ~ BM", 
               `Cyto ~ T[EM]` = "5 ~ Cyto ~ T[EM] ~ BM", `P[RE] ~ EX` = "6 ~ P[RE] ~ EX ~ BM", `T[CM]` = "2 ~ T[CM] ~ BM")
all_10x_PB <- 
  RenameIdents(all_10x_PB, `T[EM]` = "3 ~ T[EM] ~ PB", `T[TE]` = "4 ~ T[TE] ~ PB", `T[N]` = "1 ~ T[N] ~ PB", 
               `Cyto ~ T[EM]` = "5 ~ Cyto ~ T[EM] ~ PB", `P[RE] ~ EX` = "6 ~ P[RE] ~ EX ~ PB", `T[CM]` = "2 ~ T[CM] ~ PB")

# Rename Idents
all_10x_BM <- 
  RenameIdents(all_10x_BM, `T[EM]` = "T[EM] ~ BM", `T[TE]` = "T[TE] ~ BM", `T[N]` = "T[N] ~ BM", 
               `Cyto ~ T[EM]` = "Cyto ~ T[EM] ~ BM", `P[RE] ~ EX` = "P[RE] ~ EX ~ BM", `T[CM]` = "T[CM] ~ BM")
all_10x_PB <- 
  RenameIdents(all_10x_PB, `T[EM]` = "T[EM] ~ PB", `T[TE]` = "T[TE] ~ PB", `T[N]` = "T[N] ~ PB", 
               `Cyto ~ T[EM]` = "Cyto ~ T[EM] ~ PB", `P[RE] ~ EX` = "P[RE] ~ EX ~ PB", `T[CM]` = "T[CM] ~ PB")


# Merge Seurat objects back together
all_10x_merge <- merge(all_10x_BM, all_10x_PB)

# Normalise the data again after merge
all_10x_merge <- NormalizeData(all_10x_merge)

# Fix Idents of all objects to analyse
all_10x$Cluster <- Idents(all_10x)
all_10x_BM$Cluster <- Idents(all_10x_BM)
all_10x_PB$Cluster <- Idents(all_10x_PB)
all_10x_merge$Cluster <- Idents(all_10x_merge)

# Re-level data for pretty plotting
all_10x_merge$Cluster <- 
  factor(all_10x_merge$Cluster, levels = c("T[N] ~ BM","T[N] ~ PB", "T[CM] ~ BM", "T[CM] ~ PB", 
                                "T[EM] ~ BM", "T[EM] ~ PB", "T[TE] ~ BM", "T[TE] ~ PB", 
                                "Cyto ~ T[EM] ~ BM", "Cyto ~ T[EM] ~ PB",
                                "P[RE] ~ EX ~ BM", "P[RE] ~ EX ~ PB"))

# Extract metadata for additional analyses
# NB: 20220626 This appears necessary due to changes in scRepertoire
merged_meta <- expression2List(all_10x_merge, split.by = "Cluster")

merged_meta <- merged_meta[c("T[N] ~ BM","T[N] ~ PB", "T[CM] ~ BM", "T[CM] ~ PB", 
                             "T[EM] ~ BM", "T[EM] ~ PB", "T[TE] ~ BM", "T[TE] ~ PB", 
                             "Cyto ~ T[EM] ~ BM", "Cyto ~ T[EM] ~ PB",
                             "P[RE] ~ EX ~ BM", "P[RE] ~ EX ~ PB")]


# 4. Analyse ---------------------------------------------------------------------------------------

# 4.1 General "pseudobulk" analysis ----------------------------------------------------------------

# 4.1.A Quantify Clonotypes ------------------------------------------------------------------------

# Calculate the number of unique clonotypes output figure for all samples and paired samples
quantContig(all_vdj, cloneCall="gene+nt", scale = T)
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Unique_clones_sample.png"))
quantContig(all_vdj, cloneCall="gene+nt", scale = T, group = "sample")
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Unique_clones_patient.png"))


# 4.1.B Clonotype Abundance  -----------------------------------------------------------------------

# Calculate the relative abundance of clonotypes output figure for all samples and paired samples
abundanceContig(all_vdj, cloneCall = "gene", scale = F)
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Clonal_abundance_sample.png"))
abundanceContig(all_vdj, cloneCall = "gene", scale = F, group = "sample")
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Clonal_abundance_patient.png"))


# 4.1.C Length of Clonotypes -----------------------------------------------------------------------

# Calculate the length of the CDR3 region in both single cell and "pseudobulk" and output figures
lengthContig(all_vdj, cloneCall = "aa", chains = "combined")
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/CDR3_length_sc.png"))
lengthContig(all_vdj, cloneCall = "aa", chains = "single")
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/CDR3_length_bulk.png"))


# 4.1.D Clonal Space Homeostasis -------------------------------------------------------------------

# Extract the clonal bins as a proportion of total
temp <- clonalHomeostasis(all_vdj, cloneCall = "gene+nt", 
                         cloneTypes = c(Small = 0.001, Medium = 0.01, 
                                        Large = 0.1, Expanded = 1))
temp[["data"]][["value"]] <- temp[["data"]][["value"]]*100
temp[["data"]][["Var2"]] <- 
  as.factor(c(rep("Small" , 8) , rep("Medium" , 8), rep("Large" , 8) , rep("Expanded" , 8)))

# Move ggplot2 data into a dataframe to plot
temp <- as.data.frame(temp$data)

# Plot data
ggplot(temp, 
       aes(fill = Var2, y = value, x = Var1))  + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = rev(plot_cols), name = "Clonal bin", 
                    labels = c("Small", "Medium", "Large", "Expanded")) +
  labs(x = "Sample", y = "Relative Abundance (%)", 
       title = "Distribution of clonal bins") + theme_bw() + 
  theme(text=element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Distribution_of_clonal_bins.png"), 
       units = "cm", height = 30, width = 18, bg = "transparent")


# 4.2 Clonal overlap analysis ----------------------------------------------------------------------

# 4.2.A Determine the level of clonal overlap between paired samples -------------------------------

# Create an empty list of results
overlap <- list()

# Run a loop to determine overlap between paired samples
for (i in 1:4) {
  for (j in c(10,20,50,100,1378)) {
    x <- head(immdata_ds_BM[[`i`]], n = j) %>% select(CDR3.aa, V.name, J.name)
    y <- head(immdata_ds_PB[[`i`]], n = j) %>% select(CDR3.aa, V.name, J.name)
    temp <- intersect(x,y)
    overlap[[length(overlap) +1]] <- temp
  }
} 

# 4.2.A.i Plot results as Stacked barplots - this is easier through Prisim - Done

# 4.2.A.ii Plot results as Venn Diagrams:

# Create an empty list to store data
venn = list()

# Run a loop to extract the top 10 clones from all samples
for (i in names(immdata_ds)) {
  for (j in c(10,20,50,100,1378)){
    temp <- head(immdata$data[[`i`]][["Sequence"]], n = j)
    venn[[length(venn) +1]] <- temp
  }
}

# Fix the names for the venn list
names(venn) <- apply(expand.grid(c("10", "20", "50", "100", "1378"), 
                                c("BM13", "BM31", "BM43", "BM63", 
                                  "PB13", "PB31", "PB43", "PB63")),1, paste, collapse = "_")

# Split the venn list into BM and PB lists
venn_BM <- venn[1:20]
venn_PB <- venn[21:40]

# Run a loop to produce Venn Diagrams of clonal overlap
for (i in 1:20) {
  venn.diagram(x = list(venn_BM[[`i`]], venn_PB[[`i`]]),
               category.names = c("BM", "PB"), imagetype = "svg", height = 10, width = 10, resolution = 500,
               filename = paste0("Output/Figures/VDJ/10x/Final/Venn/", names(venn)[[`i`]], ".svg"),
               output = TRUE, disable.logging = TRUE,
               col=c("blue", 'red'),
               fill = c('blue', "red"), print.mode = c('percent', 'raw'))
}


# 4.2B Clonal overlap of dominant clones # NB: This function will take the top n from both samples.
# If overlap is shared, they will be ignored and therefore n needs to be increased accordingly.

temp <- compareClonotypes(all_vdj, cloneCall = "gene+nt", samples = c("BM13", "PB13"), 
                          numbers = 12) +  NoLegend()
temp[["data"]][["Proportion"]] <- temp[["data"]][["Proportion"]]*100
temp$labels$y <- "Proportion of all clonotypes (%)"
temp + theme(text=element_text(size=24), panel.border = element_blank()) + 
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  scale_x_discrete(expand = c(0.1,0.1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/PT13_Top_10_overlap.png"), 
       units = "cm", height = 25, width = 20, bg = "transparent")

temp <- compareClonotypes(all_vdj, cloneCall = "gene+nt", samples = c("BM31", "PB31"), 
                          numbers = 10) +  NoLegend()
temp[["data"]][["Proportion"]] <- temp[["data"]][["Proportion"]]*100
temp$labels$y <- "Proportion of all clonotypes (%)"
temp + theme(text=element_text(size=24)) + 
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  scale_x_discrete(expand = c(0.1,0.1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/PT31_Top_10_overlap.png"), 
       units = "cm", height = 25, width = 20, bg = "transparent")

temp <- compareClonotypes(all_vdj, cloneCall = "gene+nt", samples = c("BM43", "PB43"), 
                          numbers = 12) +  NoLegend()
temp[["data"]][["Proportion"]] <- temp[["data"]][["Proportion"]]*100
temp$labels$y <- "Proportion of all clonotypes (%)"
temp + theme(text=element_text(size=24)) + 
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  scale_x_discrete(expand = c(0.1,0.1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/PT43_Top_10_overlap.png"), 
       units = "cm", height = 25, width = 20, bg = "transparent")

temp <- compareClonotypes(all_vdj, cloneCall = "gene+nt", samples = c("BM63", "PB63"), 
                          numbers = 11) +  NoLegend()
temp[["data"]][["Proportion"]] <- temp[["data"]][["Proportion"]]*100
temp$labels$y <- "Proportion of all clonotypes (%)"
temp + theme(text=element_text(size=24)) + 
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  scale_x_discrete(expand = c(0.1,0.1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/PT63_Top_10_overlap.png"), 
       units = "cm", height = 25, width = 20, bg = "transparent")


# 4.2C Calculate the morisita index and output figure for paired samples ---------------------------

# Plot overlap across samples
clonalOverlap(all_vdj, cloneCall = "gene+nt", method = "morisita") + 
  theme_bw(base_size = 18) + ggtitle(label = "Morisita Index by Sample") + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5)) + 
  labs(x = NULL, y = NULL) + NoLegend() + theme(axis.text.x=element_text(angle = 45, hjust = 0.8))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Morisita.png"), 
       units = "cm", height = 24, width = 20, bg = "transparent")

# Plot overlap across cluster and tissue
clonalOverlap(merged_meta, cloneCall = "strict", method = "morisita", chain = "both") + 
  theme_bw(base_size = 18) + ggtitle(label = "Morisita Index across Cluster and Tissue") + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5)) + 
  labs(x = NULL, y = NULL) + theme(axis.text.x=element_text(angle = 45, hjust = 0.8)) +
  scale_x_discrete(labels = label_parse()) + scale_y_discrete(labels = label_parse()) + NoLegend()
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Morisita_cluster_tissue.png"), 
       units = "cm", height = 24, width = 38, bg = "transparent")


# 4.2D Plot Clonal Network

# Plot BM interactions
clonalNetwork(all_10x_BM, cloneCall = "strict") + 
  theme_bw() + ggtitle(label = "Clonal Network across BM") + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5)) + 
  NoLegend()
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/BM_Network.png"), 
       units = "cm", height = 25, width = 25, bg = "transparent")

# Plot PB interactions
clonalNetwork(all_10x_PB, cloneCall = "strict") + 
  theme_bw() + ggtitle(label = "Clonal Network across PB") + 
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5)) + 
  NoLegend()
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/PB_Network.png"), 
       units = "cm", height = 25, width = 25, bg = "transparent")


# 4.3 Clonal diversity analysis --------------------------------------------------------------------

# 4.3.A Determine the diversity level of clones across all samples ---------------------------------

Diverity_metrics_sample <- 
  clonalDiversity(all_vdj, cloneCall = "gene+nt", n.boots = 100, exportTable = T)
Diverity_metrics_cluster <- 
  clonalDiversity(merged_meta, cloneCall = "gene+nt", n.boots = 100, exportTable = T)
clonalDiversity(all_vdj, cloneCall = "gene+nt", n.boots = 100)
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Diversity.png"))

# Plot the Inverse Simpson Index per sample
plot <- Diverity_metrics_sample[c(2,6)]
plot$Group <- 
  factor(plot$Group, levels = c("BM13", "PB13", "BM31", "PB31", "BM43", "PB43", "BM63", "PB63"))

ggplot(plot, aes(x = Group, y = Inv.Simpson)) + 
  geom_bar(stat = "identity", fill=rep(c("blue","red"), 4)) + 
  theme_bw(base_size = 18) + ggtitle(label = "Inverse Simpson by Sample") +
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5)) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 0.8))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Inverse_Simpson.png"), 
       units = "cm", height = 30, width = 16, bg = "transparent")

# Plot the Inverse Simpson Index per cluster
plot <- Diverity_metrics_cluster[c(2,6)]
plot$Group <- gsub("_", " ~ ", plot$Group)
plot$Group <- 
  factor(plot$Group, levels = c("T[N] ~ BM","T[N] ~ PB", "T[CM] ~ BM", "T[CM] ~ PB", 
                                "T[EM] ~ BM", "T[EM] ~ PB", "T[TE] ~ BM", "T[TE] ~ PB", 
                                "Cyto ~ T[EM] ~ BM", "Cyto ~ T[EM] ~ PB",
                                "P[RE] ~ EX ~ BM", "P[RE] ~ EX ~ PB"))

ggplot(plot, aes(x = Group, y = Inv.Simpson)) + 
  geom_bar(stat = "identity", fill=rep(c("blue","red"), 6)) + 
  theme_bw(base_size = 18) +  ggtitle(label = "Inverse Simpson by Cluster") +
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5)) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + scale_x_discrete(labels = label_parse())
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Inverse_Simpson_by_cluster.png"), 
       units = "cm", height = 30, width = 24, bg = "transparent")


# 4.4 Gene usage analysis --------------------------------------------------------------------------

# N/A - This data is available directly from vloupe


# 4.5 General single cell analysis -----------------------------------------------------------------

# Plot the clonotype data over UMAP
DimPlot(all_10x, group.by = "cloneType", cols = c("#FF4B20", "#FFB433", "#7AC5FF", "#0348A6"), 
        na.value = "grey") +  ggtitle(label = "Distribution of clonal bins") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + labs(color = "Clonal bin") +
  theme(legend.text = element_text(size = 18), plot.title = element_text(size = 18, hjust = 0.5), 
        legend.title = element_text(size = 18), legend.text.align = 0)
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Distribution_of_clonal_bins_UMAP.png"), 
       units = "cm", height = 20, width = 25, bg = "transparent")

# Quantify clonal bins seperated by tissue and cluster
occupiedscRepertoire(all_10x, x.axis = "Cluster") + 
  scale_x_discrete(labels = label_parse()) + theme_bw() +
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Distribution_of_clonal_bins_barplot.png"), 
       units = "cm", height = 20, width = 10, bg = "transparent")

# Quantify clonal bins across the BM
occupiedscRepertoire(all_10x_BM, x.axis = "Cluster") +
  scale_x_discrete(labels = label_parse()) + theme_bw() +
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Distribution_of_clonal_bins_BM.png"), 
       units = "cm", height = 20, width = 10, bg = "transparent")

# Quantify clonal bins across the PB
occupiedscRepertoire(all_10x_PB, x.axis = "Cluster") + 
  scale_x_discrete(labels = label_parse()) + theme_bw() +
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Distribution_of_clonal_bins_PB.png"), 
       units = "cm", height = 20, width = 10, bg = "transparent")

# Output Startrac diversity
StartracDiversity(all_10x, type = "Tissue", 
                  sample = "Patient", by = "overall") + 
  scale_x_discrete(labels = label_parse()) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/Startrac.png"))


# 5. Save out data and print stats -----------------------------------------------------------------

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '1. V(D)J - General analysis' script")