# Script information -----------------------------------------------------------

# title: Analysis of dominant clonotypes
# author: James Favaloro
# date: 2021-04-19
# description: In this script we perform DE on dominant clonotypes from pt43 &
# pt 63

# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'Figure 1' script")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(scRepertoire)
  library(pheatmap)
  library(escape)
  library(dittoSeq)
  library(DoMultiBarHeatmap)
})

# 2. Load data, define functions and colours -------------------------------------------------------

# Load data

integrated.list = readRDS("Data/integrated.list.rds")
black_list = readRDS("Data/black_list.rds")
annotations = read.csv("Data/annotations.csv")

# Define functions

# Define %!in% function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Define colours

# Define colour scheme for scRepertoire
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
# Define colour scheme for Seurat
seurat_colours = c("0" = "#F8766D", "1" = "#D39200", "2" = "#93AA00", "3" = "#00BA38", "4" = "#00C19F", "5" = "#619CFF", "6" = "#FF61C3")

# 3. Pre-processing  -------------------------------------------------------------------------------

# Pull out the seurat object to analyse
x = integrated.list$all_ds.integrated

# Set the cluster resolution back to something more reasonable
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"
DefaultAssay(x) = "RNA"

# Fix clonotype order - not sure why but current code orders single and small in reverse(?)
slot(x, "meta.data")$cloneType <- factor(slot(x, "meta.data")$cloneType, 
                                         levels = c("Hyperexpanded (628 < X <= 6275)", 
                                                    "Large (63 < X <= 628)", 
                                                    "Medium (6 < X <= 63)", 
                                                    "Small (1 < X <= 6)", 
                                                    "Single (0 < X <= 1)", NA))

# Backup x and remove the integrated list to free up memory
x_backup = x
rm(integrated.list)

# Split the seurat object into BM and PB
x_bm = subset(x, subset = tissue == "bm")
x_pb = subset(x, subset = tissue == "pb")

# Create directories if required
dir.create("Output/Figures/General/Seurat/")
dir.create("Output/Figures/General/scRepertoire/")
dir.create("Output/Figures/General/Vlnplots/")
dir.create("Output/Figures/General/Ridgeplots/")
dir.create("Output/Figures/General/Featureplots/")


# 4. User input ------------------------------------------------------------------------------------

# Store genes we're interested in plotting as a character vector we can iterate over
# Follow the convention: vector_of_genes = c("GENE1", "GENE2", "GENE3)

exhaustion = c("PDCD1", "HAVCR2", "LAG3", "CD244", "ENTPD1", "TIGIT", "CD160", "BTLA")
activation = c("GZMA")
effector = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "EOMES", "TBET")
regulators = c("EOMES", "TBX21", "TCF7")
residency = c("CD69", "KLF2", "KLF3")

#Set genes of interest (GOI) to whichever set of genes we want to plot
# i.e. if you want to plot effector genes: GOI = effector
GOI = effector

# Set Subset of interest (SOI) to whichever subset of cells we want to plot # Not yet implemented
SOI = x


# 5. Analyse  --------------------------------------------------------------------------------------

# Collate dominant clonotypes
clonotypes = table(x$CTstrict, x$orig.ident) %>% as.data.frame.matrix() %>% rownames_to_column()

# Subset dominant clonotype
x_dom = subset(x, subset = CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT")


# Perform DE testing of differences observed between bm and pb of most dominant
# clonotypes in the dataset - hash out as required.

# 1. Most dominant in pb43 & bm43
x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT" & x$orig.ident == "bm43"] = "bm"
# vs. dominant clone in the pb
x$DE_group[x$CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT" & x$orig.ident == "pb43"] = "pb"
# vs. all other cluster 1 clones in the bm
x$DE_group[!x$CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT" & x$orig.ident == "bm43" & x$seurat_clusters == "1"] = "other_bm"

# Set Idents to DE_group,  calculate average expression of genes & perform DE testing
Idents(x) <- "DE_group"
avg.x <- log1p(AverageExpression(x, verbose = FALSE)$RNA)
avg.x = as.data.frame(avg.x)
avg.x$gene <- rownames(avg.x)

# vs. dominant clones in the pb
DE_markers = FindMarkers(x, ident.1 = "bm", ident.2 = "pb")

# vs. all other cluster 1 clones in the bm
DE_markers = FindMarkers(x, ident.1 = "bm", ident.2 = "other_bm")

DE_markers$gene = rownames(DE_markers)

# Remove black_list genes from DE_markers

DE_markers = subset(DE_markers, subset = DE_markers$gene %!in% black_list)

# Pull out genes to annotate

bm_genes = DE_markers[order(DE_markers$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
bm_genes = as.vector(bm_genes$gene)

pb_genes = DE_markers[order(DE_markers$avg_log2FC, decreasing = FALSE), ] %>% head(n = 20)
pb_genes = as.vector(pb_genes$gene)

genes = c(bm_genes, pb_genes)

# Graph the Average expression of genes and label the top 20 differentially expressed between the BM and PB
p1 = ggplot(avg.x, aes(bm,other_bm)) + geom_point(size = 0.1) + ggtitle("DE genes - vs. other cluster 1 cells in the BM")
p1 = LabelPoints(plot = p1, points = bm_genes, repel = TRUE, color = "red")
p1 = LabelPoints(plot = p1, points = pb_genes, repel = TRUE, color = "blue")
plot(p1)
ggsave(filename = paste0("Output/Figures/DE/Dominant_clonotype_bm_vs_other_cluster_1_bm.png"), device = "png")

# Look at data as a split bubble plot
DotPlot(x_dom, features = genes, assay = "SCT") + theme(axis.text.x = element_text(angle = 90))

# Prep some Ridge plots for exhaustion/activation markers etc.
# Copy code from Feature script



GOI = regulators

RidgePlot(x_dom, "PDCD1", cols = c("red", "blue"))
RidgePlot(x_dom, "TIGIT", cols = c("red", "blue"))
RidgePlot(x_dom, "HAVCR2", cols = c("red", "blue")) #TIM3
RidgePlot(x_dom, "LAG3", cols = c("red", "blue"))
RidgePlot(x_dom, "CD160", cols = c("red", "blue"))
RidgePlot(x_dom, "CD244", cols = c("red", "blue"))

exhaustion = c("PDCD1", "HAVCR2", "LAG3", "CD244", "ENTPD1", "TIGIT", "CD160", "BTLA")



# Ridge plots
for (i in 1:length(GOI)) {
  temp = (RidgePlot(x_dom, features = GOI[`i`], cols = c("red", "blue")))
  ggsave(plot = temp, filename = paste0("Output/Figures/General/Ridgeplots/Ridge_", GOI[i], ".png"))
}

# List out top 5 dominant clonotypes by sample
# Top line is frequency in sample compared to paired sample

# bm 43 1. 335:1308, 2. 162:713, 3. 148:5, 4. 77:246, 5. 76:1
# TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT	
# TRAV21.TRAJ49.TRAC_TGTGCTGAGGGCACCGGTAACCAGTTCTATTTT_TRBV6-5.TRBJ2-7.None.TRBC2_TGTGCCAGCAGGGGGACAGGGCCTCCCTACGAGCAGTACTTC
# TRAV8-6.TRAJ20.TRAC_TGTGCTGTGAGGGTCGACTACAAGCTCAGCTTT_TRBV7-3.TRBJ1-2.None.TRBC1_TGTGCCAGCAGCCCTAATCCAGGAGATAACTATGGCTACACCTTC	
# TRAV1-2.TRAJ20.TRAC_TGTGCTCTCAGGGGCGACTACAAGCTCAGCTTT_TRBV24-1.TRBJ2-1.None.TRBC2_TGTGCCACCAGTGATTCCCCCCGGACTAGCGGGAACAATGAGCAGTTCTTC
# TRAV13-1.TRAJ20.TRAC_TGTGCAGCATACGACTACAAGCTCAGCTTT_TRBV19.TRBJ2-3.None.TRBC2_TGTGCCAGTACAGGGGGACGCACAGATACGCAGTATTTT	

# bm 63 1. 246:45, 2. 230:202, 179:17, 4. 125:724, 5. 99:88 
# TRAV19.TRAJ41.TRAC_TGTGCTCTGAGGACCCCGACGAACTCAAATTCCGGGTATGCACTCAACTTC_TRBV20-1.TRBJ2-1.None.TRBC2_TGCAGTGCTAGTGGACCGGCGGAGATCAATGAGCAGTTCTTC	
# TRAV12-1.TRAJ20.TRAC_TGTGTGGTGAACTGGAGATCTAACGACTACAAGCTCAGCTTT_TRBV28.TRBJ2-3.None.TRBC2_TGTGCCAGCAGTTTTCCTAGCGGGGGGGTTAGTACAGATACGCAGTATTTT
# TRAV12-1.TRAJ43.TRAC_TGTGTGGTGGGCGATTGGTTCGGGGACATGCGCTTT_TRBV20-1.TRBJ2-1.None.TRBC2_TGCAGTGCCCTAAAGCCGGGGACGAGCTCCTACAATGAGCAGTTCTTC	
# TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC	
# TRAV14/DV4.TRAJ39.TRAC_TGTGCAATGAGAGAGGTGGAAAGCAACATGCTCACCTTT_TRBV7-6.TRBJ2-7.None.TRBC2_TGTGCCAGCAGCACCTTCTCCTACGAGCAGTACTTC	

# pb43 1. 1308:335, 2. 713:162, 3. 246:77, 4. 235:69, 5. 187:65
# TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT	
# TRAV21.TRAJ49.TRAC_TGTGCTGAGGGCACCGGTAACCAGTTCTATTTT_TRBV6-5.TRBJ2-7.None.TRBC2_TGTGCCAGCAGGGGGACAGGGCCTCCCTACGAGCAGTACTTC
# TRAV1-2.TRAJ20.TRAC_TGTGCTCTCAGGGGCGACTACAAGCTCAGCTTT_TRBV24-1.TRBJ2-1.None.TRBC2_TGTGCCACCAGTGATTCCCCCCGGACTAGCGGGAACAATGAGCAGTTCTTC	
# NA_NA_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT	
# TRAV17.TRAJ29.TRAC_TGTGCTACGGACGCGGGTAGTTCAGGAAACACACCTCTTGTCTTT_TRBV4-1.TRBJ2-1.None.TRBC2_TGCGCCAGCAGCCAGGCGGGTAGGGGGACAACCTACAATGAGCAGTTCTTC	

# pb 63 1. 724:125, 2. 331:93, 3. 245:70, 4. 202:230, 5. 193:24
# TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC	
# TRAV23/DV6.TRAJ48.TRAC_TGTGCAGCAAGCATAGGGAACTTTGGAAATGAGAAATTAACCTTT_TRBV4-3.TRBJ1-1.None.TRBC1_TGCGCCAGCAGCCCTTCGAGGAACACTGAAGCTTTCTTT
# TRAV13-1.TRAJ48.TRAC_TGTGCAGCTCCCGTGGATTTTGGAAATGAGAAATTAACCTTT_TRBV28.TRBJ1-3.None.TRBC1_TGTGCCAGCAGCCAAATCTCAGGGGGAAACACCATATATTTT	
# TRAV12-1.TRAJ20.TRAC_TGTGTGGTGAACTGGAGATCTAACGACTACAAGCTCAGCTTT_TRBV28.TRBJ2-3.None.TRBC2_TGTGCCAGCAGTTTTCCTAGCGGGGGGGTTAGTACAGATACGCAGTATTTT	
# TRAV19.TRAJ30.TRAC_TGTGCTCTGTGGTACGGGAGAGATGACAAGATCATCTTT_TRBV20-1.TRBJ1-1.None.TRBC1_TGCAGTGCAAAAGTTAACACTGAAGCTTTCTTT	


# Perform DE testing of differences observed between bm and pb of most dominant
# clonotypes in the dataset - hash out as required.

# 1. Most dominant in pb43 & bm43
x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT" & x$orig.ident == "bm43"] = "bm"
x$DE_group[x$CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT" & x$orig.ident == "pb43"] = "pb"

# 2. 2nd most dominant in pt43

x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV21.TRAJ49.TRAC_TGTGCTGAGGGCACCGGTAACCAGTTCTATTTT_TRBV6-5.TRBJ2-7.None.TRBC2_TGTGCCAGCAGGGGGACAGGGCCTCCCTACGAGCAGTACTTC" & x$orig.ident == "bm43"] = "bm"
x$DE_group[x$CTstrict == "TRAV21.TRAJ49.TRAC_TGTGCTGAGGGCACCGGTAACCAGTTCTATTTT_TRBV6-5.TRBJ2-7.None.TRBC2_TGTGCCAGCAGGGGGACAGGGCCTCCCTACGAGCAGTACTTC" & x$orig.ident == "pb43"] = "pb"

# 3. Most dominant in bm63

x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV19.TRAJ41.TRAC_TGTGCTCTGAGGACCCCGACGAACTCAAATTCCGGGTATGCACTCAACTTC_TRBV20-1.TRBJ2-1.None.TRBC2_TGCAGTGCTAGTGGACCGGCGGAGATCAATGAGCAGTTCTTC" & x$orig.ident == "bm63"] = "bm"
x$DE_group[x$CTstrict == "TRAV19.TRAJ41.TRAC_TGTGCTCTGAGGACCCCGACGAACTCAAATTCCGGGTATGCACTCAACTTC_TRBV20-1.TRBJ2-1.None.TRBC2_TGCAGTGCTAGTGGACCGGCGGAGATCAATGAGCAGTTCTTC" & x$orig.ident == "pb63"] = "pb"

# 4. 2nd most dominant in bm63

x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV12-1.TRAJ20.TRAC_TGTGTGGTGAACTGGAGATCTAACGACTACAAGCTCAGCTTT_TRBV28.TRBJ2-3.None.TRBC2_TGTGCCAGCAGTTTTCCTAGCGGGGGGGTTAGTACAGATACGCAGTATTTT" & x$orig.ident == "bm63"] = "bm"
x$DE_group[x$CTstrict == "TRAV12-1.TRAJ20.TRAC_TGTGTGGTGAACTGGAGATCTAACGACTACAAGCTCAGCTTT_TRBV28.TRBJ2-3.None.TRBC2_TGTGCCAGCAGTTTTCCTAGCGGGGGGGTTAGTACAGATACGCAGTATTTT" & x$orig.ident == "pb63"] = "pb"

# 5. Most dominant in pb63

x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC" & x$orig.ident == "bm63"] = "bm"
x$DE_group[x$CTstrict == "TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC" & x$orig.ident == "pb63"] = "pb"

# 6. 2nd most dominant in pb63

x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV23/DV6.TRAJ48.TRAC_TGTGCAGCAAGCATAGGGAACTTTGGAAATGAGAAATTAACCTTT_TRBV4-3.TRBJ1-1.None.TRBC1_TGCGCCAGCAGCCCTTCGAGGAACACTGAAGCTTTCTTT" & x$orig.ident == "bm63"] = "bm"
x$DE_group[x$CTstrict == "TRAV23/DV6.TRAJ48.TRAC_TGTGCAGCAAGCATAGGGAACTTTGGAAATGAGAAATTAACCTTT_TRBV4-3.TRBJ1-1.None.TRBC1_TGCGCCAGCAGCCCTTCGAGGAACACTGAAGCTTTCTTT" & x$orig.ident == "pb63"] = "pb"

# Set Idents to DE_group,  calculate average expression of genes & perform DE testing
Idents(x) <- "DE_group"
avg.x <- log1p(AverageExpression(x, verbose = FALSE)$RNA)
avg.x = as.data.frame(avg.x)
avg.x$gene <- rownames(avg.x)

DE_markers = FindMarkers(x, ident.1 = "bm", ident.2 = "pb")
DE_markers$gene = rownames(DE_markers)

# Remove black_list genes from DE_markers

DE_markers = subset(DE_markers, subset = DE_markers$gene %!in% black_list)

# Pull out genes to annotate

bm_genes = DE_markers[order(DE_markers$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
bm_genes = as.vector(bm_genes$gene)

pb_genes = DE_markers[order(DE_markers$avg_log2FC, decreasing = FALSE), ] %>% head(n = 20)
pb_genes = as.vector(pb_genes$gene)

# Graph the Average expression of genes and label the top 20 differentially expressed between the BM and PB
p1 = ggplot(avg.x, aes(bm,pb)) + geom_point(size = 0.1) + ggtitle("Global differences")
p1 = LabelPoints(plot = p1, points = bm_genes, repel = TRUE, color = "red")
p1 = LabelPoints(plot = p1, points = pb_genes, repel = TRUE, color = "blue")
plot(p1)

# Highlight clonotypes

# Collate dominant clonotypes
clonotypes = table(x$CTstrict, x$orig.ident) %>% as.data.frame.matrix() %>% rownames_to_column()

# Demonstrate clonal distribution using an alluvial graph highlighting the most dominant clonotype per patient (for now)
dom_pt43 = c("TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT", "TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC")

dom_pt63 = c("TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC")

# Output figure 3.5

alluvialClonotypes(x, cloneCall ="gene+nt", 
                   y.axes = c("patient", "tissue", "seurat_clusters"), 
                   color = c("TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT", "TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC")) + 
  scale_fill_manual(values = c("grey", colorblind_vector(1)))

# Try just using orig.ident
x$DE_group = NA
x$DE_group[x$CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT"] = "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT"
x$DE_group[x$CTstrict == "TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC"] = "TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC"

clones_to_plot = c("TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT", "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT")

alluvialClonotypes(x, cloneCall ="gene+nt", 
                   y.axes = c("orig.ident", "seurat_clusters"), 
                   color = clones_to_plot) +
  scale_fill_manual(values = c("grey", colorblind_vector(2)))
