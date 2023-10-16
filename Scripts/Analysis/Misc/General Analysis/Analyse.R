# Script information -------------------------------------------------------------------------------

# title: General Analysis
# author: James Favaloro
# date: 2021-12-19
# description: In this script we list the analyses we can perform on cleaned data

# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'General analysis' script")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(scRepertoire)
  library(pheatmap)
  library(escape)
  library(dittoSeq)
})

rmarkdown::render("Scripts/General Analysis/Analyse.R", output_format = c("html_document", "word_document"))
knitr::opts_current$set(eval = FALSE, echo = TRUE)

# 2. Load data -----------------------------------------------------------------

integrated.list = readRDS("Data/integrated.list.rds")
black_list = readRDS("black_list")



# 3. Analyse  ------------------------------------------------------------------

# Let's focus on the combined downsampled dataset. This allows us to query
# both intra tissue and inter patient differences. Downsampling means we can
# utilise the clonotype data correctly standardising things.
# Have hashed out saving for now and manually save each image
# Will work with Walter to tidy this up so that analysis can be automated
#%>% ggsave(filename = "Figures/1.1_Dimplot.svg", device = "svg", dpi = "print", units = "cm", width =, height =)

# Save single UMAP as 750x750
# 2x2 UMAP as 1225x750
# 4x2 UMAP as 2350x750

# pull out the seurat object to analyse
x = integrated.list$all_ds.integrated

# Set the cluster resolution back to something more reasonable
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"

# 3.1 Look at cluster disribution
# Output images for figure 1 - cluster distribution and quantification

# Output figure 1.1 - clusters
DimPlot(x)

# Output figure 1.2 - clusters split by/ grouped by tissue
DimPlot(x, split.by = "tissue")
DimPlot(x, group.by = "tissue")

# Output figure 1.3 - clusters split by/ grouped by sample
DimPlot(x, split.by = "orig.ident")
DimPlot(x, group.by = "orig.ident")

# Output figure 1.4 - cluster quantification split by tissue as heatmap/barplot
y = table(x@active.ident, x@meta.data$tissue)
pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 20, number_format = "%.0f")
p_col = rainbow(7)
barplot(y, col = p_col, legend.text = TRUE)

# Output figure 1.5 - cluster quantification split by sample as heatmap/barplot
y = table(x@active.ident, x@meta.data$orig.ident)
pheatmap(y, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 20, number_format = "%.0f")
barplot(y, col = p_col, legend.text = TRUE)

#3.2 Run DE to determine what clusters are
# Maybe we should do this first so we can annotate cluster plots?

# Set correct assay for DE testing
DefaultAssay(x) = "RNA"

# Run FindAllMarkers to determine genes that define our clusters
markers = FindAllMarkers(x)

# Annotate the clusters using annotation sourced from WHERE
# NB - this will pull out genes we have excluded from clustering
annotations = read.csv("annotation.csv")
annotated_markers <- inner_join(x = markers, 
                                   y = annotations[, c("gene_name", "description")],
                                   by = c("gene" = "gene_name")) %>%
  unique()

# Pull out the top 10 genes (based on level of expression) that define each cluster so we can visualise them

top_10_markers = annotated_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique = unique(top_10_markers$gene)

# Remove genes that aren't in the SCT data
# Define %!in% function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Remove genes that don't exist in scaled data
top_10_unique_blacklisted = top_10_unique[top_10_unique %!in% black_list]

# Output figure 2.1 - Genes that define clusters in a heatmap
DoHeatmap(x, assay = "SCT", features = top_10_unique_blacklisted)

# Output figure 2.2 - Genes that define clusters 
DotPlot(x, assay = "SCT", features = top_10_unique_blacklisted, split.by = "tissue", cols = c("red", "blue")) + theme(axis.text.x = element_text(angle = 90))

# Use 'knowledge' of biology to name our clusters
# Hash out for now because no knowledge

#new.cluster.ids <- c("Stressed", "TEM", "TCM", "KLRB1", "CCL4", "TN")
# May need to reconsider ID of 2 and 6 - based on clonal diversity it appears
# Cluster 2 is more diverse
#names(new.cluster.ids) <- levels(x)
#x <- RenameIdents(x, new.cluster.ids)

# 3.3 Have a look at clonotypic distribution

# Fix clonotype order - not sure why we need to do this but current code orders
# Single and small in reverse(?)
# BUT DOES BINNING ON A PER SAMPLE LEVEL NOT A PER PATIENT LEVEL!
# May have to merge seurat objects before integrating to fix...

slot(x, "meta.data")$cloneType <- factor(slot(x, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (628 < X <= 6275)", 
                                                         "Large (63 < X <= 628)", 
                                                         "Medium (6 < X <= 63)", 
                                                         "Small (1 < X <= 6)", 
                                                         "Single (0 < X <= 1)", NA))

# Set default colour scheme
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))
# Visualise distribution of clonotypes on per patiet basis
DimPlot(x, group.by = "patient")  +
  scale_color_manual(values=colorblind_vector(2))

# Output figure 3.1 - expansion of clonotypes
# Have used NoLegend for time being until I can figure out how to customise the labels
DimPlot(x, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + NoLegend()

# Output figure 3.2 - expansion of clonotypes by tissue
DimPlot(x, group.by = "cloneType", split.by = "tissue") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + NoLegend()

# Output figure 3.3 - expansion of clonotypes by sample
DimPlot(x, group.by = "cloneType", split.by = "orig.ident") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + NoLegend()

# Visualise distrubtion of clonotypes - not necessary to repeat as part of the 
# TCR distribution analysis
y = table(x$cloneType, x$orig.ident)
p_col = rainbow(5)
barplot(y, col = p_col)

# Output figure 3.4 Distribution of clonotypes by cluster, as well output
# the data as a table for quantification if required
occupiedscRepertoire(x, x.axis = "cluster")
z = occupiedscRepertoire(x, x.axis = "cluster", exportTable = TRUE)


# There's a whole heap of new tools in the updated scRepertoire package :)

# Create new object to allow for metrics to be calculated based on clusters
# to further probe clonal architecture

y = expression2List(x, group = "cluster")

# Save diversity metrics based on clonotype distribtion across clusters and
# graph with prism as per previous analysis on clonotype data alone - temp
# use barplot
z = clonalDiversity(y, cloneCall = "nt", exportTable = TRUE)
barplot(z$Inv.Simpson)

# Look at clonotype overlap between clusters
clonalOverlap(y, cloneCall="aa", method="overlap")

# Demonstrate clonal distribution using an alluvial graph highlighting the most 
# dominant clonotype per patient (for now)
# Need to figure out how to extract this information - shouldn't be this hard...

# Dominant clonotype pt.43 = TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT
# Dominant clonotype pt.63 = TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC

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



# 3.4 DE testing

#x$DE_group[x$tissue == 'bm' & GetAssayData(x, slot = "data")["CD69",] > 1] = "a"

# Set a blank metadata value for DE testing and set object ident for DE testing
x$DE_group = NA
x$DE_group[x$tissue == 'bm'] = "bm"
x$DE_group[x$tissue == 'pb'] = "pb"


# Do DE testing first across tissue
# This has come from https://satijalab.org/seurat/archive/v3.1/immune_alignment.html
# Modified as appropriate - not sure what's changed (aside from Seurat switching to Log2FC)
# But the code that worked before is broken... Seemingly works now...
# Maybe a volcano plot would be better?
# NB - not currently blacklisting any genes here

# Set Idents to DE_group,  calculate average expression of genes & perform DE testing
Idents(x) <- "DE_group"
avg.x <- log1p(AverageExpression(x, verbose = FALSE)$RNA)
avg.x = as.data.frame(avg.x)
avg.x$gene <- rownames(avg.x)

DE_markers = FindAllMarkers(x)

bm_genes = DE_markers[DE_markers$cluster == 'bm',]
up_regulated_genes = bm_genes[order(bm_genes$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
up_regulated_genes = as.vector(up_regulated_genes$gene)

pb_genes = DE_markers[DE_markers$cluster == 'pb',]
down_regulated_genes = pb_genes[order(pb_genes$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
down_regulated_genes = as.vector(down_regulated_genes$gene)

# Graph the Average expression of genes and label the top 20 differentially expressed between the BM and PB
p1 = ggplot(avg.x, aes(bm,pb)) + geom_point(size = 0.1) + ggtitle("Global differences")
p1 = LabelPoints(plot = p1, points = up_regulated_genes, repel = TRUE, color = "red")
p1 = LabelPoints(plot = p1, points = down_regulated_genes, repel = TRUE, color = "blue")
plot(p1)

# Do DE testing across patients
x$DE_group = NA
x$DE_group[x$patient == '43'] = "a"
x$DE_group[x$patient == '63'] = "b"
# Set Idents to DE_group,  calculate average expression of genes & perform DE testing
Idents(x) <- "DE_group"
avg.x <- log1p(AverageExpression(x, verbose = FALSE)$RNA)
avg.x = as.data.frame(avg.x)
avg.x$gene <- rownames(avg.x)

DE_markers = FindAllMarkers(x)

a_genes = DE_markers[DE_markers$cluster == 'a',]
a_genes = a_genes[order(a_genes$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
a_genes = as.vector(a_genes$gene)

b_genes = DE_markers[DE_markers$cluster == 'b',]
b_genes = b_genes[order(b_genes$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
b_genes = as.vector(b_genes$gene)

# Graph the Average expression of genes and label the top 20 differentially expressed between the BM and PB
p1 = ggplot(avg.x, aes(a,b)) + geom_point(size = 0.1) + ggtitle("Global differences")
p1 = LabelPoints(plot = p1, points = up_regulated_genes, repel = TRUE, color = "red")
p1 = LabelPoints(plot = p1, points = down_regulated_genes, repel = TRUE, color = "blue")
plot(p1)

# Do DE testing across patients
x$DE_group = NA
x$DE_group[x$patient == '43' & x$tissue == 'pb'] = "a"
x$DE_group[x$patient == '63' & x$tissue == 'pb'] = "b"
# Set Idents to DE_group,  calculate average expression of genes & perform DE testing
Idents(x) <- "DE_group"
avg.x <- log1p(AverageExpression(x, verbose = FALSE)$RNA)
avg.x = as.data.frame(avg.x)
avg.x$gene <- rownames(avg.x)

DE_markers = FindAllMarkers(x)

a_genes = DE_markers[DE_markers$cluster == 'a',]
a_genes = a_genes[order(a_genes$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
a_genes = as.vector(a_genes$gene)

b_genes = DE_markers[DE_markers$cluster == 'b',]
b_genes = b_genes[order(b_genes$avg_log2FC, decreasing = TRUE), ] %>% head(n = 20)
b_genes = as.vector(b_genes$gene)

# Graph the Average expression of genes and label the top 20 differentially expressed between the BM and PB
p1 = ggplot(avg.x, aes(a,b)) + geom_point(size = 0.1) + ggtitle("Global differences")
p1 = LabelPoints(plot = p1, points = up_regulated_genes, repel = TRUE, color = "red")
p1 = LabelPoints(plot = p1, points = down_regulated_genes, repel = TRUE, color = "blue")
plot(p1)





# X. scGSEA

# Impoart gene sets from Molsig database
GS <- getGeneSets(library = "H")

# Run GSEA against our dataset

ES <- enrichIt(x, gene.sets = GS, groups = 1000, cores = 8)

# Add GSEA back to seurat object as metadata

x = AddMetaData(x, ES)

# Convert seurat object to SCE to allow interaction with DittoSeq

x.sce = as.SingleCellExperiment(x)

# Define colors
colors <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

# Output results as heatmap demonstrating 
dittoHeatmap(x, genes = NULL, metas = names(ES), 
             annot.by = c("seurat_clusters", "tissue"), 
             fontsize = 7, 
             cluster_cols = FALSE,
             heatmap.colors = rev(colors(50)))

# This is a mess. Try selecting only a hand full or one at a time

dittoHeatmap(x, genes = NULL, metas = c("HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS"), 
             annot.by = c("seurat_clusters", "tissue"), 
             fontsize = 7, 
             cluster_cols = FALSE,
             heatmap.colors = rev(colors(50)))

# Looks good - try a handful more and save some out.

dittoHeatmap(x, genes = NULL, metas = c("HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"), 
             annot.by = c("seurat_clusters", "tissue"), 
             fontsize = 7, 
             cluster_cols = FALSE,
             heatmap.colors = rev(colors(50)))

# Come back to this - let's see about plotting x.sce - beautiful - works

# List things' we can plot:
#[1] "HALLMARK_ADIPOGENESIS"                      "HALLMARK_ALLOGRAFT_REJECTION"               "HALLMARK_ANDROGEN_RESPONSE"                
#[4] "HALLMARK_ANGIOGENESIS"                      "HALLMARK_APICAL_JUNCTION"                   "HALLMARK_APICAL_SURFACE"                   
#[7] "HALLMARK_APOPTOSIS"                         "HALLMARK_BILE_ACID_METABOLISM"              "HALLMARK_CHOLESTEROL_HOMEOSTASIS"          
#[10] "HALLMARK_COAGULATION"                       "HALLMARK_COMPLEMENT"                        "HALLMARK_DNA_REPAIR"                       
#[13] "HALLMARK_E2F_TARGETS"                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" "HALLMARK_ESTROGEN_RESPONSE_EARLY"          
#[16] "HALLMARK_ESTROGEN_RESPONSE_LATE"            "HALLMARK_FATTY_ACID_METABOLISM"             "HALLMARK_G2M_CHECKPOINT"                   
#[19] "HALLMARK_GLYCOLYSIS"                        "HALLMARK_HEDGEHOG_SIGNALING"                "HALLMARK_HEME_METABOLISM"                  
#[22] "HALLMARK_HYPOXIA"                           "HALLMARK_IL2_STAT5_SIGNALING"               "HALLMARK_IL6_JAK_STAT3_SIGNALING"          
#[25] "HALLMARK_INFLAMMATORY_RESPONSE"             "HALLMARK_INTERFERON_ALPHA_RESPONSE"         "HALLMARK_INTERFERON_GAMMA_RESPONSE"        
#[28] "HALLMARK_KRAS_SIGNALING_DN"                 "HALLMARK_KRAS_SIGNALING_UP"                 "HALLMARK_MITOTIC_SPINDLE"                  
#[31] "HALLMARK_MTORC1_SIGNALING"                  "HALLMARK_MYC_TARGETS_V1"                    "HALLMARK_MYC_TARGETS_V2"                   
#[34] "HALLMARK_MYOGENESIS"                        "HALLMARK_NOTCH_SIGNALING"                   "HALLMARK_OXIDATIVE_PHOSPHORYLATION"        
#[37] "HALLMARK_P53_PATHWAY"                       "HALLMARK_PANCREAS_BETA_CELLS"               "HALLMARK_PEROXISOME"                       
#[40] "HALLMARK_PI3K_AKT_MTOR_SIGNALING"           "HALLMARK_PROTEIN_SECRETION"                 "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"  
#[43] "HALLMARK_SPERMATOGENESIS"                   "HALLMARK_TGF_BETA_SIGNALING"                "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          
#[46] "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"         "HALLMARK_UV_RESPONSE_DN"                    "HALLMARK_UV_RESPONSE_UP"                   
#[49] "HALLMARK_WNT_BETA_CATENIN_SIGNALING"        "HALLMARK_XENOBIOTIC_METABOLISM"   

multi_dittoPlot(x.sce, vars = c("CD8_lung_up1", "CD8_spleen_up1"), 
                group.by = "seurat_clusters", split.by = "tissue", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

# Beautiful

dittoScatterHex(x.sce, x.var = "HALLMARK_APOPTOSIS", 
                y.var = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                do.contour = TRUE,
                split.by = "tissue")  + 
  scale_fill_gradientn(colors = rev(colors(11))) 

dittoScatterHex(x.sce, x.var = x.sce$HALLMARK_XENOBIOTIC_METABOLISM, 
                y.var = x.sce$HALLMARK_WNT_BETA_CATENIN_SIGNALING, 
                do.contour = TRUE,
                split.by = "groups")  + 
  scale_fill_gradientn(colors = rev(colors(11))) 


# Do Addmodule score - from memory module needs to be a list not a vector
# Tidy this up into a seperate analysis - currently using geneset imported into TRM signature

# Extract the top 200 genes (by level of expression) defining both lund and spleen TRM
CD8_lung_up = CD8_lung[order(CD8_lung$logFC, decreasing = TRUE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames() %>% list()
CD8_lung_down = CD8_lung[order(CD8_lung$logFC, decreasing = FALSE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames() %>% list()

CD8_spleen_up = CD8_spleen[order(CD8_lung$logFC, decreasing = TRUE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames() %>% list()
CD8_spleen_down = CD8_spleen[order(CD8_lung$logFC, decreasing = FALSE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames() list()

# Convert
CD8_lung_up.list = list(CD8_lung_up)
CD8_spleen_up.list = list(CD8_spleen_up)

x = AddModuleScore(x, features = CD8_lung_up.list, name = "lung_TRM")
x = AddModuleScore(x, features = CD8_spleen_up.list, name = "spleen_TRM")

FeaturePlot(x, features = c("lung_TRM1", "CD69"), blend = TRUE, split.by = "tissue", cols = c("grey95", "red", "green"))
FeaturePlot(x, features = c("spleen_TRM1", "CD69"), blend = TRUE, split.by = "tissue", cols = c("grey95", "red", "green"))

  # Let's start by looking at interpatient heterogeneity in the pb samples

pb = integrated.list$pb.integrated

# Let's see what markers define our clusters
# This should be run on normalised RNA data - so let's do that. For some reason
# Data needs to be normalised again? Maybe not? Looks OK for now...

DefaultAssay(pb) = "RNA"

# Run Find all markers to determine genes that define our clusters

pb.markers = FindAllMarkers(pb)

# Annotate the markers we found - why does this drop when we annotate?
annotations = read.csv("annotation.csv")
annotated_pb.markers <- inner_join(x = pb.markers, 
                                 y = annotations[, c("gene_name", "description")],
                                 by = c("gene" = "gene_name")) %>%
  unique()

# Pull out the top 10 genes that define each cluster so we can visualise them

top_10_pb.markers = annotated_pb.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_pb_unique = unique(top_10_pb.markers$gene)

# Define %!in% function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Remove genes that don't exist in scaled data
top_10_unique_blacklisted = top_10_unique[top_10_unique %!in% black_list]

# Plot a heatmap of genes that define clusters - WHY IS THIS NOT WORKING with cells argument?
DoHeatmap(pb, assay = "SCT", features = top_10_unique_blacklisted)

# Do the same thing but black list some genes
x = pb
DefaultAssay(x) = "RNA"
counts <- GetAssayData(x, assay = "RNA")
black_list = grep(pattern = "^MT-|RPS|RPL|TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC", x = rownames(x = x@assays$RNA@counts), value = TRUE)
counts <- counts[-(which(rownames(counts) %in% c(black_list))),]
x = subset(x, features = rownames(counts))
pb.markers_blacklisted = FindAllMarkers(x)
annotated_pb.markers_blacklisted <- inner_join(x = pb.markers_blacklisted, 
                                   y = annotations[, c("gene_name", "description")],
                                   by = c("gene" = "gene_name")) %>%
  unique()

# Pull out the top 10 genes that define each cluster so we can visualise them

top_10.markers = annotated_pb.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top_10_unique = unique(top_10.markers$gene)

# Let's look at distribution of the clusters
# z = table(pb@active.ident, pb@meta.data$tissue)
zz = table(pb@active.ident, pb@meta.data$orig.ident)
# Define a colour scheme for the graphs
p_col = rainbow(5)
z_name = c("bm", "pb")
barplot(z, col = p_col, legend.text = TRUE)
barplot(zz, col = p_col, legend.text = TRUE)
pheatmap(z, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 20)
pheatmap(zz, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, fontsize = 20)

# Subclustering

DefaultAssay(x) = "integrated"

c0 = subset(x, subset = seurat_clusters == '0') %>% FindNeighbors(reduction = "pca", dims = 1:10) %>% FindClusters(resolution = 0.2)
c1 = subset(x, subset = seurat_clusters == '1') %>% FindNeighbors(reduction = "pca", dims = 1:10) %>% FindClusters(resolution = 0.2)
c2 = subset(x, subset = seurat_clusters == '2') %>% FindNeighbors(reduction = "pca", dims = 1:10) %>% FindClusters(resolution = 0.2)
c3 = subset(x, subset = seurat_clusters == '3') %>% FindNeighbors(reduction = "pca", dims = 1:10) %>% FindClusters(resolution = 0.2)
c4 = subset(x, subset = seurat_clusters == '4') %>% FindNeighbors(reduction = "pca", dims = 1:10) %>% FindClusters(resolution = 0.2)
c5 = subset(x, subset = seurat_clusters == '5') %>% FindNeighbors(reduction = "pca", dims = 1:10) %>% FindClusters(resolution = 0.2)
c6 = subset(x, subset = seurat_clusters == '6') %>% FindNeighbors(reduction = "pca", dims = 1:10) %>% FindClusters(resolution = 0.2)

# Is GZMx expression clonally restricted?

# Subset dominant clonotypes then query GZMx expression
y = subset(x, subset = CTstrict == 'TRAV8-6.TRAJ45.TRAC;TRAV8-6.TRAJ10.TRAC_TGTGCTGTGAGTGACCGTTCAGGAGGAGGTGCTGACGGACTCACCTTT;TGTGCTGTTCCCCCGACGGGAGGAGGAAACAAACTCACCTTT_TRBV28.TRBJ2-7.None.TRBC2_TGTGCCAGCAGTTTGGGGTTACACTACGAGCAGTACGTC')
VlnPlot(y, features = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1"), pt.size = 0, split.by = "tissue")


# Try DE testing ? Impact of environment on dominant clonotype

DefaultAssay(x) = "RNA"
y = subset(x, subset = patient == "43" & CTstrict == "TRAV26-2.TRAJ27.TRAC_TGCATCCTGAGAGACAACGGAGGCAAATCAACCTTT_TRBV19.TRBJ1-5.None.TRBC1_TGTGCCAGTAGTATATGGGGGACCTCCAATCAGCCCCAGCATTTT")
y$DE_group = NA
y$DE_group[y$tissue == 'bm'] = "bm"
y$DE_group[y$tissue == 'pb'] = "pb"
Idents(y) = y$DE_group
dominant_clone_markers = FindMarkers(y, ident.1 = "bm", ident.2 = "pb")
