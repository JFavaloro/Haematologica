# Script information -------------------------------------------------------------------------------

# title: TRM Signature
# author: James Favaloro
# date: 28/06/2021
# description: In this script we will utilise the AddModuleScore feature of Seurat to determine if BM-resident CD8+CD69+T-cells hold similar gene expression signatures to those as determined by Kumar et al. (2017) (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5646692/). The gene expression profiles (NIHMS903158-supplement-2.xlsx) for CD8+CD69+T-cells extracted from both Lung and Spleen (and compared to CD8+CD69-T-cells) have been downloaded from (DOI: 10.1016/j.celrep.2017.08.078) and saved to the directory "Data/Public". NB: As 10x data cannot reliably be segregated based on the expression of a gene, this analysis will focus on the differences between BM and PB where we have previously demonstrated a near complete absence of CD69 expression on CD8+T-cells from the PB. This is however different when assessing gene expression, as cells from PB do express some level of CD69 mRNA...


# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'TRM signature' script")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(Seurat)
  library(scRepertoire)
  library(immunarch)
  library(escape)
})


# 2. Load data and define functions ----------------------------------------------------------------

# 2.1 Load our data and set the combined Seurat object to the appropriate resolution
integrated.list = readRDS("Data/integrated.list.rds")
x = integrated.list$all_ds.integrated
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"

# 2.2 Load data from Kumar et al. (2017)
CD8_lung = read_xlsx("Data/Public/NIHMS903158-supplement-2.xlsx", skip = 1, range = "G3:K2003")
CD8_spleen = read_xlsx("Data/Public/NIHMS903158-supplement-2.xlsx", skip = 1, range = "S3:W2003") 


# Something

# Extract the top 200 up and down regulated genes based on LogFC and convert to character string for use in heatmap
CD8_lung_up = CD8_lung[order(CD8_lung$logFC, decreasing = TRUE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames()
CD8_lung_down = CD8_lung[order(CD8_lung$logFC, decreasing = FALSE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames()

CD8_spleen_up = CD8_spleen[order(CD8_lung$logFC, decreasing = TRUE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames()
CD8_spleen_down = CD8_spleen[order(CD8_lung$logFC, decreasing = FALSE), ] %>% column_to_rownames("Gene") %>% head(n = 200) %>% rownames()

# Set appropriate parameters for our Seurat object
DefaultAssay(x) = "RNA"
Idents(x) = "tissue"

# Let's check to see if we've correctly selected our a & b
VlnPlot(x, features = "CD69")

# Call a heatmap to look at the proportion of bm vs. pb cells that express TRM signature genes
# NB: Gives a number of errors due to the low expression of genes in our dataset - not sure what the eamhb... is in the heatmap...
# Lung: 32/200
# Spleen: 63/200
DoHeatmap(x, features = CD8_lung_up, assay = "SCT")
DoHeatmap(x, features = CD8_spleen_up, assay = "SCT")

# Use Seurat's AddModuleScore to highlight cells from our dataset that share similarities to other TRM. NB:this function needs genes of interest in a list...

# Move our genes of interest into lists
CD8_lung_up = list(CD8_lung_up)
CD8_spleen_up = list(CD8_spleen_up)

# Append this information to our Seurat object
x = AddModuleScore(x, features = CD8_lung_up, name = "CD8_lung_up")
x = AddModuleScore(x, features = CD8_spleen_up, name = "CD8_spleen_up")

# Visualise the data
FeaturePlot(x, features = "CD8_lung_up1", split.by = "tissue")
FeaturePlot(x, features = "CD8_spleen_up1", split.by = "tissue")

# Visualise the data using merged image with CD69
FeaturePlot(x, features = c("CD8_lung_up1", "CD69"), split.by = "tissue", blend = T)
FeaturePlot(x, features = c("CD8_spleen_up1", "CD69"), split.by = "tissue", blend = T)

# Visualise the data using merged image with GZMB
FeaturePlot(x, features = c("CD8_lung_up1", "GZMB"), split.by = "tissue", blend = T)
FeaturePlot(x, features = c("CD8_spleen_up1", "GZMB"), split.by = "tissue", blend = T)



Idents(x) = "seurat_clusters"
DimPlot(x)
FeaturePlot(x, features = "GZMB")

# Trial escape

GS <- getGeneSets(library = "H")

# Prep Seurat object
x = integrated.list$all_ds.integrated
x$seurat_clusters = x$integrated_snn_res.0.2
Idents(x) = "seurat_clusters"
DimPlot(x)
DefaultAssay(x) = "RNA"

# Run scGSEA across all cells in our combined seurat object
ES <- enrichIt(x, gene.sets = GS, groups = 1000, cores = 8)

# Append the results back to the Seurat object as metadata
x <- AddMetaData(x, ES)

# Save out metadata
yy = x@meta.data

# Add TRM signature to ES
ES = cbind(ES, yy$CD8_lung_up1)
ES = cbind(ES, yy$CD8_spleen_up1)

# Plot results as heatmap
dittoHeatmap(x, genes = NULL, metas = names(ES), 
             annot.by = "tissue", 
             fontsize = 7, 
             cluster_cols = TRUE,
             heatmap.colors = rev(colors(50)))

# This is a little hard to comprehend - try splitting into clusters
y = subset(x, subset = seurat_clusters == '0')
dittoHeatmap(y, genes = NULL, metas = names(ES), 
             annot.by = "tissue", 
             fontsize = 7, 
             cluster_cols = TRUE,
             heatmap.colors = rev(colors(50)))

# Still hard to comprehend! Convert Seurat object to SingleCellExperiment
# and plot as Violin plot - selecting gene sets that may show more helpful

x.sce = as.SingleCellExperiment(x)

# I think we may need to repeat the GSEA to get this into the right format
#ES.sce <- enrichIt(obj = x.sce, gene.sets = GS, groups = 1000, cores = 8)
# For future reference, no - output is the same; easier to simply carry over the
# ES results with: 
ES.sce = ES
# Using the above seurat object where we've copied over the ES - looks as if
# it's copied across everything - why then is it not plotting the scatterHex?
# Looks like it's an issue of not correctly directing where metadata is stored

met.data <- merge(colData(x.sce), ES.sce, by = "row.names", all=TRUE)
row.names(met.data) <- met.data$Row.names
# met.data$Row.names <- NULL
colData(x.sce) <- met.data

multi_dittoPlot(x.sce, vars = c("HALLMARK_GLYCOLYSIS", "HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY"), 
                group.by = "tissue", split.by = "seurat_clusters", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

xx = multi_dittoPlot(x.sce, vars = c("CD8_lung_up1", "CD8_spleen_up1"), group.by = "tissue", split.by = "seurat_clusters", list.out = TRUE)

multi_dittoPlot(x.sce, vars = c("CD8_lung_up1", "CD8_spleen_up1"), 
                group.by = "tissue", split.by = "seurat_clusters", plot = "vlnplot", 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

dittoScatterHex(x.sce, x.var = x.sce@colData@listData[["HALLMARK_APOPTOSIS"]], 
                y.var = x.sce@colData@listData[["HALLMARK_HYPOXIA"]], 
                do.contour = TRUE,
                split.by = "seurat_clusters")  + 
  scale_fill_gradientn(colors = rev(colors(11))) 

# Do this for clusters 0-6
y = subset(x, subset = seurat_clusters == 0)
ES_cluster0 <- enrichIt(y, gene.sets = GS, groups = 1000, cores = 8)
# Let's try extract data split by tissue in a little bit but first see if we can
# visualise the output first
x = AddMetaData(x, ES_cluster0)
dittoHeatmap(x, assay = "RNA", genes = NULL, metas = names(ES), 
             annot.by = "tissue", 
             fontsize = 7, 
             cluster_cols = TRUE,
             heatmap.colors = rev(colors(50)))



### WIP ###

# Do something a bit fancier

PrctCellExpringGene(x, genes = CD8_lung_up, group.by = "tissue")

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

PrctCellExpringGene(x, genes = genes_to_plot, group.by = "tissue")


# Add column to metadata which defines our two populations of interest
# NB: Probably not advisable to do this - just because cell expressess low level
# of a gene doesn't mean it's not expressing it - likely technical reasons
# After discussion with Walter we won't do this. Instead, we'll just use bm & pb
# BUT - here's code that allows selection of cell based on both metadata and
# feature level data
#x$DE_group = NA
#x$DE_group[x$tissue == 'bm' & GetAssayData(x, slot = "data")["CD69",] > 1] = "a"
#x$DE_group[x$tissue == 'bm' & GetAssayData(x, slot = "data")["CD69",] > 2] = "b"


# Let look at the numbers of cells that express CD69 and the level at which they express it
# Thank you Walter. Gotta figure out why I can't do this all within R.
y = as.matrix(GetAssayData(x, assay = "RNA", slot = "data")["CD69",])
write.csv(y, "y.csv")
expression_stats = read_csv("y.csv", col_names = c("cb", "CD69"), skip = 1) %>% 
mutate(sample_name = as.factor(str_split_fixed(cb, "_", 2) [,1])) %>%
mutate(sample_type = substr(sample_name, start = 1, stop = 2)) %>% 
mutate(rank = row_number())
quantile(expression_stats$CD69[expression_stats$sample_type=="bm"])
quantile(expression_stats$CD69[expression_stats$sample_type=="pb"])
median(expression_stats$CD69[expression_stats$sample_type=="bm"])
median(expression_stats$CD69[expression_stats$sample_type=="pb"])  

# Set identity to allow DE and visualisation
x = SetIdent(x, value = "DE_group")

# Perform DE
DE = FindMarkers(x, ident.1 = "a", ident.2 = "b")

# OLD METHOD
# Select cells we're interested in performing DE on

a = subset(all.integrated, subset = CTaa == 'CILRDNGGKSTF_CASSIWGTSNQPQHF' & patient == '43' & tissue == 'pb' & seurat_clusters == "1")
b = subset(all.integrated, subset = seurat_clusters == "2")
b = subset(b, subset = CTaa == 'CILRDNGGKSTF_CASSIWGTSNQPQHF', invert = TRUE)
b = subset(b, subset = tissue == 'pb')

# Extract the cellbarcodes for those cells

cellnames.a = rownames(a@meta.data)
cellnames.b = rownames(b@meta.data)

# Perform DE

DE = FindMarkers(x, ident.1 = cellnames.a, ident.2 = cellnames.b)

# Show where the cells cluster

DimPlot(x, cells.highlight = c(cellnames.a, cellnames.b))


# Add P values to comparisons https://www.biostars.org/p/458261/

# Define genes of interest
gene_sig = c("CD69")

# Define function for plotting and 
vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
     plot_case1 <- function(signature){
         VlnPlot(x, features = signature,
                                     pt.size = 0.1, 
                                       group.by = "tissue", 
                                       y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
                               ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
         }
       purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
       file_name <- paste0(file_name, "_r.png")
       ggsave(file_name, width = 14, height = 8)
}

# Run
vp_case1(gene_signature = gene_sig, file_name = "x", test_sign = comparisons, y_max = 7)


# Old code

# Let's try seeing what spread of CD69 is across BM samples
a = subset(x, CD69 >1 & tissue == 'bm')

# Set the identity to tissue to look at sample as a whole & work out average CD69 expression in the BM
Idents(a) = "tissue"
AverageExpression(a, features = "CD69", slot = "data")
b = subset(a, CD69 >2)


x$DE_group[x$CTaa == 'CILRDNGGKSTF_CASSIWGTSNQPQHF' & x$patient == '43' & x$tissue == 'pb' & x$seurat_clusters != "1"] = "b"
