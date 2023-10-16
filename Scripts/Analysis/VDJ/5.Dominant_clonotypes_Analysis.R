# Script information -------------------------------------------------------------------------------

# Title: Analysis of dominant clonotypes
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we perform DE on dominant clonotypes, specifically from PT43 & PT63
# as well as the dominant clone identified from the ERGO-II analysis.

# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '5. Analysis of dominant clonotypes' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, immunarch, scatterpie, openxlsx, EnhancedVolcano, scales)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load processed 10x data and extract the all_10x object to the global environment
all_list <- readRDS("Data/R_out/10x/all_list.rds")
all_10x <- all_list$all_10x

# Load V(D)J data in immunarch format 
# NB: These are renamed filtered_contig_annotation *.CSV files for ease of loading into R
immdata <- repLoad("Data/10x/TCR", .mode = "paired")

# Load additional resources
black_list <- readRDS("Data/R_out/black_list.rds")
annotations <- read.csv("Data/Public_data/annotations.csv") 

# Define function for annotating markers
# NB: Pass to unique at the end removes redundancy but also seems to remove some genes (?)
annotate_markers <- function(input_markers){
  input_markers <- inner_join(x = input_markers, 
                              y = annotations[, c("gene_name", "description")],
                              by = c("gene" = "gene_name")) %>% unique()
  return(input_markers)
}

# Modify black_list to retain TCR genes
black_list_TCR <- 
  black_list[!grepl("TRAV|TRAJ|TRBV|TRBJ|TRGV|TRGJ|TRDJ|TRDV|TRAC|TRBC|TRGC|TRDC", black_list)]

# Define function to remove black listed genes from discriminant genes (excluding TCR genes)
keep_TCR <- function(x){
  x <- x[!((x$gene) %in% black_list_TCR),]
  return(x)
}

# Create character vectors for naming top 10 lists
all_samples <- apply(expand.grid(c("BM", "PB"), c("13","31","43","63")), 1, paste, collapse="")
top_10_names <- expand_grid(all_samples, seq(1,10))
top_10_names$names <- str_c(top_10_names$all_samples, "_", top_10_names$`seq(1, 10)`)
top_10_names <- top_10_names$names

# Create character vector of cluster names
cluster_names <- c("T[EM]", "T[TE]","T[N]", "Cyto ~ T[EM]", "P[RE] ~ EX", "T[CM]")

# Create character vector for plot colours
plot_cols <- 
  c("T[EM]" = "#F8766D", "T[TE]" = "#B79F00", "T[N]" = "#00BA38", 
    "Cyto ~ T[EM]" = "#00BFC4", "P[RE] ~ EX" = "#619CFF", "T[CM]" = "#F564E3")

# Ensure correct slice is used downstream
slice <- dplyr::slice


# 3 Pre-processing ---------------------------------------------------------------------------------

# 3.1 Fix Seurat data ------------------------------------------------------------------------------

# 3.1.A Fix aesthetics and remove blacklisted genes for DE analysis --------------------------------

# Switch default assay to "RNA" for differential expression testing
all_list <- map(all_list, `DefaultAssay<-`, value = "RNA")

# Remove non helpful genes from DE results #NB: This will remove genes from the object
for (i in names (all_list)) {
  counts <- GetAssayData(all_list[[i]], assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c(black_list_TCR))),]
  all_list[[i]] <- subset(all_list[[i]], features = rownames(counts))
}

# 3.1.B Identify the top 10 dominant clones from the all_10x object --------------------------------

# Create an empty list to store the top 10 clones from the Seurat object
seurat_top_10_list <- list()

# Isolate the top 10 clones by Sample from the all_10x object
seurat_top_10_list$clones_43_BM <- subset(all_10x@meta.data, 
                                          subset = Patient == 43 & Tissue == "BM") %>% 
  group_by(CTaa) %>% arrange(desc(Frequency), CTaa) %>% 
  slice(1) %>% arrange(desc(Frequency)) %>% head(10)

seurat_top_10_list$clones_43_PB <- subset(all_10x@meta.data, 
                                          subset = Patient == 43 & Tissue == "PB") %>% 
  group_by(CTaa) %>% arrange(desc(Frequency), CTaa) %>% 
  slice(1) %>% arrange(desc(Frequency)) %>% head(10)

seurat_top_10_list$clones_63_BM <- subset(all_10x@meta.data, 
                                          subset = Patient == 63 & Tissue == "BM") %>% 
  group_by(CTaa) %>% arrange(desc(Frequency), CTaa) %>% 
  slice(1) %>% arrange(desc(Frequency)) %>% head(10)

seurat_top_10_list$clones_63_PB <- subset(all_10x@meta.data, 
                                          subset = Patient == 63 & Tissue == "PB") %>% 
  group_by(CTaa) %>% arrange(desc(Frequency), CTaa) %>% 
  slice(1) %>% arrange(desc(Frequency)) %>% head(10)

# Move top 10 into a single list
seurat_top_10 <- data.frame(Reduce(rbind, seurat_top_10_list))

# # 3.1.C Extract top 10 clones from the all_10x Seurat object -------------------------------------

# Create and empty list to store results
top_10_list <- list()

# Run loops to extract cells
for (i in 1:10) {
  temp <- subset(all_10x, subset = CTaa == seurat_top_10$CTaa[i] & orig.ident == "BM43")
  top_10_list[[length(top_10_list) +1]] <- temp
}

for (i in 11:20) {
  temp <- subset(all_10x, subset = CTaa == seurat_top_10$CTaa[i] & orig.ident == "PB43")
  top_10_list[[length(top_10_list) +1]] <- temp
}

for (i in 21:30) {
  temp <- subset(all_10x, subset = CTaa == seurat_top_10$CTaa[i] & orig.ident == "BM63")
  top_10_list[[length(top_10_list) +1]] <- temp
}

for (i in 31:40) {
  temp <- subset(all_10x, subset = CTaa == seurat_top_10$CTaa[i] & orig.ident == "PB63")
  top_10_list[[length(top_10_list) +1]] <- temp
}

# Fix the names of the list
names(top_10_list) <- top_10_names[c(41:80)]


# 3.2 Fix Immunarch data --------------------------------------------------------------------------- 

# Add "Sample" Column to immunarch data
immdata$data$BM13$Sample <- "BM13"
immdata$data$BM43$Sample <- "BM43"
immdata$data$BM31$Sample <- "BM31"
immdata$data$BM63$Sample <- "BM63"
immdata$data$PB13$Sample <- "PB13"
immdata$data$PB31$Sample <- "PB31"
immdata$data$PB43$Sample <- "PB43"
immdata$data$PB63$Sample <- "PB63"

# Copy the immdata$data list to the global environment to allow culling
immdata_top_10 <- immdata$data

# Reorder the list to keep consitant
immdata_top_10 <- immdata_top_10[all_samples]

# Restrict data to the top 10 clonotypes observable per sample
for (i in names(immdata_top_10)) {
  immdata_top_10[[`i`]] <- head(immdata_top_10[[`i`]], n = 10)
}

# Discard unhelpful columns
for (i in names(immdata_top_10)) {
  immdata_top_10[[`i`]] <- immdata_top_10[[`i`]][-c(8:14,18,19)]
}

# Move top 10 into a single list
immdata_top_10 <- data.frame(Reduce(rbind, immdata_top_10))


# 4. Analyse data ----------------------------------------------------------------------------------

# 4.1 Query which cluster each of the dominant clones falls within ---------------------------------

# Create an empty list to store results
cluster_list <- list()

# Run a loop to query distribution of cells
for (i in 1:length(top_10_list)){
  temp <- table(top_10_list[[`i`]][["Cluster"]])
  cluster_list[[length(cluster_list) +1]] <- temp
}

# Fix names of the list
names(cluster_list) <- top_10_names[c(41:80)]

# Convert list to a tibble for plotting
cluster_list <- as_tibble(cluster_list)

# Add rownames to tibble
rownames(cluster_list) <- cluster_names

# Transpose cluster_list
cluster_list <- t(cluster_list) %>% as.data.frame()

# Move rownames to column and fix names
cluster_list <- rownames_to_column(cluster_list, "Sample")
cluster_list <- cluster_list %>% separate(Sample, into = "Sample", sep = "_")  

# Add 1:10 into column "Rank"
cluster_list$Rank <- seq(rep(1,10))

# Add 1:40 into column "Sequential"
cluster_list$Sequential <- seq(1,40)

# Add 1:4 into column "Sample2" (geom_pie requires numerical apparently...)
cluster_list$Sample2 <- c(rep(1,10), rep(2,10), rep(3,10), rep(4,10))

# Plot data
# NB: use + coord_flip() if we want the graph rotated 90 degrees
ggplot() + geom_scatterpie(aes(x=Sample2, y=Rank, group=Sequential,r =0.25), 
                           data=cluster_list, cols=cluster_names, legend_name = "Cluster") + 
  coord_fixed() + theme_bw(base_size = 18) + xlab("Sample") + ylab("Rank") +
  scale_x_continuous(breaks = 1:4, 
                     labels = c("1" = "BM43", "2" = "PB43", "3" = "BM63", "4" = "PB63")) +
  scale_y_continuous(breaks = seq(0,10,1)) +
  scale_fill_manual(values = plot_cols, 
                    labels = label_parse()) + 
  theme(text = element_text(size = 24), axis.text = element_text(size = 24), legend.text.align = 0) + coord_flip()
ggsave(bg = "transparent", 
       filename = paste0("Output/Figures/VDJ/10x/Final/Top_10_pie_chart", ".png"), scale = 2:1)


# 4.2 Assess transcriptional differences between clones present in the BM and PB -------------------

# 4.2.A Pre-processing -----------------------------------------------------------------------------

# Identify overlapping clones in the top 10 of each sample
overlap_43 <- intersect(seurat_top_10_list$clones_43_BM$CTstrict, 
                        seurat_top_10_list$clones_43_PB$CTstrict)
overlap_63 <- intersect(seurat_top_10_list$clones_63_BM$CTstrict, 
                        seurat_top_10_list$clones_63_PB$CTstrict)

# Create empty lists to store results
top_10_DE <- list()
FindMarkers_top_10_DE_43 <- list()
FindMarkers_top_10_DE_63 <- list()
FindConservedMarkers_top_10_DE_43 <- list()
FindConservedMarkers_top_10_DE_63 <- list()
top_10_overlapping <- list()

# Run a loop to extract all clones of interest
for (i in 1:40) {
  temp <- subset(all_10x, subset = CTaa == seurat_top_10$CTaa[i])
  top_10_overlapping[[length(top_10_overlapping) +1]] <- temp
}

# Fix names
names(top_10_overlapping) = top_10_names[c(41:80)]

# Extract overlapping clones for PT43
top_10_43_overlapping = list("1" = top_10_overlapping$BM43_1, "2" = top_10_overlapping$BM43_2,
                             "3" = top_10_overlapping$BM43_5, "4" = top_10_overlapping$BM43_6,
                             "5" = top_10_overlapping$BM43_7, "6" = top_10_overlapping$BM43_8)

# Extract overlapping clones for PT63
top_10_63_overlapping = list("1" = top_10_overlapping$BM63_2, "2" = top_10_overlapping$BM63_4,
                             "3" = top_10_overlapping$BM63_7, "4" = top_10_overlapping$BM63_9)

# 4.2.B Perform FindMarkers DE for PT 43 -----------------------------------------------------------

# 4.2.A.i Run a loop to perform DE tests between the same clone in different tissues in PT43
for (i in 1:length(overlap_43)) {
  all_list$all_10x@meta.data$DE_group = NA
  all_list$all_10x@meta.data$DE_group[all_list$all_10x$CTstrict == overlap_43[i] & 
                                        all_list$all_10x$Tissue == "BM"] = "BM"
  all_list$all_10x@meta.data$DE_group[all_list$all_10x$CTstrict == overlap_43[i] & 
                                        all_list$all_10x$Tissue == "PB"] = "PB"
  Idents(all_list$all_10x) <- "DE_group"
  temp <- FindMarkers(all_list$all_10x, ident.1 = "BM", ident.2 = "PB")
  FindMarkers_top_10_DE_43[[length(FindMarkers_top_10_DE_43) +1]] <- temp %>% 
    rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("gene" = "gene_name")) %>%
    unique()
}

# Fix the names of the top_10_DE_list
names(FindMarkers_top_10_DE_43) = names(top_10_43_overlapping)

# Save out the DE results
write.xlsx(FindMarkers_top_10_DE_43, 
           "Output/DE/10x/Final/FindMarkers_top_10_DE_43.xlsx")

# 4.2.A.i # Output Volcano plots for all tests -----------------------------------------------------
for (i in 1:length(FindMarkers_top_10_DE_43)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_top_10_DE_43[[`i`]], 
                     lab = FindMarkers_top_10_DE_43[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title =names(FindMarkers_top_10_DE_43)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 16, 
                     pointSize = 2, labSize = 8, legendLabSize = 16, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x/Final/Top_10/PT_43/Different_",
                                    names(FindMarkers_top_10_DE_43)[i], ".png"))
}

# 4.2.B Perform FindMarkers DE for PT 63 -----------------------------------------------------------

# Run a loop to perform DE tests between the same clone in different tissues in PT63
for (i in 1:length(overlap_63)) {
  all_list$all_10x@meta.data$DE_group = NA
  all_list$all_10x@meta.data$DE_group[all_list$all_10x$CTstrict == overlap_63[i] & 
                                        all_list$all_10x$Tissue == "BM"] = "BM"
  all_list$all_10x@meta.data$DE_group[all_list$all_10x$CTstrict == overlap_63[i] & 
                                        all_list$all_10x$Tissue == "PB"] = "PB"
  Idents(all_list$all_10x) <- "DE_group"
  temp = FindMarkers(all_list$all_10x, ident.1 = "BM", ident.2 = "PB")
  FindMarkers_top_10_DE_63[[length(FindMarkers_top_10_DE_63) +1]] <- temp %>% 
    rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], 
               by = c("gene" = "gene_name")) %>%
    unique()
}

# Fix the names of the top_10_DE_list
names(FindMarkers_top_10_DE_63) = names(top_10_63_overlapping)

# Save out the DE results
write.xlsx(FindMarkers_top_10_DE_63, 
           "Output/DE/10x/Final/FindMarkers_top_10_DE_63.xlsx")

# 4.2.B.ii # Output Volcano plots for all tests ----------------------------------------------------
for (i in 1:length(FindMarkers_top_10_DE_63)) {
  temp <- 
    (EnhancedVolcano(FindMarkers_top_10_DE_63[[`i`]], 
                     lab = FindMarkers_top_10_DE_63[[`i`]]$gene, 
                     x = 'avg_log2FC', y = 'p_val', 
                     title = names(FindMarkers_top_10_DE_63)[i], 
                     subtitle = "Differential Expression",
                     drawConnectors = TRUE, axisLabSize = 12, 
                     pointSize = 2, labSize = 4, legendLabSize = 12, 
                     legendPosition = "bottom", max.overlaps = Inf))
  ggsave(plot = 
           temp, filename = paste0("Output/Figures/Volcanoplots/10x/Final/Top_10/PT_63/Different_",
                                    names(FindMarkers_top_10_DE_63)[i], ".png"), scale = 2:1)
}

# 4.2.C Perform FindConservedMarkers DE for PT 43 --------------------------------------------------

# Run for loop to iterate across all samples
# NB: There is insufficient data to compare between clusters for each clone, instead we will
# compare global similarities and if required, subset specific clusters later
top_10_43_overlapping <- map(top_10_43_overlapping, `Idents<-`, value = "Tissue")
for (i in names(top_10_43_overlapping)) {
  temp <- FindConservedMarkers(top_10_43_overlapping[[i]], ident.1 = "BM", ident.2 = "PB", 
                               grouping.var = "Patient")
  FindConservedMarkers_top_10_DE_43[[length(FindConservedMarkers_top_10_DE_43) +1]] <- temp %>% 
    rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>% 
    unique()
}

# Fix the names of the top_10_DE_list
names(FindConservedMarkers_top_10_DE_43) = names(top_10_43_overlapping)

# Save out the DE results
write.xlsx(FindConservedMarkers_top_10_DE_43, 
           "Output/DE/10x/Final/FindConservedMarkers_top_10_DE_43.xlsx")

# 4.2D Perform FindConservedMarkers DE for PT 63 ---------------------------------------------------

# Run for loop to iterate across all samples
# NB: There is insufficient data to compare between clusters for each clone
top_10_63_overlapping <- map(top_10_63_overlapping, `Idents<-`, value = "Tissue")
for (i in names(top_10_63_overlapping)) {
  temp <- FindConservedMarkers(top_10_63_overlapping[[i]], ident.1 = "BM", ident.2 = "PB", 
                               grouping.var = "Patient")
  FindConservedMarkers_top_10_DE_63[[length(FindConservedMarkers_top_10_DE_63) +1]] <- temp %>% 
    rownames_to_column %>% rename(gene = rowname) %>% 
    inner_join(y = annotations[, c("gene_name", "description")], by = c("gene" = "gene_name")) %>% 
    unique()
}

# Fix the names of the top_10_DE_list
names(FindConservedMarkers_top_10_DE_63) = names(top_10_63_overlapping)

# Save out the DE results
write.xlsx(FindConservedMarkers_top_10_DE_63, 
           "Output/DE/10x/Final/FindConservedMarkers_top_10_DE_63.xlsx")


# 5. Save out data ---------------------------------------------------------------------------------

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean up the environment
rm(all_10x, all_list, counts, temp, 
   temp_list, top_10_markers_all_10x, top_10_markers_BM_10x, top_10_markers_PB_10x, black_list, 
   ccolours_5, ccolours_6, ccolours_7, cluster_names_5, cluster_names_6, cluster_names_7, 
   filename, i, j, top_10_unique_all_10x, top_10_unique_BM_10x, top_10_unique_PB_10x)

# Save out the DE_results for further analysis
save.image("Data/R_out/Controls_1/DE_Workup_Controls_1.rds")

gc()

print("#># Finished running '5. Dominant clonotype analysis")