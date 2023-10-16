# Script information -------------------------------------------------------------------------------

# Title: ReactomeGSA - Analysis
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will use ReactomeGSA to determine enrichment across a number of 
# cellular processes using the vignette:
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html
# in conjunction with reactome.org. This script is reasonably modular - the pathways_to_plot list 
# can be amended if required.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '5. ReactomeGSA - Analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, ReactomeGSA, fs, scales, pheatmap)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load GSVA results
all_10x_gsva_results <- readRDS("Data/R_out/GSEA/all_10x_gsva_results.rds")
merged_gsva_results <- readRDS("Data/R_out/GSEA/merged_gsva_results.rds")

# Load additional resources
pathways_to_plot <- read.csv("Data/Public_data/Pathways_for_plotting.csv", 
                             na.strings = "", sep = ",") %>% lapply(na.omit)
pathway_names <- read.csv("Data/Public_data/Pathways_for_plotting.csv", 
                          na.strings = "", sep = ",") %>% lapply(na.omit) %>% unlist() 

# Define aesthetics for plotting
all_10x_plot_cols <- 
  c("T[EM]~BM" = "#F8766D", "T[EM]~PB" = "#F8766D", "T[TE]~BM" = "#D39200", "T[TE]~PB" = "#D39200",
    "T[N]~BM" = "#93AA00", "T[N]~PB" = "#93AA00", "Cyto~T[EM]~BM" = "#00BA38", "Cyto~T[EM]~PB" = "#00BA38", 
    "P[RE]~EX~BM" = "#00C19F", "P[RE]~EX~PB" = "#00C19F", "T[CM]~BM" = "#619CFF", "T[CM]~PB" = "#619CFF")

# Not sure this will work
merged_plot_cols <- 
  c("T[N]~control" = "#F8766D", "T[N]~MM" = "#F8766D", 
    "T[CM]~control" = "#D39200", "T[CM]~MM" = "#D39200",
    "T[EM]~control" = "#93AA00", "T[EM]~MM" = "#93AA00", 
    "T[TE]~control" = "#619CFF", "T[TE]~MM" = "#619CFF")

all_10x_plot_names <- 
  c("Cyto~T[EM]~BM", "Cyto~T[EM]~PB", "P[RE]~EX~BM", "P[RE]~EX~PB", 
    "T[CM]~BM", "T[CM]~PB", "T[EM]~BM", "T[EM]~PB", "T[N]~BM", "T[N]~PB", "T[TE]~BM", "T[TE]~PB")

merged_plot_names_BM <- 
  c("T[CM]~control", "T[CM]~MM", "T[EM]~control", "T[EM]~MM", "T[N]~MM",
    "T[TE]~control",  "T[TE]~MM")

merged_plot_names_PB <- 
  c("T[CM]~control", "T[CM]~MM", "T[EM]~control", "T[EM]~MM", 
    "T[N]~control", "T[N]~MM", "T[TE]~control",  "T[TE]~MM")

# Rename columns for all_10x_gsva_results
colnames(all_10x_gsva_results@results$Seurat$pathways) <- 
  c("Pathway", "Name", all_10x_plot_names)
colnames(all_10x_gsva_results@results$Seurat$fold_changes) <- 
  c("Identifier", all_10x_plot_names)

# Rename columns for the merged gsva results
for (i in names(merged_gsva_results)[1:2]) {
  colnames(merged_gsva_results[[`i`]]@results$Seurat$pathways) <- 
    c("Pathway", "Name", merged_plot_names_BM)
}
for (i in names(merged_gsva_results)[1:2]) {
  colnames(merged_gsva_results[[`i`]]@results$Seurat$fold_changes) <- 
    c("Identifier", merged_plot_names_BM)
}

for (i in names(merged_gsva_results)[3:6]) {
  colnames(merged_gsva_results[[`i`]]@results$Seurat$pathways) <- 
    c("Pathway", "Name", merged_plot_names_PB)
}
for (i in names(merged_gsva_results)[3:6]) {
  colnames(merged_gsva_results[[`i`]]@results$Seurat$fold_changes) <- 
    c("Identifier", merged_plot_names_PB)
}

# Create directories required for output of figures
dir.create("temp/Activated BM/Top_10", recursive = TRUE)
dir.create("temp/Resting BM/Top_10", recursive = TRUE)
dir.create("temp/Activated PB/Top_10", recursive = TRUE)
dir.create("temp/Resting PB/Top_10", recursive = TRUE)
dir.create("temp/Age-matched PB/Top_10", recursive = TRUE)
dir.create("temp/Young PB/Top_10", recursive = TRUE)
dir.create("Output/Figures/Barplots/10x/Final/GSVA/Top_10/", recursive = TRUE)
dir.create("Output/Figures/Barplots/Controls_1/Final/GSVA/Top_10", recursive = TRUE)
dir.create("Output/Figures/Barplots/Controls_3/Final/GSVA/Top_10", recursive = TRUE)


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 all_10x pre-processing -----------------------------------------------------------------------

# Obtain pathway level expression and fix names
all_10x_pathway_expression <- pathways(all_10x_gsva_results)
colnames(all_10x_pathway_expression) <- gsub("\\.Seurat", "", colnames(all_10x_pathway_expression))


# 3.2 Controls pre-processing ----------------------------------------------------------------------

# Create an empty list to store results
merged_pathway_expression <- list()

# Run loop to obtain pathway level expression and fix names
for (i in names(merged_gsva_results)) {
  temp <- pathways(merged_gsva_results[[`i`]])
  merged_pathway_expression[[length(merged_pathway_expression) +1]] <- temp
}

# Fix the names in the merged_pathway_expression object
names(merged_pathway_expression) <- names(merged_gsva_results)
for (i in names(merged_pathway_expression)) {
  colnames(merged_pathway_expression[[`i`]]) <- 
    gsub("\\.Seurat", "", colnames(merged_pathway_expression[[`i`]]))
}


# 4. Analysis  -------------------------------------------------------------------------------------

# 4.1 Find pathways with greatest differences between clusters for the all_10x object --------------
all_10x_max_difference <- do.call(rbind, apply(all_10x_pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

# Reorder and sort results based on the difference
all_10x_max_difference$diff <- all_10x_max_difference$max - all_10x_max_difference$min
all_10x_max_difference <- 
  all_10x_max_difference[order(all_10x_max_difference$diff, decreasing = T), ]


# 4.2 Find pathways with greatest differences between clusters for Controls ------------------------

# Create an empty list to store results
merged_max_difference <- list()

# Run a loop to determine max differences
for (i in names(merged_pathway_expression)) {
  temp <- do.call(rbind, apply(merged_pathway_expression[[`i`]], 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
  }))
  merged_max_difference[[length(merged_max_difference) +1]] <- temp
}

# Rename objects in the merged_max_difference list
names(merged_max_difference) <- names(merged_gsva_results)

# Run a loop to reorder and sort results based on the difference
for (i in names(merged_max_difference)) {
  merged_max_difference[[`i`]][["diff"]] <- 
    merged_max_difference[[`i`]][["max"]] - merged_max_difference[[`i`]][["min"]]
  merged_max_difference[[`i`]] <- 
    merged_max_difference[[`i`]][order(merged_max_difference[[`i`]][["diff"]], decreasing = T), ]
}


# 4.3 Plot pathways in the all_10x object ----------------------------------------------------------

# 4.3.1 Plot the top 10 most differential pathways
for (i in 1:10) {
  temp <- plot_gsva_pathway(all_10x_gsva_results, pathway_id = rownames(all_10x_max_difference)[i]) + 
    geom_bar(aes(fill = cluster_id), stat = "identity") +
    scale_fill_manual(values = all_10x_plot_cols,
                      labels = label_parse()) +
    scale_x_discrete(labels = label_parse()) +
    guides(fill = guide_legend("Cluster")) +
    theme(legend.text.align = 0)
  ggsave(plot = temp, filename = paste0("Output/Figures/Barplots/10x/Final/GSVA/Top_10/", 
                                        path_sanitize(all_10x_max_difference$name[[`i`]]), ".png"))
}


# 4.3.2 Plot pathways of interest
for (i in names(pathway_names)) {
  temp <- plot_gsva_pathway(all_10x_gsva_results, pathway_id = pathway_names[i]) + 
    geom_bar(aes(fill = cluster_id), stat = "identity") +
    scale_fill_manual(values = all_10x_plot_cols, 
                      labels = label_parse()) + 
    scale_x_discrete(labels = label_parse()) +
    guides(fill = guide_legend("Cluster")) +
    theme(legend.text.align = 0)
  ggsave(plot = temp, filename = paste0("Output/Figures/Barplots/10x/Final/GSVA/", 
                                        pathway_names[i], ".png"))
}


# 4.4 Plot pathways in the for Controls  -----------------------------------------------------------

# 4.4.1 Plot the top 10 most differential pathways
for (i in names(merged_gsva_results)) {
  for (j in 1:10) {
    temp <- plot_gsva_pathway(merged_gsva_results[[`i`]], 
                              pathway_id = rownames(merged_max_difference[[`i`]])[j]) + 
      geom_bar(aes(fill = cluster_id), stat = "identity") +
      scale_fill_manual(values = merged_plot_cols, 
                        labels = label_parse()) + 
      scale_x_discrete(labels = label_parse()) +
      guides(fill = guide_legend("Cluster")) +
      theme(legend.text.align = 0)
    ggsave(plot = temp, filename = 
             paste0("temp/", names(merged_gsva_results[`i`]), "/Top_10/", 
                    path_sanitize(merged_max_difference[[`i`]][["name"]][[`j`]]), ".png"))
  }
}


# 4.4.2 Plot pathways of interest
for (i in names(merged_gsva_results)) {
  for (j in names(pathway_names)) {
    temp <- plot_gsva_pathway(merged_gsva_results[[`i`]], pathway_id = pathway_names[j]) + 
      geom_bar(aes(fill = cluster_id), stat = "identity") +
      scale_fill_manual(values = merged_plot_cols, 
                        labels = label_parse()) + 
      scale_x_discrete(labels = label_parse()) + 
      guides(fill = guide_legend("Cluster")) +
      theme(legend.text.align = 0)
    ggsave(plot = temp, filename = 
             paste0("temp/", names(merged_gsva_results[`i`]), "/",
                    pathway_names[j], ".png"))
  }
}


# 5. Save out data and print stats -----------------------------------------------------------------

# Move folders to appropriate locations
file.rename("temp/Activated BM/", 
            "Output/Figures/Barplots/Controls_1/Final/GSVA/Activated BM/")
file.rename("temp/Resting BM/", 
            "Output/Figures/Barplots/Controls_1/Final/GSVA/Resting BM/")
file.rename("temp/Activated PB/", 
            "Output/Figures/Barplots/Controls_1/Final/GSVA/Activated PB/")
file.rename("temp/Resting PB/", 
            "Output/Figures/Barplots/Controls_1/Final/GSVA/Resting PB/")
file.rename("temp/Age-matched PB/", 
            "Output/Figures/Barplots/Controls_3/Final/GSVA/Age-matched PB/")
file.rename("temp/Young PB/", 
            "Output/Figures/Barplots/Controls_3/Final/GSVA/Young PB/")

# Clean the temp directory
unlink("temp", recursive = TRUE)

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '5. ReactomeGSA - Analysis' script")