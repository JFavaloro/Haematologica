# Script information -------------------------------------------------------------------------------

# Title: ERGO-II analysis
# Author: James Favaloro and Samuel Gardiner
# Date: 2023-08-21
# Description: In this script we will analyse our ERGO-II hits against our clonotype data. If 
# expansions within our clonotype data are as a result of an expansion of anti-MM clones, we may see 
# disproportionate representation within the BM. We will first examine the overall results but 
# restrict analysis to clones >= Medium expansion (>=0.1%)


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running 'ERGO-II analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, immunarch, VennDiagram, scales, ggalluvial)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load previously saved workspace
load("Data/R_out/ERGO/ERGO_results.RData")

# Load processed 10x data and bring the all_10x object into the global environment
all_list <- readRDS("Data/R_out/10x/all_list.rds")
all_10x <- all_list$all_10x

# Define plotting vectors
plot_colours <- c("blue", "red")
tissue <- c(rep("BM", 4), rep("PB", 4))
cluster_names <- c(expression("T[EM]", "T[TE]", "T[N]", "Cyto ~ T[EM]", "P[RE] ~ EX", "T[CM]"))
cluster_names <- c("T[EM]", "T[TE]", "T[N]", "Cyto ~ T[EM]", "P[RE] ~ EX", "T[CM]")


# 3. Pre-processing --------------------------------------------------------------------------------

# Rename metadata to allow prettier plotting
all_10x@meta.data <- rename(all_10x@meta.data, Sample = orig.ident)


# 4. Analyse data ----------------------------------------------------------------------------------

# 4.1 Global analysis of data - how many hits do we see? -------------------------------------------

# 4.1.1A Return the number of clones we observe against each of the 197 peptides -------------------
peptide_results <- table(ERGO_results_long$Peptide) %>% as.data.frame()
# We see hits against 145 of the 197 peptides

# Restrict further analysis to Medium expansion or greater (>=0.1%)
peptide_results <- subset(ERGO_results_long, subset = Proportion >= 0.001)
peptide_results <- table(peptide_results$Peptide) %>% as.data.frame()
colnames(peptide_results) <- c("Peptide", "Frequency")
peptide_results <- arrange(peptide_results, desc(Frequency))
# We see hits against 29 of the 197 peptides

# Restrict further analysis to highest scoring Peptide result per input clone
peptide_results <- subset(ERGO_results_long, subset = Proportion >= 0.001)
peptide_results <- peptide_results %>% arrange(desc(Score)) %>% group_by(Barcode) %>% slice(1)
peptide_results <- table(peptide_results$Peptide) %>% as.data.frame()
colnames(peptide_results) <- c("Peptide", "Frequency")
peptide_results <- arrange(peptide_results, desc(Frequency))
# We see hits against 6 of 197 peptides

# 4.1.1B Return the number of clones we observe against each of the 58 proteins --------------------
protein_results <- table(ERGO_results_long$Protein) %>% as.data.frame()
# We see hits against 54 of the 58 proteins

# Restrict further analysis to Medium expansion or greater (>=0.1%)
protein_results <- subset(ERGO_results_long, subset = Proportion >= 0.001)
protein_results <- table(protein_results$Protein) %>% as.data.frame()
colnames(protein_results) <- c("Protein", "Frequency")
protein_results <- arrange(protein_results, desc(Frequency))
# We see hits against 24 of the 58 proteins

# Restrict further analysis to highest scoring Protein result per input clone
protein_results <- subset(ERGO_results_long, subset = Proportion >= 0.001)
protein_results <- protein_results %>% arrange(desc(Score)) %>% group_by(Barcode) %>% slice(1)
protein_results <- table(protein_results$Protein) %>% as.data.frame()
colnames(protein_results) <- c("Protein", "Frequency")
protein_results <- arrange(protein_results, desc(Frequency))
# We see hits against 6 of the 58 proteins

# 4.1.1C Tabulate the number of hits per sample based on clone bin sizes ---------------------------

# Tabulate all results
table(ERGO_results_wide$Sample)

# Tabulate the number of hits against hyperexpanded clones (i.e. >10% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.1, ERGO_results_wide$Sample)

# Tabulate the number of hits against large clones (i.e. between 1-10% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.01 & ERGO_results_wide$Proportion < 0.1, 
      ERGO_results_wide$Sample)

# Tabulate the number of hits against medium clones (i.e. between 0.1-1% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.001 & ERGO_results_wide$Proportion < 0.01, 
      ERGO_results_wide$Sample)

# Tabulate the number of hits against small clones (i.e. between 0.01-0.01% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.0001 & ERGO_results_wide$Proportion < 0.001, 
      ERGO_results_wide$Sample)

# Tabulate the number of hits against single clones (i.e. single reads)
table(ERGO_results_wide$Clones == 1, ERGO_results_wide$Sample)


# 4.2 Tissue specific analysis of data - do we see disproportionate representation in one tissue? --

# 4.2.1A Peptide analysis - Barplot ----------------------------------------------------------------

# Prep the data into a format suitable for plotting
peptide_results_tissue <- subset(ERGO_results_long, subset = Proportion >= 0.001)
peptide_results_tissue <- peptide_results_tissue %>% arrange(desc(Score)) %>% 
  group_by(Barcode) %>% slice(1)
peptide_results_tissue <- table(peptide_results_tissue$Peptide, peptide_results_tissue$Tissue) %>% 
  as.data.frame()
colnames(peptide_results_tissue) <- c("Peptide", "Tissue", "Frequency")
peptide_results_tissue <- arrange(peptide_results_tissue, desc(Frequency))

# Graph the results
ggplot(peptide_results_tissue, aes(x = reorder(Peptide, Frequency), y = Frequency, fill = Tissue)) + 
  geom_bar(stat = "identity") +
  xlab("Peptide") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(breaks = seq(0,8,1)) +
  scale_fill_manual(values = plot_colours) +
  labs(title = "Number of highest scoring hits against peptides in medium or greater expansions") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(text = element_text(size = 16))
ggsave("Output/Figures/ERGO/peptide_results_tissue.png", width = 12, height = 15)

# 4.2.1B Peptide analysis - Venn Diagram -----------------------------------------------------------

# Prep the data into a format suitable for plotting
results <- subset(ERGO_results_long, subset = Proportion >= 0.001)
results <- results %>% arrange(desc(Score)) %>% 
  group_by(Barcode) %>% slice(1)
results = table(results$Peptide, results$Tissue) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results, category.names = c("BM", "PB"), 
             imagetype = "svg", height = 10, width = 10, resolution = 500,
             filename = "Output/Figures/ERGO/peptide_venn_diagram.svg",
             output = TRUE, disable.logging = TRUE,
             col=c("blue", 'red'),
             fill = c('blue', "red"), print.mode = c('percent', 'raw'))


# 4.2.2A Protein analysis - Barplot ----------------------------------------------------------------

# Prep the data into a format suitable for plotting
protein_results_tissue <- subset(ERGO_results_long, subset = Proportion >= 0.001)
protein_results_tissue <- protein_results_tissue %>% arrange(desc(Score)) %>% 
  group_by(Barcode) %>% slice(1)
protein_results_tissue = table(protein_results_tissue$Protein, protein_results_tissue$Tissue) %>%
  as.data.frame()
colnames(protein_results_tissue) = c("Protein", "Tissue", "Frequency")
protein_results_tissue = arrange(protein_results_tissue, desc(Frequency))

# Graph the results
ggplot(protein_results_tissue, aes(x = reorder(Protein, Frequency), y = Frequency, fill = Tissue)) + 
  geom_bar(stat = "identity") +
  xlab("Protein") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(breaks = seq(0,8,1)) +
  scale_fill_manual(values = plot_colours) +
  labs(title = "Number of highest scoring hits against proteins in medium or greater expansions") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("Output/Figures/ERGO/protein_results_tissue.png", width = 10, height = 5)


# 4.2.2B Protein analysis - Venn Diagram -----------------------------------------------------------

# Prep the data into a format suitable for plotting
results <- subset(ERGO_results_long, subset = Proportion >= 0.001)
results <- results %>% arrange(desc(Score)) %>% 
  group_by(Barcode) %>% slice(1)
results = table(results$Protein, results$Tissue) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results, category.names = c("BM", "PB"), 
             imagetype = "svg", height = 10, width = 10, resolution = 500,
             filename = "Output/Figures/ERGO/protein_venn_diagram.svg",
             output = TRUE, disable.logging = TRUE,
             col=c("blue", 'red'),
             fill = c('blue', "red"), print.mode = c('percent', 'raw'))


# 4.3 Transcriptional analysis of identified clones ------------------------------------------------
# NB: Analysis is restricted to PT43 and PT63 as identified clones are predominantly found in the BM

# Isolate clones from PT 43 and PT63 along with their binding information
seurat_clones <- subset(ERGO_results_long, subset = Proportion >= 0.001)
seurat_clones <- filter(seurat_clones, Sample %in% c("BM43", "PB43", "BM63", "PB63"))
seurat_clones <- seurat_clones %>% arrange(desc(Score)) %>% group_by(Barcode) %>% slice(1)
seurat_clones <- seurat_clones[c(2,5,6)]

# Define function to plot all 7 clonotypes to which we have paired transcriptomic data

# Have run trace(alluvialClonotypes, edit=TRUE) and removed lines 50:51 to remove original geom_text:
#plot <- plot + geom_text(stat = ggalluvial::StatStratum, 
#infer.label = FALSE, reverse = TRUE, size = 2)

alluvial_subset <- function(TRB, TRBV, TRBJ) {
  y <- all_10x@meta.data %>% filter(str_detect(all_10x@meta.data$CTaa, TRB)) 
  # Get all clonotypes with a specific CDR3B
  y <- y %>% filter(str_detect(y$CTgene, TRBV)) # self explanatory
  y <- y %>% filter(str_detect(y$CTgene, TRBJ)) # self explanatory
  y <- subset(all_10x, subset = CTaa %in% y$CTaa)
  
  # Use scRepertoire to show cluster distribution of dominant clonotype
  alluvialClonotypes(y, cloneCall = "strict", 
                     y.axes = c("Sample", "Cluster"), 
                     color = "Cluster") +
    guides(fill = guide_legend("Cluster")) +
    labs(y = "Number of cells", x = NA,
         title = "Transcriptional state of identified clone") + theme_bw() +
    theme(axis.title.x=element_blank()) +
    scale_fill_manual(values = c("T[EM]" = "#F8766D", "T[TE]" = "#C49A00", 
                                 "T[N]" = "#53B400", "Cyto ~ T[EM]" = "#00C094", 
                                 "P[RE] ~ EX" = "#00B6EB", "T[CM]" = "#A58AFF"), 
                      labels = label_parse()) + 
    theme(legend.text.align = 0, text = element_text(size = 16)) +
    scale_x_discrete(expand = c(0.1,0.1)) +
    geom_text(stat = "stratum", size = 12, parse = T) + theme(axis.text = element_text(size = 12))
    paste("Output/Figures/ERGO/", TRBV, "_", TRB, "_", TRBJ, ".png") %>% 
      ggsave(scale = 1:2)
}

# Run function
pwalk(seurat_clones, alluvial_subset)


# 7. Save out data and print stats -----------------------------------------------------------------

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '4. ERGO-II analysis' script")