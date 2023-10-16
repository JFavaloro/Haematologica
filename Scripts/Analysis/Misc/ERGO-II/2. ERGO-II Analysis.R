# Script information -------------------------------------------------------------------------------

# title: ERGO-II analysis
# author: James Favaloro and Samuel Gardiner
# date: 28/06/2021
# description: In this script we will analyse our ERGO-II hits against our clonotype data. If expansions within our clonotype data are as a result of an expansion of anti-MM clones, we may see disproportionate representation within the BM.


# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'ERGO-II analysis' script")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scRepertoire)
  library(immunarch)
  library(VennDiagram)
  library(UpSetR)
  library(ggalluvial)
})


# 2. Load data and define functions ----------------------------------------------------------------

# 2.1 Load previously saved workspace
load("Output/ERGO/ERGO_results.RData")

# 2.2 Define functions
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))

# 3. Adjust environment for plotting data ----------------------------------------------------------

plot_colours = c("red", "blue")
tissue = c("BM", "BM", "BM", "BM", "PB", "PB", "PB", "PB")

# Show colours we're using to match alluvial graphs to clusters viewed by UMAP
scales::show_col(scales::hue_pal()(7))


# 4. Analyse data ----------------------------------------------------------------------------------

# 4.1 Global analysis of data - how many hits do we see?

# 4.1.1a Return the number of clones we observe against each of the 197 peptides
peptide_results = table(ERGO_results_long$Peptide) %>% as.data.frame()
# peptide_results = colSums(ERGO_results_wide[9:205] >= 0.9) %>% as.data.frame() %>% rownames_to_column()  # NB: This will include 0 results - might be helpful
colnames(peptide_results) = c("peptide", "frequency")
peptide_results = arrange(peptide_results, frequency)
# We see hits against 145 of the 197 peptides

# 4.1.1b Return the number of clones we observe against each of the 58 proteins
protein_results = table(ERGO_results_long$Protein) %>% as.data.frame()
colnames(protein_results) = c("protein", "frequency")
protein_results = arrange(protein_results, frequency)
# We see hits against 54 of the 58 proteins

# 4.1.1c Tabulate the number of hits per sample based on clone bin sizes
# NB: This can work based on proportion as it is based on 'pseudo-bulk' TCR data
# Have copied this data into prism to graph as it's easier to do so

# Tabulate all results
table(ERGO_results_wide$Sample)

# Tabulate the number of hits against hyperexpanded clones (i.e. >10% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.1, ERGO_results_wide$Sample)

# Tabulate the number of hits against large clones (i.e. between 1-10% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.01 & ERGO_results_wide$Proportion < 0.1, ERGO_results_wide$Sample)

# Tabulate the number of hits against medium clones (i.e. between 0.1-1% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.001 & ERGO_results_wide$Proportion < 0.01, ERGO_results_wide$Sample)

# Tabulate the number of hits against small clones (i.e. between 0.01-0.01% of total repertoire)
table(ERGO_results_wide$Proportion >= 0.0001 & ERGO_results_wide$Proportion < 0.001, ERGO_results_wide$Sample)

# Tabulate the number of hits against single clones (i.e. single reads)
table(ERGO_results_wide$Clones == 1, ERGO_results_wide$Sample)

# 4.2 Tissue specific analysis of data - do we see disproportionate representation in one tissue?

# 4.2.1a Peptide analysis - Barplot

# Prep the data into a format suitable for plotting
peptide_results_tissue = table(ERGO_results_long$Peptide, ERGO_results_long$Tissue) %>% as.data.frame()
colnames(peptide_results_tissue) = c("peptide", "tissue", "frequency")
peptide_results_tissue = arrange(peptide_results_tissue, peptide)

# Graph the results
ggplot(peptide_results_tissue, aes(x=reorder(peptide, frequency), y=frequency, fill=tissue)) + 
  geom_bar(stat="identity") +
  xlab("peptide") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = plot_colours)
ggsave("Output/Figures/ERGO/peptide_results_tissue.png", width = 10, height = 5)

# 4.2.1b Peptide analysis - Venn Diagram

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long$Peptide, ERGO_results_long$Tissue) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
logi_results = rbind(logi_results, tissue)
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results,
             filename = "Output/Figures/ERGO/peptide_venn_diagram.png",
             col = c("red", "blue"),
             fill = c("red", "blue"),
             output=TRUE
)

# 4.2.2a Protein analysis - Barplot

# Prep the data into a format suitable for plotting
protein_results_tissue = table(ERGO_results_long$Protein, ERGO_results_long$Tissue) %>% as.data.frame()
colnames(protein_results_tissue) = c("protein", "tissue", "frequency")
protein_results_tissue = arrange(protein_results_tissue, protein)

# Graph the results
ggplot(protein_results_tissue, aes(x=reorder(protein, frequency), y=frequency, fill=tissue)) + 
  geom_bar(stat="identity") +
  xlab("protein") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = plot_colours)
ggsave("Output/Figures/ERGO/protein_results_tissue.png", width = 10, height = 5)

# 4.2.2b Protein analysis - Venn Diagram

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long$Protein, ERGO_results_long$Tissue) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
logi_results = rbind(logi_results, tissue)
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results,
             filename = "Output/Figures/ERGO/protein_venn_diagram.png",
             col = c("red", "blue"),
             fill = c("red", "blue"),
             output=TRUE
)

# 4.3 Sample specific analysis of data - do we see disproportionate representation in any particular sample/s?

# 4.3.1a Peptide analysis

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long$Peptide, ERGO_results_long$Sample) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
logi_results = rbind(logi_results, tissue)
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
temp = upset(fromList(logi_results), order.by = "freq", nsets = 8)
pdf(file = "Output/Figures/ERGO/peptide_overlap.pdf")
temp
dev.off()

# 4.3.1b Protein analysis

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long$Protein, ERGO_results_long$Sample) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
logi_results = rbind(logi_results, tissue)
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
temp = upset(fromList(logi_results), order.by = "freq", nsets = 8)
pdf(file = "Output/Figures/ERGO/protein_overlap.pdf")
temp
dev.off()


# 4.4 in-depth analysis of the 10 most overly represented peptides/proteins

# 4.4a Peptides

# Extract the names of the top 10 most represented peptides
top_10_peptides = peptide_results %>% arrange(desc(frequency)) %>% head(10)
top_10_peptides = as.character(top_10_peptides$peptide)

# Filter our data to only look at clonotypes that are potentially reactive against these peptides
top_10_peptides = subset(ERGO_results_long, subset = ERGO_results_long$Peptide %in% top_10_peptides)
top_10_peptides_list = split(top_10_peptides, f = top_10_peptides$Peptide)

# Summarise the number of clones specific on a per sample basis and plot
top_10_peptides_plot = top_10_peptides %>% group_by(Peptide, Sample) %>% summarise(Clones = sum(Clones)) %>% mutate(rank = sum(Clones)) %>% ungroup() %>% arrange(rank) %>% mutate(Peptide=factor(Peptide), Peptide=fct_reorder(Peptide, rank)) 

ggplot(top_10_peptides_plot, mapping = aes(x = Clones, y = Peptide, fill = Sample)) +
geom_bar(stat = "identity")
ggsave("Output/Figures/ERGO/top_10_peptides.png", width = 10, height = 5)


# 4.4b Proteins

# Extract the names of the top 10 most represented proteins
top_10_proteins = protein_results %>% arrange(desc(frequency)) %>% head(10)
top_10_proteins = as.character(top_10_proteins$protein)

# Filter our data to only look at clonotypes that are potentially reactive against these proteins
top_10_proteins = subset(ERGO_results_long, subset = ERGO_results_long$Protein %in% top_10_proteins)
top_10_proteins_list = split(top_10_proteins, f = top_10_proteins$Protein)

# Summarise the number of clones specific on a per sample basis and plot
top_10_proteins_plot = top_10_proteins %>% group_by(Protein, Sample) %>% summarise(Clones = sum(Clones)) %>% mutate(rank = sum(Clones)) %>% ungroup() %>% arrange(rank) %>% mutate(Protein=factor(Protein), Protein=fct_reorder(Protein, rank)) 

ggplot(top_10_proteins_plot, mapping = aes(x = Clones, y = Protein, fill = Sample)) +
  geom_bar(stat = "identity")
ggsave("Output/Figures/ERGO/top_10_proteins.png", width = 10, height = 5)


# 4.5 More in-depth analysis - focus on clones binned medium or greater. This will allow us to query the transcriptional state of our VDJ data using scRepertoire

# Let's have a look at overlap of medium & greater expansions:
med_expansions = ERGO_results_wide %>% filter(Proportion >1e-3)
med_expansions = as.character(med_expansions$TRB)
med_expansions_wide = filter(ERGO_results_wide, TRB %in% med_expansions) %>% select_if(~ !any(is.na(.)))
med_expansions_long = filter(ERGO_results_long, TRB %in% med_expansions) %>% select_if(~ !any(is.na(.)))

# Strip all binding information just keep clone size and sample information
med_expansions_trim = med_expansions_wide[-c(4:200)]
med_expansions_trim = med_expansions_trim[-c(6:12)]
# Looks like we have binding data for 14 clonotypes - 11 of which are shared between tissue for individual patients - as we only have paired transcriptome data for 2 patients, focus on those: Leaves 7 clonotypes, 1 unique to the bm

# Limit our med_expansions_trim to focus on patients 43 and 63
med_expansions_long_seurat = filter(med_expansions_long, Sample %in% c("BM43", "PB43", "BM63", "PB63"))
med_expansions_trim = filter(med_expansions_trim, Sample %in% c("BM43", "PB43", "BM63", "PB63"))
med_expansions_trim = distinct(med_expansions_trim, TRB, TRBV, TRBJ)

# Fix clonal bin sizes in our Seurat object NB: not sure why this needs to be repeated but without it "Single" and "Small" appear out of order
slot(x, "meta.data")$cloneType <- factor(slot(x, "meta.data")$cloneType, 
                                         levels = c("Hyperexpanded (628 < X <= 6275)", 
                                                    "Large (63 < X <= 628)", 
                                                    "Medium (6 < X <= 63)", 
                                                    "Small (1 < X <= 6)", 
                                                    "Single (0 < X <= 1)", NA))

# Backup x
x_backup = x

# Subset out clonotypes of interest:



# Subset x to remove Single clones and non-matched alpha data for prettier alluvial graphs
# This is simply because we're trying to link back pseudo-bulk to single cell
x = subset(x, subset = cloneType != "Single (0 < X <= 1)")
x = subset(x, subset = CTaa != str_detect(x$CTaa, "NA_")) #This isn't currently working - Syntax?

# Pull out the Seurat object metadata so we can link back our pseudobulk data back to single cell
temp = x@meta.data

# Define function to plot all 7 clonotypes to which we have paired transcriptomic data
alluvial_subset = function(TRB, TRBV, TRBJ) {
  y = temp %>% filter(str_detect(x@meta.data$CTaa, TRB)) # Get all clonotypes with a specific CDR3B
  y = y %>% filter(str_detect(y$CTgene, TRBV)) # self explanatory
  y = y %>% filter(str_detect(y$CTgene, TRBJ))# self explanatory
  y = subset(x, subset = CTaa %in% y$CTaa)
  
  # Use scRepertoire to show cluster distribution of dominant clonotype
  alluvialClonotypes(y, cloneCall = "aa", 
                     y.axes = c("orig.ident", "cluster", "cloneType"), 
                     color = "cluster") +
    scale_fill_manual(values = c("0" = "#F8766D", "1" = "#C49A00", "2" = "#53B400", "3" = "#00C094", "4" = "#00B6EB", "5" = "#A58AFF", "6" = "#FB61D7")) # NB: Current version of scRepertoire (v1.0.2)  needs specifically to call library(ggalluvial) to provide labels - author will update package to negate this in a future version.
  
  paste("Output/Figures/ERGO/alluvial", TRB, TRBV, TRBJ, ".png") %>% ggsave()
}

# Run function
pwalk(med_expansions_trim, alluvial_subset)


# 5. Write out results for further analysis --------------------------------------------------------

# save.image("Output/ERGO/ERGO_analysis.RData")

print("#># Finished running 'ERGO-II analysis' script")

# Addendum

# Pullout individual clones
# PT43
CASSIWGTSNQPQHF = subset(x, subset = CTaa == "CILRDNGGKSTF_CASSIWGTSNQPQHF")
CASSIWGTSNQPQHF_bm = subset(CASSIWGTSNQPQHF, subset = tissue == "bm")
CASSIWGTSNQPQHF_pb = subset(CASSIWGTSNQPQHF, subset = tissue == "pb")

CASNDRGTDTQYF = subset(x, subset = CTaa == "CAVRSLYNFNKFYF_CASNDRGTDTQYF")
CASNDRGTDTQYF_bm = subset(CASNDRGTDTQYF, subset = tissue == "bm")
CASNDRGTDTQYF_pb = subset(CASNDRGTDTQYF, subset = tissue == "pb")

CASRPSSGRGYNEQFF = subset(x, subset = CTaa == "CLVGDNTNAGKSTF_CASRPSSGRGYNEQFF")
CASRPSSGRGYNEQFF_bm = subset(CASRPSSGRGYNEQFF, subset = tissue == "bm")
CASRPSSGRGYNEQFF_pb = subset(CASRPSSGRGYNEQFF, subset = tissue == "pb")

# PT63
CASSPRDFWETQYF = subset(x, subset = CTaa == "CAYRSATHDMRF_CASSPRDFWETQYF")

CASIGGTYLL = subset(x, subset = CTaa == "CAASSPSNDYKLSF_CASIGGTYLL")
CASIGGTYLL_bm = subset(CASIGGTYLL, subset = tissue == "bm")
CASIGGTYLL_pb = subset(CASIGGTYLL, subset = tissue == "pb")

CASRVAEGAYDPAFF = subset(x, subset = CTaa == "CAGAFRGSNYKLTF_CASRVAEGAYDPAFF")
CASRVAEGAYDPAFF_bm = subset(CASRVAEGAYDPAFF, subset = tissue == "bm")
CASRVAEGAYDPAFF_pb = subset(CASRVAEGAYDPAFF, subset = tissue == "pb")

# This appears to be an alpha/beta/beta?
CASSSLDRPGDGYTF = subset(x, subset = CTaa == "CASSSLDRPGDGYTF_CASSSLDRPGDGYTF;CASSATRPGQLNSPLHF")
CASSSLDRPGDGYTF_bm = subset(CASSSLDRPGDGYTF, subset = tissue == "bm")
CASSSLDRPGDGYTF_pb = subset(CASSSLDRPGDGYTF, subset = tissue == "pb")
