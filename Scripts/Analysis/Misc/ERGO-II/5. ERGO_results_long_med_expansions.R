# X.1 Adjust environment for plotting data - Ask Sam if there's a better way to do this
plot_colours = c("red", "blue")
plot_colours_br = c("blue", "red")
plot_colours_rgb = c("red", "green", "blue")
plot_colours_brg = c("blue", "red", "green")
pt_colours = scale_fill_manual(values = c("BM31" = "#F8766D", "BM43" = "#D39200", "BM63" = "#93AA00", "HD02" = "#00BA38", "HD04" = "#00C19F", "Jurkat" = "#00B9E3", "PB31" = "#619CFF", "PB43" = "#DB72FB", "PB63" = "#FF61C3")) # Ask sam if we can set this as an environment - pretty sure we can


# 1. Peptide analysis

# Set our global analysis to allow plotting later
peptide_results = table(ERGO_results_long_med_expansions$Peptide) %>% as.data.frame()
colnames(peptide_results) = c("peptide", "frequency")
peptide_results = arrange(peptide_results, frequency)

# 1.1 Peptide analysis - Barplot

# 1.1a Barplot - Tissue

# Prep the data into a format suitable for plotting
peptide_results_tissue = table(ERGO_results_long_med_expansions$Peptide, ERGO_results_long_med_expansions$Tissue) %>% as.data.frame()
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

# 1.1b Barplot - Disease

# Prep the data into a format suitable for plotting
peptide_results_disease = table(ERGO_results_long_med_expansions$Peptide, ERGO_results_long_med_expansions$Disease) %>% as.data.frame()
colnames(peptide_results_disease) = c("peptide", "disease", "frequency")
peptide_results_disease = arrange(peptide_results_disease, peptide)

# Graph the results
ggplot(peptide_results_disease, aes(x=reorder(peptide, frequency), y=frequency, fill=disease)) + 
  geom_bar(stat="identity") +
  xlab("peptide") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = plot_colours_br)
ggsave("Output/Figures/ERGO/peptide_results_disease.png", width = 10, height = 5)


# 1.2 Peptide analysis - Venn Diagram

# 1.2a Tissue - Venn Diagram

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long_med_expansions$Peptide, ERGO_results_long_med_expansions$Tissue) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results,
             filename = "Output/Figures/ERGO/peptide_venn_diagram_tissue.png",
             col = plot_colours,
             fill = plot_colours,
             output=TRUE
)

# 1.2b Disease - Venn Diagram

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long_med_expansions$Peptide, ERGO_results_long_med_expansions$Disease) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results,
             filename = "Output/Figures/ERGO/peptide_venn_diagram_disease.png",
             col = plot_colours_br,
             fill = plot_colours_br,
             output=TRUE
)

# 1.3 Peptide analysis - UpSetR overlap

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long_med_expansions$Peptide, ERGO_results_long_med_expansions$Sample) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
temp = upset(fromList(logi_results), order.by = "freq", nsets = 9)
pdf(file = "Output/Figures/ERGO/peptide_overlap.pdf")
temp
dev.off()


# 1.4 Top 10 Peptides by frequency

# Extract the names of the top 10 most represented peptides
top_10_peptides = peptide_results %>% arrange(desc(frequency)) %>% head(10)
top_10_peptides = as.character(top_10_peptides$peptide)

# Filter our data to only look at clonotypes that are potentially reactive against these peptides
top_10_peptides = subset(ERGO_results_long_med_expansions, subset = ERGO_results_long_med_expansions$Peptide %in% top_10_peptides)
top_10_peptides_list = split(top_10_peptides, f = top_10_peptides$Peptide)

# Summarise the number of clones specific on a per sample basis and plot
top_10_peptides_plot = top_10_peptides %>% group_by(Peptide, Sample) %>% summarise(Clones = sum(Clones)) %>% mutate(rank = sum(Clones)) %>% ungroup() %>% arrange(rank) %>% mutate(Peptide=factor(Peptide), Peptide=fct_reorder(Peptide, rank)) 

ggplot(top_10_peptides_plot, mapping = aes(x = Clones, y = Peptide, fill = Sample)) +
  geom_bar(stat = "identity")
ggsave("Output/Figures/ERGO/top_10_peptides.png", width = 10, height = 5)


------------------------
# 2. Protein analysis
  
# Set our global analysis to allow plotting later
protein_results = table(ERGO_results_long_med_expansions$Protein) %>% as.data.frame()
colnames(protein_results) = c("protein", "frequency")
protein_results = arrange(protein_results, frequency)

# 2.1 Protein analysis - Barplot

# 2.1a Barplot - Tissue

# Prep the data into a format suitable for plotting
protein_results_tissue = table(ERGO_results_long_med_expansions$Protein, ERGO_results_long_med_expansions$Tissue) %>% as.data.frame()
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

# 2.1b Barplot - Disease

# Prep the data into a format suitable for plotting
protein_results_disease = table(ERGO_results_long_med_expansions$Protein, ERGO_results_long_med_expansions$Disease) %>% as.data.frame()
colnames(protein_results_disease) = c("protein", "disease", "frequency")
protein_results_disease = arrange(protein_results_disease, protein)

# Graph the results
ggplot(protein_results_disease, aes(x=reorder(protein, frequency), y=frequency, fill=disease)) + 
  geom_bar(stat="identity") +
  xlab("protein") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = plot_colours_br)
ggsave("Output/Figures/ERGO/protein_results_disease.png", width = 10, height = 5)

# 2.2 Protein analysis - Venn Diagram

# 2.2a Tissue - Venn Diagram

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long_med_expansions$Protein, ERGO_results_long_med_expansions$Tissue) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results,
             filename = "Output/Figures/ERGO/protein_venn_diagram_tissue.png",
             col = plot_colours,
             fill = plot_colours,
             output=TRUE
)

# 2.2b Disease - Venn Diagram

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long_med_expansions$Protein, ERGO_results_long_med_expansions$Disease) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results,
             filename = "Output/Figures/ERGO/protein_venn_diagram_disease.png",
             col = plot_colours_br,
             fill = plot_colours_br,
             output=TRUE
)


# 2.3 Protein analysis - UpSetR overlap

# Prep the data into a format suitable for plotting
results = table(ERGO_results_long_med_expansions$Protein, ERGO_results_long_med_expansions$Sample) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
temp = upset(fromList(logi_results), order.by = "freq", nsets = 9)
pdf(file = "Output/Figures/ERGO/protein_overlap.pdf")
temp
dev.off()


# 2.4 Top 10 Proteins by frequency

# Extract the names of the top 10 most represented proteins
top_10_proteins = protein_results %>% arrange(desc(frequency)) %>% head(10)
top_10_proteins = as.character(top_10_proteins$protein)

# Filter our data to only look at clonotypes that are potentially reactive against these proteins
top_10_proteins = subset(ERGO_results_long_med_expansions, subset = ERGO_results_long_med_expansions$Protein %in% top_10_proteins)
top_10_proteins_list = split(top_10_proteins, f = top_10_proteins$Protein)

# Summarise the number of clones specific on a per sample basis and plot
top_10_proteins_plot = top_10_proteins %>% group_by(Protein, Sample) %>% summarise(Clones = sum(Clones)) %>% mutate(rank = sum(Clones)) %>% ungroup() %>% arrange(rank) %>% mutate(Protein=factor(Protein), Protein=fct_reorder(Protein, rank)) 

ggplot(top_10_proteins_plot, mapping = aes(x = Clones, y = Protein, fill = Sample)) +
  geom_bar(stat = "identity")
ggsave("Output/Figures/ERGO/top_10_proteins.png", width = 10, height = 5)