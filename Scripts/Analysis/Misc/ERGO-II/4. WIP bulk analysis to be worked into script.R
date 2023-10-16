# Write out a function to loop through analysis for all filtered datasets

scale_fill_manual(values = c("BM31" = "#F8766D", "BM43" = "#D39200", "BM63" = "#93AA00", "HD02" = "#00BA38", "HD04" = "#00C19F", "Jurkat" = "#00B9E3", "PB31" = "#619CFF", "PB43" = "#DB72FB", "PB63" = "#FF61C3"))

# 4.X Disease specific analysis of data - do we see more hits in MM than HD? #NB Need to fix this up

# List the dataframes we wish to analyse
df = list(ERGO_results_long = ERGO_results_long, ERGO_results_long_filtered = ERGO_results_long_filtered, ERGO_results_long_pb = ERGO_results_long_pb, ERGO_results_long_pb_filtered = ERGO_results_long_pb_filtered, ERGO_results_long_med_expansions = ERGO_results_long_med_expansions, ERGO_results_long_pb_med_expansions = ERGO_results_long_pb_med_expansions)

# Define function to plot all data including various filters
plot_func = function(df) {
# 4.2.1a Peptide analysis - Barplot

# Prep the data into a format suitable for plotting
peptide_results_disease = table(df$Peptide, df$Disease) %>% as.data.frame()
colnames(peptide_results_disease) = c("peptide", "disease", "frequency")
peptide_results_disease = arrange(peptide_results_disease, peptide)

# Graph the results
ggplot(peptide_results_disease, aes(x=reorder(peptide, frequency), y=frequency, fill=disease)) + 
  geom_bar(stat="identity") +
  xlab("peptide") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = plot_colours) %>% paste("Output/Figures/ERGO/peptide_bulk_", ".png", width = 10, height = 5) %>% ggsave()

# 4.2.1b Peptide analysis - Venn Diagram

# Graph the results - Ask Sam
venn.diagram(logi_results,
             filename = paste("Output/Figures/ERGO/peptide_venn_diagram_bulk_.png"),
             col = c("blue", "red"),
             fill = c("blue", "red"),
             output=TRUE
)

# 4.2.2a Protein analysis - Barplot

# Prep the data into a format suitable for plotting
df = filter(df, df$Tissue == "PB")
protein_results_disease = table(df$Protein, df$Disease) %>% as.data.frame()
colnames(protein_results_disease) = c("protein", "disease", "frequency")
protein_results_disease = arrange(protein_results_disease, protein)

# Graph the results
ggplot(protein_results_disease, aes(x=reorder(protein, frequency), y=frequency, fill=disease)) + 
  geom_bar(stat="identity") +
  xlab("protein") +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = plot_colours) %>% paste("Output/Figures/ERGO/protein_bulk_", ".png", width = 10, height = 5) %>% ggsave()

# 4.2.1b Protein analysis - Venn Diagram

# Prep the data into a format suitable for plotting
results = table(df$Protein, df$Disease) %>% as.data.frame.array()
logi_results = results >= 1
logi_results = as.data.frame(logi_results) *1
logi_results = rbind(logi_results, tissue)
# Adjust the logical results to a list to allow plotting
logi_results = lapply(logi_results, function(x) rownames(logi_results)[x==1])

# Graph the results
venn.diagram(logi_results,
             filename = paste("Output/Figures/ERGO/protein_venn_diagram_bulk_.png"),
             col = c("blue", "red"),
             fill = c("blue", "red"),
             output=TRUE
)


# 4.4 in-depth analysis of the 10 most overly represented peptides/proteins

# 4.4a Peptides

# Extract the names of the top 10 most represented peptides
top_10_peptides = peptide_results %>% arrange(desc(frequency)) %>% head(10)
top_10_peptides = as.character(top_10_peptides$peptide)

# Filter our data to only look at clonotypes that are potentially reactive against these peptides
top_10_peptides = subset(df, subset = df$Peptide %in% top_10_peptides)
top_10_peptides_list = split(top_10_peptides, f = top_10_peptides$Peptide)

# Summarise the number of clones specific on a per sample basis and plot
top_10_peptides_plot = top_10_peptides %>% group_by(Peptide, Sample) %>% summarise(Clones = sum(Clones)) %>% mutate(rank = sum(Clones)) %>% ungroup() %>% arrange(rank) %>% mutate(Peptide=factor(Peptide), Peptide=fct_reorder(Peptide, rank)) 

ggplot(top_10_peptides_plot, mapping = aes(x = Clones, y = Peptide, fill = Sample)) +
  geom_bar(stat = "identity") %>% paste("Output/Figures/ERGO/top_10_peptides_bulk_.png", width = 10, height = 5) %>% ggsave()


# 4.4b Proteins

# Extract the names of the top 10 most represented proteins
top_10_proteins = protein_results %>% arrange(desc(frequency)) %>% head(10)
top_10_proteins = as.character(top_10_proteins$protein)

# Filter our data to only look at clonotypes that are potentially reactive against these proteins
top_10_proteins = subset(df, subset = df$Protein %in% top_10_proteins)
top_10_proteins_list = split(top_10_proteins, f = top_10_proteins$Protein)

# Summarise the number of clones specific on a per sample basis and plot
top_10_proteins_plot = top_10_proteins %>% group_by(Protein, Sample) %>% summarise(Clones = sum(Clones)) %>% mutate(rank = sum(Clones)) %>% ungroup() %>% arrange(rank) %>% mutate(Protein=factor(Protein), Protein=fct_reorder(Protein, rank)) 

ggplot(top_10_proteins_plot, mapping = aes(x = Clones, y = Protein, fill = Sample)) +
  geom_bar(stat = "identity") %>% paste("Output/Figures/ERGO/top_10_protein_bulk_.png", width = 10, height = 5) %>% ggsave()
}

# Run function
lapply(df, plot_func)
