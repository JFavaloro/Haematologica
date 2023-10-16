# Script information -------------------------------------------------------------------------------

# Title: V(D)J - Comparative analysis and TCRMatch
# Author: James Favaloro and Samuel Gardiner
# Date: 2023-08-21
# Description: In this script we will attempt to identify public clonotypes in our dataset by 
# comparing against VDJdb and McPAS, two public databases of known antigen specificity:
# (https://vdjdb.cdr3.net/search & http://friedmanlab.weizmann.ac.il/McPAS-TCR/. As these databases 
# are based predominantly off bulk datasets, we will treat our data as if it were bulk. This work 
# includes manipulations not performed in R, chiefly, the VDJdb and McPAS datasets have been 
# manually curated to exclude non human data and split into alpha and beta datasets for ease of 
# manipulation. TCRMatch has been run manually by uploading the top 10 CDR3 Beta chain data for each
# sample and downloading the results as a *.tsv file - placed in Data/Public_data/VDJ.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '2. V(D)J - Comparative analysis' script")

suppressPackageStartupMessages({
  librarian::shelf(tidyverse, immunarch, VennDiagram)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load V(D)J data in immunarch format 
# NB: These are renamed filtered_contig_annotation *.CSV files for ease of loading into R
immdata <- repLoad("Data/10x/TCR", .mode = "single")

# Load the VDJdb public datasets
VDJdb_tra <- dbLoad("Data/Public_data/VDJ/20201112_VDJdb_TRA.txt", .db = "vdjdb")
VDJdb_trb <- dbLoad("Data/Public_data/VDJ/20201112_VDJdb_TRB.txt", .db = "vdjdb")

# Load McPAS public dataset
McPAS_tra <- dbLoad("Data/Public_data/VDJ/20201030_McPAS_TCR_filtered_modified_TRA.tsv", 
                    .db = "vdjdb")
McPAS_trb <- dbLoad("Data/Public_data/VDJ/20201030_McPAS_TCR_filtered_modified_TRB.tsv", 
                    .db = "vdjdb")

# Load the TCRMatch results
# NB: These results are from database query dated 27MAY2022
TCRMatch_names <- list.files("Data/Public_data/VDJ/TCRMatch", pattern = "*.tsv", full.names = TRUE)
TCRMatch_results <- lapply(TCRMatch_names,read_tsv)
TCRMatch_names <- 
  gsub("Data/Public_data/VDJ/TCRMatch/tcrmatch_result_", "", gsub(".tsv", "", TCRMatch_names))
names(TCRMatch_results) <- TCRMatch_names
TCRMatch_results <- map(TCRMatch_results, ~ (.x %>% select(-last_col())))

# Define a colour scheme for plotting TCRMatch results
plot_cols = c("red", "blue", "green", "purple", "black")


# 3. Pre-processing --------------------------------------------------------------------------------

# Remove multi and non productive clones
immdata_top_10 <- immdata$data[-c(1,4,7,8,11,14,17,20,23)]

# Restrict data to the top 10 clonotypes observable per sample
for (i in names(immdata_top_10)) {
  immdata_top_10[[`i`]] <- head(immdata_top_10[[`i`]], n = 10)
}

# Split the top 10 list into alpha and beta lists
immdata_top_10_alpha <- immdata_top_10[-c(2,4,6,8,10,12,14,16)]
immdata_top_10_beta <- immdata_top_10[-c(1,3,5,7,9,11,13,15)]

# Create a character vector of the top 10 CDR3 for beta chain data for feeding TCRMatch
top_10_CDR3B <- list()

for (i in names(immdata_top_10_beta)) {
  temp = immdata_top_10_beta[[`i`]][["CDR3.aa"]]
  top_10_CDR3B[[length(top_10_CDR3B) +1]] <- temp
}

top_10_CDR3B <- top_10_CDR3B %>% flatten() %>%as.character() %>% as_tibble()
write_csv(top_10_CDR3B, "top_10_CDR3.csv", col_names = FALSE)

# Remove SARS-CoV2 entries from the TCRMatch_results list
for (i in names(TCRMatch_results)) {
  TCRMatch_results[[`i`]] <- 
    TCRMatch_results[[`i`]][ - grep("SARS", TCRMatch_results[[`i`]][["organism"]]), ]
}


# 4. Analyse ---------------------------------------------------------------------------------------

# 4.1 Compare to VDJdb dataset ---------------------------------------------------------------------
VDJdb_alpha_results <- 
  dbAnnotate(immdata_top_10_alpha, VDJdb_tra, 
             c("CDR3.aa", "V.name", "J.name"), 
             c("CDR3.aa", "V.Name", "J.Name")) %>% as_tibble() %>% filter(Samples > 2)
VDJdb_beta_results <- 
  dbAnnotate(immdata_top_10_beta, VDJdb_trb, 
             c("CDR3.aa", "V.name", "J.name"), 
             c("CDR3.aa", "V.Name", "J.Name")) %>% as_tibble() %>% filter(Samples > 2)


# 4.2 Compare to MCPas dataset ---------------------------------------------------------------------
McPAS_alpha_results <- 
  dbAnnotate(immdata_top_10_alpha, McPAS_tra, 
             c("CDR3.aa", "V.name", "J.name"), 
             c("CDR3.aa", "V.Name", "J.Name")) %>% as_tibble() %>% filter(Samples > 2)
McPAS_beta_results <- 
  dbAnnotate(immdata_top_10_beta, McPAS_trb, 
             c("CDR3.aa", "V.name", "J.name"), 
             c("CDR3.aa", "V.Name", "J.Name")) %>% as_tibble() %>% filter(Samples > 2)


# 4.3 Analyse TCRMatch results ---------------------------------------------------------------------

# Create empty list to store results
TCR_match_highest_score <- list()

# Run loop to extract highest scoring hit
for (i in names(TCRMatch_results)) {
  temp <- TCRMatch_results[[`i`]] %>% arrange(desc(score)) %>% group_by(input_sequence) %>% slice(1)
  TCR_match_highest_score[[length(TCR_match_highest_score) +1]] <- temp
}

# Fix names of the TCR_match_highest_score list
names(TCR_match_highest_score) <- names(TCRMatch_results)

# Run a loop to add flank amino acids back to TCR_match_highest_score results
for (i in names(TCR_match_highest_score)) {
  TCR_match_highest_score[[`i`]][["input_sequence"]] <- 
    paste0("C", TCR_match_highest_score[[`i`]][["input_sequence"]], "F")
}

# Link the results of TCR_Match back to the top 10 list
TCR_match_highest_score <- map2(immdata_top_10_beta, TCR_match_highest_score, function(.x,.y) {
  left_join(.x,.y, by = c("CDR3.aa" = "input_sequence"), all = TRUE) %>% arrange(desc(Proportion))
})

# Adjust score value to allow plotting and move Proportion column to left and *100, 
for (i in names(TCR_match_highest_score)) {
  TCR_match_highest_score[[`i`]][["score"]] <- 
    TCR_match_highest_score[[`i`]][["score"]]-0.9
  TCR_match_highest_score[[`i`]][["score"]] <- 
    TCR_match_highest_score[[`i`]][["score"]]*10
  TCR_match_highest_score[[`i`]][["Proportion"]] <- 
    TCR_match_highest_score[[`i`]][["Proportion"]]*100 
  TCR_match_highest_score[[`i`]] <- TCR_match_highest_score[[`i`]] %>% 
    select(Proportion, everything())
}

# Add Sample column to TCR_match_highest_score
for (i in 1:length(TCR_match_highest_score)){
  TCR_match_highest_score[[`i`]] <- cbind(TCR_match_highest_score[[`i`]], Sample = TCRMatch_names[i])
}

# Add Rank column to TCR_match_highest_score
for (i in 1:length(TCR_match_highest_score)){
  TCR_match_highest_score[[`i`]] <- 
    TCR_match_highest_score[[`i`]] %>% rowid_to_column("Rank")
}

# Move all relevant data into a dataframe
TCRMatch_plot <- bind_rows(TCR_match_highest_score)
TCRMatch_plot <- TCRMatch_plot[,-c(9:21)]

# Fix NA values in the score column to 0
TCRMatch_plot$score[is.na(TCRMatch_plot$score)] <- 0

# Standardise the names in the "organism" column
TCRMatch_plot$organism[grepl("Infl", TCRMatch_plot$organism)] <- "Influenza A"
TCRMatch_plot$organism[grepl("s 5", TCRMatch_plot$organism)] <- "CMV"
TCRMatch_plot$organism[grepl("s 4", TCRMatch_plot$organism)] <- "EBV"
TCRMatch_plot$organism[grepl("B virus", TCRMatch_plot$organism)] <- "Other"
TCRMatch_plot$organism[grepl("lus", TCRMatch_plot$organism)] <- "Other"
TCRMatch_plot$organism[grepl("ens", TCRMatch_plot$organism)] <- "Other"
TCRMatch_plot$organism[grepl("type I", TCRMatch_plot$organism)] <- "Other"
TCRMatch_plot$organism[grepl("wheat", TCRMatch_plot$organism)] <- "Other"
TCRMatch_plot$organism[grepl("Yel", TCRMatch_plot$organism)] <- "Other"

# Tidy up column names for prettier plotting
names(TCRMatch_plot)[10] <- "Score"
names(TCRMatch_plot)[14] <- "Antigen"

# Rearrange data to allow paired samples to plot next to each other
pt_order <- c("BM13", "PB13", "BM31", "PB31", "BM43", "PB43", "BM63", "PB63")

TCRMatch_plot <- TCRMatch_plot %>%
  mutate(Sample = factor(Sample, levels = pt_order)) %>%
  arrange(Sample)

# Plot data as a bubble graph
ggplot(TCRMatch_plot,
       aes(x = Sample, y = Rank, size = Proportion)) +
  # Bottom layer: a coloured, solid circle with alpha
  geom_point(aes(color = Antigen, alpha = Score)) +
  # Top layer: an outer, hollow circle with a coloured outline:
  geom_point(shape = "circle open", mapping = aes(color = Antigen), stroke = 2) +
  # Explicitly set na.value with alpha 00:
  scale_colour_manual(values = plot_cols, na.value = "#000000FF") +
  scale_size(range = c(1,15), breaks = c(1, 5, 10)) +
  scale_y_continuous(breaks = 1:10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(), 
        text = element_text(size = 18), axis.text = element_text(size = 22))
ggsave(filename = paste0("Output/Figures/VDJ/10x/Final/TCRMatch.png"), 
       units = "cm", height = 30, width = 30, bg = "transparent")


# 5. Save out data and print stats -----------------------------------------------------------------

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '2. Comparative analysis' script")