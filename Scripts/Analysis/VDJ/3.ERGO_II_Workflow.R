# Script information -------------------------------------------------------------------------------

# Title: Complete ERGO-II workflow
# Author: James Favaloro and Samuel Gardiner
# Date: 2023-08-21
# Description: In this script clonotype data will be queried using the webtool form of the 
# deep-learning peptide prediction tool, ERGO-II (http://tcr2.cs.biu.ac.il/). We will use the 
# results from Walz et al (2015). (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4608392/) who found 
# 197 peptides from 58 proteins to be overexpressed on MHC class-I in both myeloma cell lines and
# primary samples.

# NB: This script is compatible with any list of peptides provided in *.csv format and any TCR 
# provided in Immunarch compatible format (https://immunarch.com/)
# NB: Using clonotypes provided by immunarch will not discriminate between multiple observations
# i.e. if a clonotype is seen in both alpha/beta and beta only these will be considered as a single
# observation - this will result in an overestimation of more prominent clones but allows for
# flexibility in the script and compatibility with both bulk and single cell data.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running 'ERGO-II workflow' script")

suppressPackageStartupMessages({
  librarian::shelf(Seurat, tidyverse, scRepertoire, immunarch, httr)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# 2.1 Load list of peptides
peptide_list <- read_csv("Data/Public_data/VDJ/peptide_list.csv")
peptides <- peptide_list$Peptide %>% sort()


# 2.2 Load 10x VDJ data as single chain data. 
# NB: this is data has been pre-filtered for productive rearrangements only
clonotypes <- repLoad("Data/10x/TCR_clean", .mode = "single")


# 2.3 Load processed 10x data
all_list <- readRDS("Data/R_out/10x/all_list.rds")


# 3. Pre-processing --------------------------------------------------------------------------------

# 3.1 Create *.csv for upload to ERGO-II webtool ---------------------------------------------------

# 3.1.A Add additional identifying data to clonotypes list and remove unnecessary data -------------
# Add "Sample" Column to immunarch data
clonotypes$data$`1_1-BM13_TRB`$Sample <- "BM13"
clonotypes$data$`1_2-PB13_TRB`$Sample <- "PB13"
clonotypes$data$`2_1-BM31_TRB`$Sample <- "BM31"
clonotypes$data$`2_2-PB31_TRB`$Sample <- "PB31"
clonotypes$data$`3_1-BM43_TRB`$Sample <- "BM43"
clonotypes$data$`3_2-PB43_TRB`$Sample <- "PB43"
clonotypes$data$`4_1-BM63_TRB`$Sample <- "BM63"
clonotypes$data$`4_2-PB63_TRB`$Sample <- "PB63"

# Add "Tissue" Column to immunarch data
clonotypes$data$`1_1-BM13_TRB`$Tissue <- "BM"
clonotypes$data$`1_2-PB13_TRB`$Tissue <- "PB"
clonotypes$data$`2_1-BM31_TRB`$Tissue <- "BM"
clonotypes$data$`2_2-PB31_TRB`$Tissue <- "PB"
clonotypes$data$`3_1-BM43_TRB`$Tissue <- "BM"
clonotypes$data$`3_2-PB43_TRB`$Tissue <- "PB"
clonotypes$data$`4_1-BM63_TRB`$Tissue <- "BM"
clonotypes$data$`4_2-PB63_TRB`$Tissue <- "PB"

# Drop TRA data
clonotypes <- clonotypes$data[-c(1, 3, 5, 7, 9, 11, 13, 15)]

# 3.1.B Create patient clonotype files in the correct format for upload to ERGO-II -----------------
# NB: rbindlist converts to data.table. Data needs to be in data.frame for mutate to function

# 3.1.B.i PT 13 & PT 31 ----------------------------------------------------------------------------
pt13pt31_clonotypes <- clonotypes
pt13pt31_clonotypes <- pt13pt31_clonotypes[-c(5,6,7,8)] %>% rbindlist()
pt13pt31_trb <- subset(pt13pt31_clonotypes, select = -c(1,2,3,6,8:21))
setnames(pt13pt31_trb, "CDR3.aa", "TRB")
setnames(pt13pt31_trb, "V.name", "TRBV")
setnames(pt13pt31_trb, "J.name", "TRBJ")
pt13pt31_trb$Tissue <- NULL
pt13pt31_trb$TRA <- NA
pt13pt31_trb$TRAV <- NA
pt13pt31_trb$TRAJ <- NA
pt13pt31_trb$`T-Cell-Type` <- NA # NB: Need `` as R will read "-" as an operator
pt13pt31_trb$Peptide <- NA
pt13pt31_trb$MHC <- NA
pt13pt31_trb <- pt13pt31_trb[,c(4,1,5,6,2,3,7,8,9)]
pt13pt31_trb <- as.data.frame(pt13pt31_trb)

# 3.1.B.ii PT 43 & PT 63 ---------------------------------------------------------------------------
pt43pt63_clonotypes <- clonotypes
pt43pt63_clonotypes <- pt43pt63_clonotypes[-c(1,2,3,4)] %>% rbindlist()
pt43pt63_trb <- subset(pt43pt63_clonotypes, select = -c(1,2,3,6,8:21))
setnames(pt43pt63_trb, "CDR3.aa", "TRB")
setnames(pt43pt63_trb, "V.name", "TRBV")
setnames(pt43pt63_trb, "J.name", "TRBJ")
pt43pt63_trb$Tissue <- NULL
pt43pt63_trb$TRA <- NA
pt43pt63_trb$TRAV <- NA
pt43pt63_trb$TRAJ <- NA
pt43pt63_trb$`T-Cell-Type` <- NA
pt43pt63_trb$Peptide <- NA
pt43pt63_trb$MHC <- NA
pt43pt63_trb <- pt43pt63_trb[,c(4,1,5,6,2,3,7,8,9)]
pt43pt63_trb <- as.data.frame(pt43pt63_trb)

# 3.1.C Copy individual peptides to allow generation of *.csv for upload to ERGO-II ----------------

# 3.1.C.i PT 13 & PT 31 ----------------------------------------------------------------------------
ERGO_pt13pt31 <- map(peptides,
                     function(peptide) {
                       peptide_df <- pt13pt31_trb %>% 
                         mutate(Peptide = peptide) # NB: Doesn't work on data.tables
                     }) %>%
  set_names(peptides)

# 3.1.C.ii PT 43 & PT 63 ---------------------------------------------------------------------------
ERGO_pt43pt63 <- map(peptides,
                     function(peptide) {
                       peptide_df <- pt43pt63_trb %>% 
                         mutate(Peptide = peptide) # NB: Doesn't work on data.tables
                     }) %>%
  set_names(peptides)

# 3.1.D Write out all *.csv for in the correct format ready for upload to ERGO ---------------------

# 3.1.D.i Create the folder for PT13 & PT31 and point R to it --------------------------------------
cidr <- getwd()
mkfolder <- "Data/R_out/ERGO/PT13_PT31/Upload/"
dir.create(file.path(cidr, mkfolder), recursive = T)

# 3.1.D.ii Write out the *.csv for PT13 & PT31 -----------------------------------------------------
imap(ERGO_pt13pt31, function(df, name){
  write_csv(df, file = paste0(mkfolder, name, "_pt13pt31.csv"), na = "")
})

# 3.1.D.iii Create the folder for PT43 & PT63 and point R to it ------------------------------------
cidr <- getwd()
mkfolder <- "Data/R_out/ERGO/PT43_PT63/Upload/"
dir.create(file.path(cidr, mkfolder), recursive = T)

# 3.1.D.iv Write out the *.csv for PT43 & PT63 -----------------------------------------------------
imap(ERGO_pt43pt63, function(df, name){
  write_csv(df, file = paste0(mkfolder, name, "_pt43pt63.csv"), na = "")
})


# 4. Run the ERGO_script_webtool script ------------------------------------------------------------

# Set webtool address
target <- "https://tcr2.cs.biu.ac.il"

# 4.1 Run web tool for patients 13 & 31 ------------------------------------------------------------

# 4.1.1 Create directory for results and point R to the exported *.csv -----------------------------
cidr <- getwd()
mkfolder <- "Data/R_out/ERGO/PT13_PT31/Results/"
dir.create(file.path(cidr, mkfolder), recursive = T)
filenames <- list.files("Data/R_out/ERGO/PT13_PT31/Upload/")

# 4.1.2 Run the webtool ----------------------------------------------------------------------------
# For every file name, make a HTTP POST request
walk(filenames, function(x) {
  Sys.sleep(5) # Wait 5 seconds before making the next request
  message(paste("Requesting", x))
  request = POST(target,
                 body = list(
                   # Emulate the manual form fields
                   input_file = upload_file(paste0("Data/R_out/ERGO/PT13_PT31/Upload/", x)), 
                   model_type = "AE", # AE: Autoencoder, LSTM: LSTM
                   dataset = "mcpas", # mcpas: McPAS, vdjdb: VDJdb
                   use_vj = "true" # Can add use_alpha, use_mhc & use_t_type if desired
                 ),
  )
  # Write out the HTTP response without parsing it,
  write_file(content(request, as = "raw"),
             file = paste0(mkfolder, x)) # Write out the ERGO-II results
}
)


# 4.2 Run web tool for patients 43 & 63 ------------------------------------------------------------

# 4.2.1 Create directory for results and point R to the exported *.csv -----------------------------
cidr <- getwd()
mkfolder <- "Data/R_out/ERGO/PT43_PT63/Results/"
dir.create(file.path(cidr, mkfolder), recursive = T)
filenames <- list.files("Data/R_out/ERGO/PT43_PT63/Upload/")

# 4.2.2 Run the webtool ----------------------------------------------------------------------------
# For every file name, make a HTTP POST request
walk(filenames, function(x) {
  Sys.sleep(5) # Wait 5 seconds before making the next request
  message(paste("Requesting", x))
  request = POST(target,
                 body = list(
                   # Emulate the manual form fields
                   input_file = upload_file(paste0("Data/R_out/ERGO/PT43_PT63/Upload/", x)),
                   model_type = "AE", # AE: Autoencoder, LSTM: LSTM
                   dataset = "mcpas", # mcpas: McPAS, vdjdb: VDJdb
                   use_vj = "true" # Can add use_alpha, use_mhc & use_t_type if desired
                 ),
  )
  # Write out the HTTP response without parsing it,
  write_file(content(request, as = "raw"),
             file = paste0(mkfolder, x)) # Write out the ERGO-II results
}
)


# 5. Read in the ERGO-II results -------------------------------------------------------------------

# 5.1 Read in results for PT 13 & 31 ---------------------------------------------------------------
mydir <- "Data/R_out/ERGO/PT13_PT31/Results/"
myfiles <- list.files(path = mydir, pattern = "*.csv", full.names = TRUE)
pt13pt31_results <- map_dfr(myfiles, read_csv)
pt13pt31_results_long <- left_join(pt13pt31_results, peptide_list, by = "Peptide")
pt13pt31_results_wide <- pt13pt31_results %>%
  pivot_wider(names_from = Peptide, values_from = Score) %>%
  unnest(cols = where(is.list))


# 5.2 Read in results for PT 43 & 63 ---------------------------------------------------------------
mydir <- "Data/R_out/ERGO/PT43_PT63/Results/"
myfiles <- list.files(path = mydir, pattern = "*.csv", full.names = TRUE)
pt43pt63_results <- map_dfr(myfiles, read_csv)
pt43pt63_results_long <- left_join(pt43pt63_results, peptide_list, by = "Peptide")
pt43pt63_results_wide <- pt43pt63_results %>%
  pivot_wider(names_from = Peptide, values_from = Score) %>%
  unnest(cols = where(is.list))


# 6. Combine the datasets to allow for analysis ----------------------------------------------------

# 6.1 Analysis of TCR data independent of transcriptome --------------------------------------------

# 6.1.1 Gather all ERGO-II results into one convenient place ---------------------------------------

# Create a Monster list of all patient results
temp <- full_join(pt13pt31_results_long, pt43pt63_results_long)
temp2 <- full_join(pt13pt31_results_wide, pt43pt63_results_wide)

# Trim the Monster list to Score >= 0.9 and round scores to 2 decimal places
temp <- temp %>% filter_at(vars(10), any_vars(. >=0.9)) %>% mutate_if(is.numeric, round, digits = 2)
temp2 <- temp2 %>% filter_at(vars(9:ncol(temp2)), any_vars(. >=0.9)) %>% 
  mutate_if(is.numeric, round, digits = 2)

# 6.1.2 Gather all clonotype data into one convenient place ----------------------------------------

# Create a Monster list of all patient results
temp3 <- rbindlist(clonotypes)

# 6.1.3 Merge the two datasets  --------------------------------------------------------------------

# Left join the two datasets to create a file for analysis
ERGO_results_long <- 
  left_join(temp, temp3, by = c("TRB" = "CDR3.aa", "TRBV" = "V.name", "TRBJ" = "J.name")) %>% 
  unique()

ERGO_results_wide <- 
  left_join(temp2, temp3, by = c("TRB" = "CDR3.aa", "TRBV" = "V.name", "TRBJ" = "J.name"))
ERGO_results_wide <- 
  ERGO_results_wide[!duplicated(ERGO_results_wide$RawClonotypeID),] 
# Drops from 1218 to 969 - should only be one row per clonotype per tissue - i.e. if a clone is 
# shared between PB and BM there should be 2 results that only differ based on sample
# NB: Immunarch has been updated - raw_clonotype_id is now RawClonotypeID


# 6.2 Analysis of TCR data with transcriptome ------------------------------------------------------

# Copy results for pt43 and pt63 into temp objects
temp <- pt43pt63_results_long
temp2 <- pt43pt63_results_wide

# Trim results to Score >= 0.9 and round scores to 2 decimal places
temp <- temp %>% filter_at(vars(10), any_vars(. >=0.9)) %>% mutate_if(is.numeric, round, digits = 2)
temp2 <- temp2 %>% 
  filter_at(vars(9:ncol(temp2)), any_vars(. >=0.9)) %>% mutate_if(is.numeric, round, digits = 2)

# 6.2.1 Extract proportion of cells from our seurat object
seurat_clonotypes <- table(all_list$all_10x$CTaa, all_list$all_10x$orig.ident) %>% 
  as.data.frame.matrix() %>% rownames_to_column()

# Extract the CDR3 from the CTaa and replace it
clonotypes_trb <- gsub(".*_", "", seurat_clonotypes$rowname) %>% as.array()
seurat_clonotypes$rowname <- clonotypes_trb

# Rename 1st column to CDR3
seurat_clonotypes <- seurat_clonotypes %>% rename(CDR3 = 1)

# Remove all cells to which we are missing TCRB
seurat_clonotypes <- subset(seurat_clonotypes, subset = seurat_clonotypes$CDR3 != "NA")

# 6.2.2 # Left join the two datasets to create a file for analysis
# NB: As this data is treated as psuedobulk, at this time we are considering all clonotypes 
# clonal irrespective of whether or not they are observed as paired or unpaired.
ERGO_results_long_seurat <- left_join(temp, seurat_clonotypes, by = c("TRB" = "CDR3")) %>% unique()
ERGO_results_wide_seurat <- left_join(temp2, seurat_clonotypes, by = c("TRB" = "CDR3")) %>% unique()


# 7. Save out data and print stats -----------------------------------------------------------------

write.csv(ERGO_results_long, "Output/ERGO/ERGO_results_long.csv")
write.csv(ERGO_results_wide, "Output/ERGO/ERGO_results_wide.csv")
write.csv(ERGO_results_long_seurat, "Output/ERGO/ERGO_results_long_seurat.csv")
write.csv(ERGO_results_wide_seurat, "Output/ERGO/ERGO_results_wide_seurat.csv")

# Report time here as environment will be cleaned prior to saving
end_time <- Sys.time()
end_time - start_time

# Clean up the environment
rm(list=ls()[! ls() %in% c("clonotypes", "ERGO_results_long", "ERGO_results_long_seurat", 
                           "ERGO_results_wide", "ERGO_results_wide_seurat", "seurat_clonotypes")])

# Save out the results for further analysis
save.image("Data/R_out/ERGO/ERGO_results.RData")

gc()

print("#># Finished running '3. ERGO-II workflow' script")