# Script information -------------------------------------------------------------------------------

# title: Complete ERGO-II workflow and analysis of bulk TCRB data
# author: James Favaloro and Samuel Gardiner
# date: 15/07/2021
# description: In this script we will trial our analysis pipeline on bulk TCRB data generated as
# part of Deniss' honours project. This data includes 3/4 of the same patients as our 10x run
# as well as Jurkats and 2x age-matched peripheral blood controls.

# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'ERGO-II workflow' script")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scRepertoire)
  library(immunarch)
  library(httr)
  library(VennDiagram)
  library(UpSetR)
  library(ggalluvial)
})



# 2. Load data -------------------------------------------------------------------------------------

# 2.1 Load list of peptides
peptide_list = read_csv("Data/peptide_list.csv")
peptides = peptide_list$Peptide %>% sort()


# 2.2 Load Bulk VDJ data. NB: this is data is NOT pre-filtered for productive rearrangements only
clonotypes = repLoad("Data/BulkVDJ/")

# 2.3 - work out how to filter and apply filter???



# 3. Create *.csv for upload to ERGO-II ------------------------------------------------------------

# Add "Sample" Column to immunarch data
clonotypes$data$BM031.clonotypes.TRB$Sample = "BM31"
clonotypes$data$BM043.clonotypes.TRB$Sample = "BM43"
clonotypes$data$BM063.clonotypes.TRB$Sample = "BM63"
clonotypes$data$HD02_10.clonotypes.TRB$Sample = "HD02"
clonotypes$data$HD04_10.clonotypes.TRB$Sample = "HD04"
clonotypes$data$Jurkat_10.clonotypes.TRB$Sample = "Jurkat"
clonotypes$data$PB031.clonotypes.TRB$Sample = "PB31"
clonotypes$data$PB043.clonotypes.TRB$Sample = "PB43"
clonotypes$data$PB063.clonotypes.TRB$Sample = "PB63"

# Add "Tissue" Column to immunarch data
clonotypes$data$BM031.clonotypes.TRB$Tissue = "BM"
clonotypes$data$BM043.clonotypes.TRB$Tissue = "BM"
clonotypes$data$BM063.clonotypes.TRB$Tissue = "BM"
clonotypes$data$HD02_10.clonotypes.TRB$Tissue = "PB"
clonotypes$data$HD04_10.clonotypes.TRB$Tissue = "PB"
clonotypes$data$Jurkat_10.clonotypes.TRB$Tissue = "Cell_line"
clonotypes$data$PB031.clonotypes.TRB$Tissue = "PB"
clonotypes$data$PB043.clonotypes.TRB$Tissue = "PB"
clonotypes$data$PB063.clonotypes.TRB$Tissue = "PB"

# Add "Disease" Column to immunarch data
clonotypes$data$BM031.clonotypes.TRB$Disease = "MM"
clonotypes$data$BM043.clonotypes.TRB$Disease = "MM"
clonotypes$data$BM063.clonotypes.TRB$Disease = "MM"
clonotypes$data$HD02_10.clonotypes.TRB$Disease = "HD"
clonotypes$data$HD04_10.clonotypes.TRB$Disease = "HD"
clonotypes$data$Jurkat_10.clonotypes.TRB$Disease = "NA"
clonotypes$data$PB031.clonotypes.TRB$Disease = "MM"
clonotypes$data$PB043.clonotypes.TRB$Disease = "MM"
clonotypes$data$PB063.clonotypes.TRB$Disease = "MM"



# 3.1 Create clonotype files in the correct format for upload to ERGO-II

# 3.1.1 PT 31
pt31_clonotypes = clonotypes$data
pt31_clonotypes = pt31_clonotypes[c(1,7)] %>% rbindlist()
pt31_trb = subset(pt31_clonotypes, select = -c(1,2,3,6,8:21))
setnames(pt31_trb, "CDR3.aa", "TRB")
setnames(pt31_trb, "V.name", "TRBV")
setnames(pt31_trb, "J.name", "TRBJ")
pt31_trb$TRA = NA
pt31_trb$TRAV = NA
pt31_trb$TRAJ = NA
pt31_trb$`T-Cell-Type` = NA # NB: Need `` as R will read "-" as an operator
pt31_trb$Peptide = NA
pt31_trb$MHC = NA
pt31_trb = pt31_trb[,c(4,1,5,6,2,3,7,8,9)]
pt31_trb = as.data.frame(pt31_trb) # NB: rbindlist converts to data.table. Data needs to be in data.frame for mutate to function in section 3.2

# 3.1.2 PT 43
pt43_clonotypes = clonotypes$data
pt43_clonotypes = pt43_clonotypes[c(2,8)] %>% rbindlist()
pt43_trb = subset(pt43_clonotypes, select = -c(1,2,3,6,8:21))
setnames(pt43_trb, "CDR3.aa", "TRB")
setnames(pt43_trb, "V.name", "TRBV")
setnames(pt43_trb, "J.name", "TRBJ")
pt43_trb$TRA = NA
pt43_trb$TRAV = NA
pt43_trb$TRAJ = NA
pt43_trb$`T-Cell-Type` = NA # NB: Need `` as R will read "-" as an operator
pt43_trb$Peptide = NA
pt43_trb$MHC = NA
pt43_trb = pt43_trb[,c(4,1,5,6,2,3,7,8,9)]
pt43_trb = as.data.frame(pt43_trb) # NB: rbindlist converts to data.table. Data needs to be in data.frame for mutate to function in section 3.2

# 3.1.3 PT 63
pt63_clonotypes = clonotypes$data
pt63_clonotypes = pt63_clonotypes[c(3,9)] %>% rbindlist()
pt63_trb = subset(pt63_clonotypes, select = -c(1,2,3,6,8:21))
setnames(pt63_trb, "CDR3.aa", "TRB")
setnames(pt63_trb, "V.name", "TRBV")
setnames(pt63_trb, "J.name", "TRBJ")
pt63_trb$TRA = NA
pt63_trb$TRAV = NA
pt63_trb$TRAJ = NA
pt63_trb$`T-Cell-Type` = NA # NB: Need `` as R will read "-" as an operator
pt63_trb$Peptide = NA
pt63_trb$MHC = NA
pt63_trb = pt63_trb[,c(4,1,5,6,2,3,7,8,9)]
pt63_trb = as.data.frame(pt63_trb) # NB: rbindlist converts to data.table. Data needs to be in data.frame for mutate to function in section 3.2

# 3.1.4 Controls
controls_clonotypes = clonotypes$data
controls_clonotypes = controls_clonotypes[c(4,5,6)] %>% rbindlist()
controls_trb = subset(controls_clonotypes, select = -c(1,2,3,6,8:21))
setnames(controls_trb, "CDR3.aa", "TRB")
setnames(controls_trb, "V.name", "TRBV")
setnames(controls_trb, "J.name", "TRBJ")
controls_trb$TRA = NA
controls_trb$TRAV = NA
controls_trb$TRAJ = NA
controls_trb$`T-Cell-Type` = NA # NB: Need `` as R will read "-" as an operator
controls_trb$Peptide = NA
controls_trb$MHC = NA
controls_trb = controls_trb[,c(4,1,5,6,2,3,7,8,9)]
controls_trb = as.data.frame(controls_trb) # NB: rbindlist converts to data.table. Data needs to be in data.frame for mutate to function in section 3.2

# 3.2 Copy individual peptides to allow generation of *.csv for upload to ERGO-II

# 3.2.1 PT 31
ERGO_pt31 <- map(peptides,
                     function(peptide) {
                       peptide_df <- pt31_trb %>% 
                         mutate(Peptide = peptide) # NB: Doesn't work on data.tables
                     }) %>%
  set_names(peptides)

# 3.2.2 PT 43
ERGO_pt43 <- map(peptides,
                     function(peptide) {
                       peptide_df <- pt43_trb %>% 
                         mutate(Peptide = peptide) # NB: Doesn't work on data.tables
                     }) %>%
  set_names(peptides)

# 3.2.3 PT 63
ERGO_pt63 <- map(peptides,
                     function(peptide) {
                       peptide_df <- pt63_trb %>% 
                         mutate(Peptide = peptide) # NB: Doesn't work on data.tables
                     }) %>%
  set_names(peptides)

# 3.2.4 Controls
ERGO_controls <- map(peptides,
                 function(peptide) {
                   peptide_df <- controls_trb %>% 
                     mutate(Peptide = peptide) # NB: Doesn't work on data.tables
                 }) %>%
  set_names(peptides)

# 3.3 Write out all *.csv for in the correct format ready for upload to ERGO

# 3.3.1a Create the folder for PT31 and point R to it
cidr = getwd()
mkfolder = "Data/ERGO/PT31/Upload/"
dir.create(file.path(cidr, mkfolder), recursive = T)

# 3.3.1b Write out the *.csv for PT31
imap(ERGO_pt31, function(df, name){
  write_csv(df, file = paste0(mkfolder, name, "_pt31.csv"), na = "")
})

# 3.3.2a Create the folder for PT43 and point R to it
cidr = getwd()
mkfolder = "Data/ERGO/PT43/Upload/"
dir.create(file.path(cidr, mkfolder), recursive = T)

# 3.3.2b Write out the *.csv for PT43
imap(ERGO_pt43, function(df, name){
  write_csv(df, file = paste0(mkfolder, name, "_pt43.csv"), na = "")
})

# 3.3.3a Create the folder for PT63 and point R to it
cidr = getwd()
mkfolder = "Data/ERGO/PT63/Upload/"
dir.create(file.path(cidr, mkfolder), recursive = T)

# 3.3.3b Write out the *.csv for PT63
imap(ERGO_pt63, function(df, name){
  write_csv(df, file = paste0(mkfolder, name, "_pt63.csv"), na = "")
})

# 3.3.4a Create the folder for Controls and point R to it
cidr = getwd()
mkfolder = "Data/ERGO/Controls/Upload/"
dir.create(file.path(cidr, mkfolder), recursive = T)

# 3.3.4b Write out the *.csv for Controls
imap(ERGO_controls, function(df, name){
  write_csv(df, file = paste0(mkfolder, name, "_controls.csv"), na = "")
})



# 4. Run the ERGO_script_webtool script ------------------------------------------------------------

# Set webtool address
target = "http://tcr2.cs.biu.ac.il"

# 4.1 Run web tool for patient 31

# 4.1.1 Create directory for results and point R to the exported *.csv
cidr = getwd()
mkfolder = "Data/ERGO/PT31/Results/"
dir.create(file.path(cidr, mkfolder), recursive = T)
filenames <- list.files("Data/ERGO/PT31/Upload/") #Point to the directory where *.csv files are located

# 4.1.2 Run the webtool
# For every file name, make a HTTP POST request
walk(filenames, function(x) {
  Sys.sleep(5) # Wait 5 seconds before making the next request
  message(paste("Requesting", x))
  request = POST(target,
                 body = list(
                   # Emulate the manual form fields
                   input_file = upload_file(paste0("Data/ERGO/PT31/Upload/", x)), #Point to the directory where *.csv files are located
                   
                   model_type = "AE", # AE: Autoencoder, LSTM: LSTM
                   dataset = "mcpas", #mcpas: McPAS, vdjdb: VDJdb
                   use_vj = "true" # Other arguments can be added here including use_alpha, use_mhc & use_t_type
                 ),
                 # Use the SSWAHS proxy for each request - this fixes the firewall issue with hospital computers. Hash out if running from outside hospital.
                 # config = (use_proxy("sswproxy5.sswahs.nsw.gov.au", port = 8080))
  )
  # Write out the HTTP response without parsing it,
  write_file(content(request, as = "raw"),
             file = paste0(mkfolder, x)) # Write out the ERGO-II results
}
)

# 4.2 Run web tool for patient 43

# 4.2.1 Create directory for results and point R to the exported *.csv
cidr = getwd()
mkfolder = "Data/ERGO/PT43/Results/"
dir.create(file.path(cidr, mkfolder), recursive = T)
filenames <- list.files("Data/ERGO/PT43/Upload/") #Point to the directory where *.csv files are located

# 4.2.2 Run the webtool
# For every file name, make a HTTP POST request
walk(filenames, function(x) {
  Sys.sleep(5) # Wait 5 seconds before making the next request
  message(paste("Requesting", x))
  request = POST(target,
                 body = list(
                   # Emulate the manual form fields
                   input_file = upload_file(paste0("Data/ERGO/PT43/Upload/", x)), #Point to the directory where *.csv files are located
                   
                   model_type = "AE", # AE: Autoencoder, LSTM: LSTM
                   dataset = "mcpas", #mcpas: McPAS, vdjdb: VDJdb
                   use_vj = "true" # Other arguments can be added here including use_alpha, use_mhc & use_t_type
                 ),
                 # Use the SSWAHS proxy for each request - this fixes the firewall issue with hospital computers. Hash out if running from outside hospital.
                 # config = (use_proxy("sswproxy5.sswahs.nsw.gov.au", port = 8080))
  )
  # Write out the HTTP response without parsing it,
  write_file(content(request, as = "raw"),
             file = paste0(mkfolder, x)) # Write out the ERGO-II results
}
)


# 4.3 Run web tool for patient 63

# 4.3.1 Create directory for results and point R to the exported *.csv
cidr = getwd()
mkfolder = "Data/ERGO/PT63/Results/"
dir.create(file.path(cidr, mkfolder), recursive = T)
filenames <- list.files("Data/ERGO/PT63/Upload/") #Point to the directory where *.csv files are located

# 4.3.2 Run the webtool
# For every file name, make a HTTP POST request
walk(filenames, function(x) {
  Sys.sleep(5) # Wait 5 seconds before making the next request
  message(paste("Requesting", x))
  request = POST(target,
                 body = list(
                   # Emulate the manual form fields
                   input_file = upload_file(paste0("Data/ERGO/PT63/Upload/", x)), #Point to the directory where *.csv files are located
                   
                   model_type = "AE", # AE: Autoencoder, LSTM: LSTM
                   dataset = "mcpas", #mcpas: McPAS, vdjdb: VDJdb
                   use_vj = "true" # Other arguments can be added here including use_alpha, use_mhc & use_t_type
                 ),
                 # Use the SSWAHS proxy for each request - this fixes the firewall issue with hospital computers. Hash out if running from outside hospital.
                 # config = (use_proxy("sswproxy5.sswahs.nsw.gov.au", port = 8080))
  )
  # Write out the HTTP response without parsing it,
  write_file(content(request, as = "raw"),
             file = paste0(mkfolder, x)) # Write out the ERGO-II results
}
)

# 4.4 Run web tool for Controls

# 4.4.1 Create directory for results and point R to the exported *.csv
cidr = getwd()
mkfolder = "Data/ERGO/Controls/Results/"
dir.create(file.path(cidr, mkfolder), recursive = T)
filenames <- list.files("Data/ERGO/Controls/Upload/") #Point to the directory where *.csv files are located

# 4.4.2 Run the webtool
# For every file name, make a HTTP POST request
walk(filenames, function(x) {
  Sys.sleep(5) # Wait 5 seconds before making the next request
  message(paste("Requesting", x))
  request = POST(target,
                 body = list(
                   # Emulate the manual form fields
                   input_file = upload_file(paste0("Data/ERGO/Controls/Upload/", x)), #Point to the directory where *.csv files are located
                   
                   model_type = "AE", # AE: Autoencoder, LSTM: LSTM
                   dataset = "mcpas", #mcpas: McPAS, vdjdb: VDJdb
                   use_vj = "true" # Other arguments can be added here including use_alpha, use_mhc & use_t_type
                 ),
                 # Use the SSWAHS proxy for each request - this fixes the firewall issue with hospital computers. Hash out if running from outside hospital.
                 # config = (use_proxy("sswproxy5.sswahs.nsw.gov.au", port = 8080))
  )
  # Write out the HTTP response without parsing it,
  write_file(content(request, as = "raw"),
             file = paste0(mkfolder, x)) # Write out the ERGO-II results
}
)



# 5. Read in the ERGO-II results -------------------------------------------------------------------

# 5.1 Read in results for PT 31
mydir = "Data/ERGO/PT31/Results/"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
pt31_results = map_dfr(myfiles,
                           read_csv)
pt31_results_long = left_join(pt31_results, peptide_list, by = "Peptide")
pt31_results_wide <- pt31_results %>%
  pivot_wider(names_from = Peptide, values_from = Score) %>%
  unnest(cols = where(is.list))

# 5.2 Read in results for PT 43
mydir = "Data/ERGO/pt43/Results/"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
pt43_results = map_dfr(myfiles,
                       read_csv)
pt43_results_long = left_join(pt43_results, peptide_list, by = "Peptide")
pt43_results_wide <- pt43_results %>%
  pivot_wider(names_from = Peptide, values_from = Score) %>%
  unnest(cols = where(is.list))

# 5.3 Read in results for PT 63
mydir = "Data/ERGO/pt63/Results/"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
pt63_results = map_dfr(myfiles,
                       read_csv)
pt63_results_long = left_join(pt63_results, peptide_list, by = "Peptide")
pt63_results_wide <- pt63_results %>%
  pivot_wider(names_from = Peptide, values_from = Score) %>%
  unnest(cols = where(is.list))

# 5.4 Read in results for Controls
mydir = "Data/ERGO/Controls/Results/"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
controls_results = map_dfr(myfiles,
                           read_csv)
controls_results_long = left_join(controls_results, peptide_list, by = "Peptide")
controls_results_wide <- controls_results %>%
  pivot_wider(names_from = Peptide, values_from = Score) %>%
  unnest(cols = where(is.list))


# 6. Combine the datasets to allow for analysis ----------------------------------------------------

# 6.1 Analysis of TCR data independent of transcriptome

# 6.1.1 Gather all ERGO-II results into one convenient place

# Create a Monster list of all patient results
temp = full_join(pt31_results_long, pt43_results_long) %>% full_join(pt63_results_long) %>% full_join(controls_results_long)
temp2 = full_join(pt31_results_wide, pt43_results_wide) %>% full_join(pt63_results_wide) %>% full_join(controls_results_wide)

# Trim the Monster list to Score >= 0.9 and round scores to 2 decimal places
temp = temp %>% filter_at(vars(10), any_vars(. >=0.9)) %>% mutate_if(is.numeric, round, digits = 2)
temp2 = temp2 %>% filter_at(vars(9:ncol(temp2)), any_vars(. >=0.9)) %>% mutate_if(is.numeric, round, digits = 2)


# 6.1.2 Gather all clonotype data into one convenient place

# Create a Monster list of all patient results
temp3 = rbindlist(clonotypes$data)


# 6.1.3 Merge the two datasets

# Left join the two datasets to create a file for analysis
ERGO_results_long = left_join(temp, temp3, by = c("TRB" = "CDR3.aa", "TRBV" = "V.name", "TRBJ" = "J.name")) %>% unique()

ERGO_results_wide = left_join(temp2, temp3, by = c("TRB" = "CDR3.aa", "TRBV" = "V.name", "TRBJ" = "J.name")) %>% unique()
# I think this might be right for this data(?) Drops from 4853 to 4071

# Apply a filter to remove single read clonotypes
ERGO_results_long_filtered = subset(ERGO_results_long, subset = ERGO_results_long$Clones != "1")
ERGO_results_wide_filtered = subset(ERGO_results_wide, subset = ERGO_results_wide$Clones != "1")

# Apply a filter to limit our analysis to PB samples
ERGO_results_long_pb = subset(ERGO_results_long, subset = ERGO_results_long$Tissue == "PB")
ERGO_results_wide_pb = subset(ERGO_results_wide, subset = ERGO_results_wide$Tissue == "PB")

# Apply a filter to limit our analysis to filtered PB samples
ERGO_results_long_pb_filtered = subset(ERGO_results_long_filtered, subset = ERGO_results_long_filtered$Tissue == "PB")
ERGO_results_wide_pb_filtered = subset(ERGO_results_wide_filtered, subset = ERGO_results_wide_filtered$Tissue == "PB")

# Apply a filter to limit our analysis to medium and greater expansions
med_expansions = ERGO_results_wide %>% filter(Proportion >1e-3)
med_expansions = as.character(med_expansions$TRB)
ERGO_results_long_med_expansions = filter(ERGO_results_long, TRB %in% med_expansions) %>% select_if(~ !any(is.na(.)))
ERGO_results_wide_med_expansions = filter(ERGO_results_wide, TRB %in% med_expansions) %>% select_if(~ !any(is.na(.)))

# Apply a filter to limit our analysis to medium and greater expansions in PB samples
ERGO_results_long_pb_med_expansions = subset(ERGO_results_long_med_expansions, subset = ERGO_results_long_med_expansions$Tissue == "PB")
ERGO_results_wide_pb_med_expansions = subset(ERGO_results_wide_med_expansions, subset = ERGO_results_wide_med_expansions$Tissue == "PB")

# 7. Write out results for safe keeping ------------------------------------------------------------

# Write out these file for safekeeping
mkfolder = "Output/ERGO/"
write.csv(ERGO_results_long, "Output/ERGO/ERGO_results_bulk_long.csv")
write.csv(ERGO_results_wide, "Output/ERGO/ERGO_results_bulk_wide.csv")
write.csv(ERGO_results_long_filtered, "Output/ERGO/ERGO_results_bulk_long_filtered.csv")
write.csv(ERGO_results_wide_filtered, "Output/ERGO/ERGO_results_bulk_wide_filtered.csv")
write.csv(ERGO_results_long_pb, "Output/ERGO/ERGO_results_bulk_long_pb.csv")
write.csv(ERGO_results_wide_pb, "Output/ERGO/ERGO_results_bulk_wide_pb.csv")
write.csv(ERGO_results_long_pb_filtered, "Output/ERGO/ERGO_results_bulk_long_pb_filtered.csv")
write.csv(ERGO_results_wide_pb_filtered, "Output/ERGO/ERGO_results_bulk_wide_pb_filtered.csv")
write.csv(ERGO_results_long_med_expansions, "Output/ERGO/ERGO_results_bulk_long_med_expansions.csv")
write.csv(ERGO_results_wide_med_expansions, "Output/ERGO/ERGO_results_bulk_wide_med_expansions.csv")
write.csv(ERGO_results_long_pb_med_expansions, "Output/ERGO/ERGO_results_bulk_long_pb_med_expansions.csv")
write.csv(ERGO_results_wide_pb_med_expansions, "Output/ERGO/ERGO_results_bulk_wide_pb_med_expansions.csv")


# Clean up the environment
rm(list=ls()[! ls() %in% c("clonotypes", "ERGO_results_long", "ERGO_results_wide", "ERGO_results_long_filtered", "ERGO_results_wide_filtered", "ERGO_results_long_pb", "ERGO_results_wide_pb", "ERGO_results_long_pb_filtered", "ERGO_results_wide_pb_filtered", "ERGO_results_long_med_expansions", "ERGO_results_wide_med_expansions", "ERGO_results_long_pb_med_expansions", "ERGO_results_wide_pb_med_expansions")])

# Save out the results for further analysis
save.image("Output/ERGO/ERGO_results_bulk.RData")

# 8. Analysis loop ---------------------------------------------------------------------------------

# NB: 16/7/21 - this is currently moved out into seperate scripts just to get the images - work with Sam to tidy this up into a function that works over a list

# 8.1 Adjust environment for plotting data
plot_colours = c("blue", "red")
plot_colours_bulk = c("red", "green", "blue")
tissue = c("BM", "BM", "BM", "BM", "PB", "PB", "PB", "PB")

# 8.2 Define analysis functions


# Run the analysis loop
lapply(df, plot_func)




# 4. Analyse data ----------------------------------------------------------------------------------

# 4.1 Global analysis of data - how many hits do we see?

# 4.1.1a Return the number of clones we observe against each of the 197 peptides
peptide_results = table(ERGO_results_long$Peptide) %>% as.data.frame()
# peptide_results = colSums(ERGO_results_wide[9:205] >= 0.9) %>% as.data.frame() %>% rownames_to_column()  # NB: This will include 0 results - might be helpful
colnames(peptide_results) = c("peptide", "frequency")
peptide_results = arrange(peptide_results, frequency)
table(peptide_results)
# We see hits against 184 of the 197 peptides

# CAN ASK IF WE SEE A DIFFERENCE HERE BETWEEN NORMAL AND MM

# 4.1.1b Return the number of clones we observe against each of the 58 proteins
protein_results = table(ERGO_results_long$Protein) %>% as.data.frame()
colnames(protein_results) = c("protein", "frequency")
protein_results = arrange(protein_results, frequency)
table(protein_results)
# We see hits against 57 of the 58 proteins

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
  scale_fill_manual(values = plot_colours_bulk)
ggsave("Output/Figures/ERGO/peptide_results_bulk_tissue.png", width = 10, height = 5)

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
             filename = "Output/Figures/ERGO/peptide_venn_diagram_bulk.png",
             col = c("red", "green", "blue"),
             fill = c("red", "green", "blue"),
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
  scale_fill_manual(values = plot_colours_bulk)
ggsave("Output/Figures/ERGO/protein_results_bulk_tissue.png", width = 10, height = 5)

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
             filename = "Output/Figures/ERGO/protein_venn_diagram_bulk.png",
             col = c("red", "green", "blue"),
             fill = c("red", "green", "blue"),
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
pdf(file = "Output/Figures/ERGO/peptide_overlap_bulk.pdf")
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
pdf(file = "Output/Figures/ERGO/protein_overlap_bulk.pdf")
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
ggsave("Output/Figures/ERGO/top_10_peptides_bulk.png", width = 10, height = 5)


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
ggsave("Output/Figures/ERGO/top_10_proteins_bulk.png", width = 10, height = 5)


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

# MM vs. HD

# 5. Write out results for further analysis --------------------------------------------------------

# save.image("Output/ERGO/ERGO_analysis.RData")

print("#># Finished running 'ERGO-II analysis' script")