# 1. Import libraries ------------------------------------------------------------------------------

print("#># Start running 'ERGO-II workflow' script")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scRepertoire)
  library(immunarch)
  library(httr)
})



# 2. Load data -------------------------------------------------------------------------------------

# 2.1 Load list of peptides
peptide_list = read_csv("Data/peptide_list.csv")
peptides = peptide_list$Peptide %>% sort()


# 2.2 Load 10x VDJ data as single chain data. NB: this is data has been pre-filtered for productive rearrangements only

clonotypes = read_csv("Data/Public/AML TCR/SRR8359192.csv")
clonotypes2 = read_csv("Data/Public/AML TCR/SRR8359193.csv")
clonotypes3 = read_csv("Data/Public/AML TCR/SRR8359194.csv")
clonotypes = list(clonotypes,clonotypes2, clonotypes3)
names(clonotypes) <- c("SRR8359192", "SRR8359193", "SRR8359194")


clonotypes$SRR8359192$sample = "SRR8359192"
clonotypes$SRR8359193$sample = "SRR8359193"
clonotypes$SRR8359194$sample = "SRR8359193"

clonotypes$SRR8359192$tissue = "BM"
clonotypes$SRR8359193$tissue = "BM"
clonotypes$SRR8359194$tissue = "BM"

# Control_BM
clonotypes = clonotypes %>% rbindlist()
setnames(clonotypes, "cdr3", "TRB")
setnames(clonotypes, "v_gene", "TRBV")
setnames(clonotypes, "j_gene", "TRBJ")
clonotypes$TRA = NA
clonotypes$TRAV = NA
clonotypes$TRAJ = NA
clonotypes$`T-Cell-Type` = NA # NB: Need `` as R will read "-" as an operator
clonotypes$Peptide = NA
clonotypes$MHC = NA
clonotypes_up = clonotypes[,c(8,1,9,10,2,4,11,12,13)]
clonotypes_up = as.data.frame(clonotypes_up) # NB: rbindlist converts to data.table. Data needs to be in data.frame for mutate to function in section 3.2

# 3.2 Copy individual peptides to allow generation of *.csv for upload to ERGO-II

ERGO_clonotypes <- map(peptides,
                     function(peptide) {
                       peptide_df <- clonotypes_up %>% 
                         mutate(Peptide = peptide) # NB: Doesn't work on data.tables
                     }) %>%
  set_names(peptides)

# 3.3 Write out all *.csv for in the correct format ready for upload to ERGO

cidr = getwd()
mkfolder = "Data/ERGO/Control_AML/Upload/"
dir.create(file.path(cidr, mkfolder), recursive = T)


imap(ERGO_clonotypes, function(df, name){
  write_csv(df, file = paste0(mkfolder, name, "_clonotypes.csv"), na = "")
})

# 4. Run the ERGO_script_webtool script ------------------------------------------------------------

# Set webtool address
target = "https://tcr2.cs.biu.ac.il"

# 4.1 Run web tool

# 4.1.1 Create directory for results and point R to the exported *.csv
cidr = getwd()
mkfolder = "Data/ERGO/Control_AML/Results/"
dir.create(file.path(cidr, mkfolder), recursive = T)
filenames <- list.files("Data/ERGO/Control_AML/Upload/") #Point to the directory where *.csv files are located

# 4.1.2 Run the webtool
# For every file name, make a HTTP POST request
walk(filenames, function(x) {
  Sys.sleep(5) # Wait 5 seconds before making the next request
  message(paste("Requesting", x))
  request = POST(target,
                 body = list(
                   # Emulate the manual form fields
                   input_file = upload_file(paste0("Data/ERGO/Control_AML/Upload/", x)), #Point to the directory where *.csv files are located
                   
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

# 5.1 Read in results
mydir = "Data/ERGO/Control_AML/Results/"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
control_results = map_dfr(myfiles,
                           read_csv)
control_results_long = left_join(control_results, peptide_list, by = "Peptide")
control_results_wide <- control_results %>%
  pivot_wider(names_from = Peptide, values_from = Score) %>%
  unnest(cols = where(is.list))


# 6.1.1 Gather all ERGO-II results into one convenient place

# Create a Monster list of all patient results
temp = control_results_long
temp2 = control_results_wide

# Trim the Monster list to Score >= 0.9 and round scores to 2 decimal places
temp = temp %>% filter_at(vars(10), any_vars(. >=0.9)) %>% mutate_if(is.numeric, round, digits = 2)
temp2 = temp2 %>% filter_at(vars(9:ncol(temp2)), any_vars(. >=0.9)) %>% mutate_if(is.numeric, round, digits = 2)


# 6.1.2 Gather all clonotype data into one convenient place

# Create a Monster list of all patient results
temp3 = clonotypes


# 6.1.3 Merge the two datasets

# Left join the two datasets to create a file for analysis
ERGO_results_long = left_join(temp, temp3, by = c("TRB" = "TRB", "TRBV" = "TRBV", "TRBJ" = "TRBJ")) %>% unique()

ERGO_results_wide = left_join(temp2, temp3, by = c("TRB" = "TRB", "TRBV" = "TRBV", "TRBJ" = "TRBJ"))

# ERGO_results_wide = ERGO_results_wide[!duplicated(ERGO_results_wide$raw_clonotype_id),] # Drops from 1218 to 963 - should only be one row per clonotype per tissue - i.e. if a clone is shared between PB and BM there should be 2 results that only differ based on sample

# Write out these file for safekeeping
mkfolder = "Output/ERGO/"
write.csv(ERGO_results_long, "Output/ERGO/control_ERGO_results_long.csv")
write.csv(ERGO_results_wide, "Output/ERGO/control_ERGO_results_wide.csv")


# 4.1 Global analysis of data - how many hits do we see?

# 4.1.1a Return the number of clones we observe against each of the 197 peptides
peptide_results = table(ERGO_results_long$Peptide.x) %>% as.data.frame() # Why have the column names changed?
# peptide_results = colSums(ERGO_results_wide[9:205] >= 0.9) %>% as.data.frame() %>% rownames_to_column()  # NB: This will include 0 results - might be helpful
colnames(peptide_results) = c("peptide", "frequency")
peptide_results = arrange(peptide_results, frequency)
# We see hits against 78 of the 197 peptides

# 4.1.1b Return the number of clones we observe against each of the 58 proteins
protein_results = table(ERGO_results_long$Protein) %>% as.data.frame()
colnames(protein_results) = c("protein", "frequency")
protein_results = arrange(protein_results, frequency)
# We see hits against 41 of the 58 proteins

# 4.1.1c Tabulate the number of hits per sample based on clone bin sizes
# NB: This can work based on proportion as it is based on 'pseudo-bulk' TCR data
# Have copied this data into prism to graph as it's easier to do so

# Tabulate all results
table(ERGO_results_wide$sample)

# Tabulate the number of hits against hyperexpanded clones (i.e. >10% of total repertoire)
table(ERGO_results_wide$proportion >= 0.1, ERGO_results_wide$sample)

# Tabulate the number of hits against large clones (i.e. between 1-10% of total repertoire)
table(ERGO_results_wide$proportion >= 0.01 & ERGO_results_wide$proportion < 0.1, ERGO_results_wide$sample)

# Tabulate the number of hits against medium clones (i.e. between 0.1-1% of total repertoire)
table(ERGO_results_wide$proportion >= 0.001 & ERGO_results_wide$proportion < 0.01, ERGO_results_wide$sample)

# Tabulate the number of hits against small clones (i.e. between 0.01-0.01% of total repertoire)
table(ERGO_results_wide$proportion >= 0.0001 & ERGO_results_wide$proportion < 0.001, ERGO_results_wide$sample)

# 4.5 More in-depth analysis - focus on clones binned medium or greater. This will allow us to query the transcriptional state of our VDJ data using scRepertoire

# Let's have a look at overlap of medium & greater expansions:
med_expansions = ERGO_results_wide %>% filter(proportion >1e-3)
med_expansions = as.character(med_expansions$TRB)
med_expansions_wide = filter(ERGO_results_wide, TRB %in% med_expansions) %>% select_if(~ !any(is.na(.)))
med_expansions_long = filter(ERGO_results_long, TRB %in% med_expansions) %>% select_if(~ !any(is.na(.)))

# Strip all binding information just keep clone size and sample information
med_expansions_trim = med_expansions_wide[-c(4:200)]
med_expansions_trim = med_expansions_trim[-c(6:12)]
# Looks like we have binding data for 14 clonotypes - 11 of which are shared between tissue for individual patients - as we only have paired transcriptome data for 2 patients, focus on those: Leaves 7 clonotypes, 1 unique to the bm