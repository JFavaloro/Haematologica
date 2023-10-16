# Script information -------------------------------------------------------------------------------

# Title: Control dataset #1 - Pre-processing
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will load a control dataset and associated metadata and prepare
# the data such that we reduce the number of genes in the dataset to those found in our 10x data.
# This dataset comprises approximate age-matched control BM and PB CD3+ T-cells from four donors:
# BM: 1 (65/M), 2 (52/M) and PB: A (51/M), B (53/M). from Szabo et al., 2019:
# https://doi.org/10.1038/s41467-019-12464-3. The authors have provided additional metadata, making
# data manipulation much more manageable, allowing us to remove CD4+ cells and create a control 
# Seurat object in the next script. Although the samples are not paired, we will treat them as such 
# for ease of downstream utility.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '1. Control dataset #1 - Pre-processing' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, data.table, openxlsx, R.utils)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# 2.1 Use fred to read the count matrices ----------------------------------------------------------
controls <- list.files("Data/Public_data/Controls_1", pattern = "*.gz", full.names = TRUE)
raw_data <- lapply(controls,fread, fill = TRUE, header = TRUE)
controls <- gsub("Data/Public_data/Controls_1/", "", gsub(".txt.gz", "", controls))
names(raw_data) <- controls
raw_data <- map(raw_data, ~ (.x %>% select(-last_col())))

# 2.2 Read in metadata file and separate out for individual donors and cell states -----------------
# NB: This metadata file has been manually curated from Szabo et al., 2019 Source Data 
# (Supplementary Figure 6) to select for BM and PB CD8+ T-cells.

# Read in metadata file and rename BL to PB for consistency
# NB: Ideally we'd read in barcodes as rownames, however as we have duplicates this is not an
# option. This will result in a warning later on but works as intended.
metadata <- read.xlsx("Data/Public_data/Controls_1/Control_CD8_barcodes.xlsx") %>%
  filter(Tissue %in% c("BL", "BM"))
metadata$Tissue <- str_replace_all(metadata$Tissue, "BL", "PB")

metadata_list <- metadata %>%
  group_by(orig.ident, Tissue, stimulation_status) %>%
  group_split()
names(metadata_list) <- controls


# 2.3 Read in 10x data for creating an accession list & create Seurat objects ----------------------

# Create a character vector of all our samples
all_samples <- apply(expand.grid(c("BM", "PB"), c("13","31","43","63")), 1, paste, collapse="")

# Import GEX data with ENSG rather than gene
all_gex <- lapply(all_samples, 
                  function(ID)  Read10X(data.dir = 
                                          str_glue('Data/10x/GEX/{ID}/filtered_feature_bc_matrix'), 
                                        gene.column = 1)) 
names(all_gex) <- all_samples

# Create Seurat Objects with default parameters:
# min.cells = gene must be observed in at least 3 cells
# min.features = cells must express at least 200 different genes
all_seurat <- lapply(1:length(all_samples),
                     function(x) CreateSeuratObject(counts = 
                                                      all_gex[[x]], 
                                                    project = str_to_lower(all_samples[x]), 
                                                    min.cells = 3, min.features = 200))
names(all_seurat) <- all_samples


# 3. Tidying up control dataset --------------------------------------------------------------------

# 3.1 Create a list of genes seen in our dataset using Accession rather than gene name -------------

# Create an empty list to store all genes observed in our 10x data
accession_list <- list()

# Retrieve all genes from 10x data into a list
for (i in names(all_seurat)){
  temp <- list(all_seurat[[i]]@assays[["RNA"]]@counts@Dimnames[[1]])
  accession_list[[length(accession_list) +1]] <- temp
}
accession_list <- accession_list %>% unlist() %>% unique() %>% 
  as_tibble() %>% select(Accession = value)

# 3.2 Remove suffix ".##" from Accession column of control dataset ---------------------------------
# NB: The ENSG nomenclature is Gene_ID.version - our dataset is missing the version and so we need
# to remove it - the Gene_ID is SET

# Define function to remove version number (suffix ".##") from ENSG in control dataset
remove_suffix <- function(raw_data) {
  mutate(raw_data, across(.cols = c(Accession),
                          ~gsub(pattern = "\\.[0-9]+", replacement = "", x = .x)))  
}

# Apply function to raw data appending the suffix "_clean" to each item in the list
# NB: We won't change the names for ease of merging the data with metadata later on
raw_data_clean = lapply(raw_data, remove_suffix) # %>% set_names(paste0(names(raw_data), "_clean"))
# Or the other way round: (prefix code) paste0("clean_", names(raw_data))


# 3.3 Reduce number of genes in control dataset to those observed in our dataset using Accession ---

# Define function for removing genes
remove_genes <- function(input_raw){
  input_raw <- input_raw[input_raw$Accession %in% accession_list$Accession,]
  return(input_raw)
}

# Apply remove_genes function to raw_data_clean
raw_data_clean <- lapply(raw_data_clean, remove_genes)


# 3.4 Remove the Accession column from the control dataset -----------------------------------------

# Define function to remove Accession column
drop_accession <- function(input_raw){
  input_raw <- input_raw[,-1]
  return(input_raw)
}

# Apply function to raw_clean_data
raw_data_clean <- lapply(raw_data_clean, drop_accession)


# 4. Save out data and print stats -----------------------------------------------------------------

saveRDS(raw_data_clean, "Data/R_out/Controls_1/raw_data_clean.rds")
saveRDS(metadata_list, "Data/R_out/Controls_1/metadata_list.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '1. Control dataset #1 - Pre-processing' script")