# Script information -------------------------------------------------------------------------------

# Title: Load data
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will load the 10x 5'GEX data as well as the V(D)J data. We will do 
# some quick pre-processing of the V(D)J data using the scRepertoire R package. As the GEX and V(D)J
# data will be required in the analysis scripts, we will save this data out in the R format.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '1. Load data' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# 2.1 Load GEX data --------------------------------------------------------------------------------

# Get all combinations of PB/BM and 13/31/43/63
all_samples <- apply(expand.grid(c("BM", "PB"), c("13","31","43","63")), 1, paste, collapse="")

# Import GEX data and rename it
all_gex <- lapply(all_samples,
                  function(ID)  Read10X(data.dir =
                                          str_glue('Data/10x/GEX/{ID}/filtered_feature_bc_matrix'))) 
names(all_gex) <- all_samples


# 2.2 Load VDJ data --------------------------------------------------------------------------------

# Define a function to import the VDJ data
fetch_VDJ <- function(ID){
  read.csv(str_glue('Data/10x/GEX/{ID}/filtered_contig_annotations.csv'), 
           stringsAsFactors = FALSE) %>% # read in csv file
    combineTCR(samples = str_glue('Donor{str_sub(ID, 3, 4)}'), 
               ID = str_to_upper(str_sub(ID, 1, 2)), cells ="T-AB") %>% 
    # Use scRepertoire::combineTCR to consolidate TCR and GEX data. 
    .[[1]] %>% # change from a list to a data.frame
    mutate(barcode = str_split_fixed(barcode, "_", 3)[,3]) # fix up the barcode names
}

# Get all of the VDJ data in one list and rename it
all_vdj <- lapply(all_samples, fetch_VDJ)
names(all_vdj) <- all_samples


# 3. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_gex, "Data/R_out/10x/all_gex.rds")
saveRDS(all_vdj, "Data/R_out/10x/all_vdj.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '1. Load data' script")