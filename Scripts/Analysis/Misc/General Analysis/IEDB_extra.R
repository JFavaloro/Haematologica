#Load McPAS-TCR data
#McPas_TCR = dbLoad("FILE_NAME.csv", .db = "mcpas") // where FILE_NAME is saved *.csv file from http://friedmanlab.weizmann.ac.il/McPAS-TCR/ & .db is database type. Loads MCPAS_TCR data as dataframe "McPas_TCR". NB: Prior to loading this file open downloaded *.csv in excel and delete mouse and CD4 data. Save as McPas_TCR_filtered. The dataset is paired and as such CDR3 is seperated into alpha and beta columns.
# It's possible to compare alpha/beta + VJ with MCpas by editing the data format to match vdjdb and loading it like that
# McPAS_TCR_A = dbLoad("Data/Public datasets/20201030_McPAS_TCR_filtered_modified_TRA.tsv", .db = "vdjdb")
#a = dbAnnotate(immdata$data, McPAS_TCR_A, c("CDR3.aa", "V.name", "J.name"), c("CDR3.aa", "V.Name", "J.Name"))

#IEDB
# Have removed all superfluous information and renamed columns (calculated) to aa, V.name and J.name for both alpha and beta
#Still need to figure out the rename command in R to remove anything *0(n) to get V.name and J.name into 'simple' nomenclature below should work but keeps *
iedb = read.csv("Data/Public datasets/20201106_IEDB trimmed.csv")
iedb = lapply(iedb, gsub, pattern = "\\*[0-9]+", replacement = "") # Doesn't seem to work?
#Regular expressions LEARN
#Add flank amino acids back
iedb$beta.CDR3.aa = paste0("C", iedb$beta.CDR3.aa, "F")
iedb$alpha.CDR3.aa = paste0("C", iedb$alpha.CDR3.aa, "F")

iedb_clean <- iedb %>%
  mutate(across(.cols = c(beta.J.Name, beta.V.Name, alpha.V.Name, alpha.J.Name),
                ~gsub(pattern = "\\*[0-9]+", replacement = "", x = .x)),
         across(.cols = c(alpha.CDR3.aa, beta.CDR3.aa),
                ~if_else(.x == "",
                         true = "",
                         false = paste0("C", .x, "F"))))



function(x) {
  gsub(pattern = "\\*[0-9]+", replacement = "", x)
}

strip_after("ABCDEF*2029")

# Remove clonotypes missing CDR3.aa information and split IEDB into TRA and TRB, then remove missing V.name and J.name for TRB (alpha looks to be complete data)
#This doesn't seem to work correctly - is there really only 400 "complete" TRA?
#I think this dataset may be inappropriate for this sort of analysis... Leaving this for now.
iedb_tra = iedb_clean[which(str_detect(iedb_clean$alpha.CDR3.aa, "")),]
iedb_trb = iedb_clean[which(str_detect(iedb_clean$beta.CDR3.aa, "")),]
iedb_trb = iedb_trb[which(str_detect(iedb_clean$beta.V.Name, "")),]
iedb_trb = iedb_trb[which(str_detect(iedb_clean$beta.J.Name, "")),]
#Write out the files as CSV and rearrange the columns as per the other datasets, then reimport them
write.csv(iedb_tra, "iedb_tra.csv")
write.csv(iedb_trb, "iedb_trb.csv")
iedb_tra = read.csv("Data/Public datasets/iedb_tra.csv")
iedb_trb = read.csv("Data/Public datasets/iedb_trb.csv")
# Junk
#Try clonotype overlap between single cell data for top 10,20,50,100,6262 cells...
# This won't work... Need to do this from # of clonotypes...
# Downsample 10x data to lowest number of observed clonotypes (1378) and repeat...
#Again, better to write a loop here but manual is fine for now...
# Overlap can be compared between what is decided as a clonotype
# CT Strict is TRUE definition, so let's try that first - can't do this with immdata...
# Do overlap on CDR3.aa instead
#Will save this code to file; 
sample1 = head(immdata$data$BM013_TRA, n = 1378)
sample2 = head(immdata$data$PB013_TRA, n = 1378)
x = intersect(sample1$CDR3.aa, sample2$CDR3.aa)
#Why is this not getting any hits? Venn diagrams work...
# Because I'm an idiot and got the names wrong in the list... Try again. Works. Barely. Huzzah.


circ = list(bm13 = immdata$data$BM013_filtered_contig_annotations, bm31 = immdata$data$BM031_filtered_contig_annotations, bm43 = immdata$data$BM043_filtered_contig_annotations, bm63 = immdata$data$BM063_filtered_contig_annotations, pb13 = immdata$data$PB013_filtered_contig_annotations, pb31 = immdata$data$PB031_filtered_contig_annotations, pb43 = immdata$data$PB043_filtered_contig_annotations, pb63 = immdata$data$PB063_filtered_contig_annotations)


# VDJ analysis script
# This is ALL the VDJ data we have as of end of 2020
# Data set currently comprises 36 files
#PT13, PT31, PT43 & PT63 10xVDJ as both single cell (PTID_filtered_contig_annotations) and pseudobulk (PTID_TRX)
#PT31, PT43 & PT63 as bulk with single read clonotypes removed (DP dataset) (PTID.Clonotypes.TRX)
# Meta information for immunarch contains Sample (Sample name), Patient (Patient ID: 13, 31, 43 or 63), Tissue (bm or pb), Chain (a, b or ab), Disease (MM) & Type (10x, 10x_bulk or bulk)
#20201116 It was discovered that without manual removal of unproductive clonotypes, Immunarch will consider these clones and results in inclusion of extra data that may not be of interest/ may complicate analysis - e.g. "creation" of a TCRaab in PB043. This has been fixed and the 10x data in the immunarch format is now saved in the "Data/Clean" folder with the suffix "_manual_removal_of_non_productive"

#10x data is complicated by the fact that some cells are missing information on either a or b chain, or contain more than one a or b chain. Literature suggests that cells containing more than one a chain may be real. Due to allelic silencing of the b chain any cell that contains more than one b chain can be considered a doublet and may be removed prior to analysis if warranted.
#Thank you, Walter
#vdj <- all_vdj[[1]]
#vdj2 <- vdj[-which(str_detect(vdj$cdr3_aa2, ";")),]


