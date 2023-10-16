# Script information -------------------------------------------------------------------------------

# Title: Data cleanup
# Author: James Favaloro
# Date: 2023-08-21
# Description: In this script we will perform some basic filtering on the Seurat objects to 
# eliminate poor quality cells and empty droplets. We will also create some plots and save them out 
# for later reference.


# 1. Import libraries ------------------------------------------------------------------------------

start_time <- Sys.time()
print("#># Start running '3. Data cleanup' script")

suppressPackageStartupMessages({
librarian::shelf(Seurat, tidyverse, scRepertoire, cowplot)
})


# 2. Load data, define functions and create character vectors --------------------------------------

# Load data
all_seurat <- readRDS("Data/R_out/10x/all_seurat.rds")
all_vdj <- readRDS("Data/R_out/10x/all_vdj.rds")

# Define a function to pull metadata columns from our list of Seurat objects
get_metadata_column <- function(input_seurat_list, column_name){
  # pull all of the column from the Seurat object
  all_cols <- lapply(input_seurat_list, function(x) x[[]][,column_name])
  return(data.frame(lapply(all_cols, "length<-", max(lengths(all_cols)))) %>% gather) 
  # convert list to data.frame, then from wide to long format
}


# 3. Sample overview -------------------------------------------------------------------------------

# Get some QC metrics for plotting
nCount_RNA <- get_metadata_column(input_seurat_list = all_seurat, column_name = "nCount_RNA")
nFeature_RNA <- get_metadata_column(input_seurat_list = all_seurat, column_name = "nFeature_RNA")
percent.mt <- get_metadata_column(input_seurat_list = all_seurat, column_name = "percent.mt")

# Define function to produce summary statistics (median and +/- 1, 1.5 & 2*sd) for plotting
for (i in c(1,1.5,2)){
  data_summary <- function(x) {
    m <- median(x)
    ymin <- m-(i*sd(x))
    ymax <- m+i*(sd(x))
    return(c(y=m,ymin=ymin,ymax=ymax))}
  
  p1 <- ggplot(nCount_RNA, aes(key, value)) + geom_violin() + ylim(NA, 1E4) +
    ggtitle("Number of UMIs detected") + xlab("Sample") + ylab("UMIs detected") +
    stat_summary(fun.data=data_summary, color = "red", geom = "pointrange", size =0.3)
  p2 <- ggplot(nFeature_RNA, aes(key, value)) + geom_violin() + ylim(NA, 3E3) +
    ggtitle("Number of genes detected") + xlab("Sample") + ylab("Genes detected") +
    stat_summary(fun.data=data_summary, color = "red", geom = "pointrange", size =0.3)
  p3 <- ggplot(percent.mt, aes(key, value)) + geom_violin() + ylim(NA, 20) +
    ggtitle("Mitochondrial genes (%)") + xlab("sample") + ylab("Mitochondrial gene (%)") +
    stat_summary(fun.data=data_summary, color = "red", geom = "pointrange", size =0.3)
  ggsave(plot = plot_grid(p1,p2,p3, ncol = 3), 
         filename = paste0("Output/QC/10x/Filter/Seurat_QC_metrics_overview_",
                           i,"SD", ".png"), width = 10, height = 3)
}


# 4. Filter data -----------------------------------------------------------------------------------

# 4.1 Filter Seurat --------------------------------------------------------------------------------

# NB: Analysis suggests +/- 1.5 SD for all QC metrics looks appropriate 
# Define function to filter Seurat objects based on the median +/- 1.5 SD for all three QC metrics
filter_seurat <- function(input_seurat){
  filtered_seurat <- subset(input_seurat, subset =
                              nFeature_RNA > median(input_seurat$nFeature_RNA) - 
                              (1.5*sd(input_seurat$nFeature_RNA)) &
                              nFeature_RNA < median(input_seurat$nFeature_RNA) + 
                              (1.5*sd(input_seurat$nFeature_RNA)) &
                              nCount_RNA > median(input_seurat$nCount_RNA) - 
                              (1.5*sd(input_seurat$nCount_RNA)) &
                              nCount_RNA < median(input_seurat$nCount_RNA) + 
                              (1.5*sd(input_seurat$nCount_RNA)) &
                              percent.mt > median(input_seurat$percent.mt) - 
                              (1.5*sd(input_seurat$percent.mt)) &
                              percent.mt < median(input_seurat$percent.mt) + 
                              (1.5*sd(input_seurat$percent.mt)))
  return(filtered_seurat)
}

# Apply filter
all_seurat_filtered <- lapply(all_seurat, filter_seurat)


# 4.2 Filter V(D)J ---------------------------------------------------------------------------------

# Define function to filter V(D)J based on missing information in both TCR1 and TCR2
# This will filter cells for which we have no V(D)J information
# NB: There is no additional pre-processing being performed here, thus, this data will also NOT have
# patient specific barcodes attached. It may or may not be used later...

filter_vdj <- function(input_vdj){
  filtered_vdj <- input_vdj[!with(input_vdj, is.na(input_vdj$TCR1) & is.na(input_vdj$TCR2)),]
  return(filtered_vdj)
}

# Apply filter
all_vdj_filtered <- lapply(all_vdj, filter_vdj)

## Look at the number of cells in the objects before and after filtering

# Before filtering Seurat objects
sapply(all_seurat, ncol)
# BM13   PB13  BM31  PB31  BM43  PB43  BM63  PB63 
# 10000  9985  3081  9632  8662  8040 10000  7512 

# After filtering Seurat objects
sapply(all_seurat_filtered, ncol)
# BM13 PB13 BM31 PB31 BM43 PB43 BM63 PB63 
# 7871 7822 2429 7558 7120 6497 7824 6275 

# Before filtering V(D)J
sapply(all_vdj, nrow)
# BM13 PB13 BM31 PB31 BM43 PB43 BM63 PB63 
# 8091 8351 6262 7233 7993 6909 9202 6747 

# After filtering V(D)J
sapply(all_vdj_filtered, nrow)
# BM13 PB13 BM31 PB31 BM43 PB43 BM63 PB63 
# 8016 8187 6161 7126 7835 6852 9041 6716 


# 5. Save out data and print stats -----------------------------------------------------------------

saveRDS(all_seurat_filtered, "Data/R_out/10x/all_seurat_filtered.rds")
saveRDS(all_vdj_filtered, "Data/R_out/10x/all_vdj_filtered.rds")

end_time <- Sys.time()
end_time - start_time
gc()

print("#># Finished running '3. Data cleanup' script")