#!/bin/bash

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
#         File: pipeline.sh
#               Execute MM CD8+ T-cell project analysis
#       Author: James Favaloro
#        email: James.Favaloro@health.nsw.gov.au / Christian.Bryant@health.nsw.gov.au
#  Modified on: 2023/08/21
#      Version: 1.0
#      Example: ./pipeline.sh
#               
#
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

echo "This analysis was successfully run on a 2017 iMac, core i7 with 64GB of 2400MHz DDR4 ram. The authors extend their gratitude to Prof Joseph Powell and Dr Walter Muskovic for their assistance in the early stages of this project"

# 1. Load 10x 5' GEX and VDJ data - 62 lines, 2,462 characters
RScript Scripts/10x_Processing/1.Load_data.R

# 2. Create Seurat objects from 10x data - 153 lines, 6,754 characters
RScript Scripts/10x_Processing/2.Create_Seurat_objects.R

# 3. Perform QC of 10x data - 140 lines, 6,175 characters
RScript Scripts/10x_Processing/3.Data_cleanup.R

# 4. Perform further QC and SCTransform of 10x data - 201 lines, 7,904 characters
RScript Scripts/10x_Processing/4.Further_QC_and_SCTransform.R

# 5. Integrate and Cluster 10x data - 247 lines, 11,045 characters
RScript Scripts/10x_Processing/5.Integration_Clustering.R

# 6. Repeat Integration and Clustering 10x data - 233 lines, 9,613 characters
RScript Scripts/10x_Processing/6.Integration_Clustering_2.R

# 7. Perform Pre-processing of Controls dataset #1 - 139 lines, 6,327 characters
RScript Scripts/Controls_Processing/1.Controls_1_Preprocessing.R

# 8. Perform QC and SCTransform of Controls dataset #1 - 311 lines, 13,072 characters
RScript Scripts/Controls_Processing/2.Controls_1_QC_SCTransform.R

# 9. Perform Integration and Clustering of Controls dataset #1 - 184 lines, 7,626 characters
RScript Scripts/Controls_Processing/3.Controls_1_Integration_Clustering.R

# 10. Perform Pre-processing of Controls dataset #2 - 274 lines, 12,460 characters
RScript Scripts/Controls_Processing/4.Controls_2_Preprocessing.R

# 11. Perform QC and SCTransform of Controls dataset #2 - 254 lines, 10,147 characters
RScript Scripts/Controls_Processing/5.Controls_2_QC_SCTransform.R

# 12. Perform Integration and Clustering of Controls dataset #2 - 143 lines, 5,435 characters
RScript Scripts/Controls_Processing/6.Controls_2_Integration_Clustering.R

# 13. Perform complete workflow for Controls dataset #3 - 163 lines, 9,299 characters
RScript Scripts/Controls_Processing/7.Controls_3_Complete_workflow.R

# 14. Perform DE Workup of 10x data - 1,027 lines, 45,403 characters
RScript Scripts/Analysis/DE/1.DE_Workup_10x.R

# 15. Perform DE Analysis of 10x data - 620 lines, 30,943 characters
RScript Scripts/Analysis/DE/2.DE_Analysis_10x.R

# 16. Perform DE Workup of Controls dataset #1 - 1,515 lines, 68,804 characters
RScript Scripts/Analysis/DE/3.DE_Workup_Controls_1.R

# 17. Perform DE Analysis of Controls dataset #1 - 477 lines, 21,479 characters
RScript Scripts/Analysis/DE/4.DE_Analysis_Controls_1.R

# 18. Perform DE Workup of Controls dataset #2 - 733 lines, 32,535 characters
RScript Scripts/Analysis/DE/5.DE_Workup_Controls_2.R

# 19. Perform DE Analysis of Controls dataset #2 - 353 lines, 16,407 characters
RScript Scripts/Analysis/DE/6.DE_Analysis_Controls_2.R

# 20. Perfrom DE Workup and Analysis of Controls dataset #3 - 363 lines, 15,939 characters
RScript Scripts/Analysis/DE/7.DE_Workup_and_Analysis_Controls_3.R

# 21. Annotate 10x data using ProjecTILs and SingleR - 148 lines, 6,300 characters
RScript Scripts/Analysis/ProjecTILs/1.Annotate_10x_data.R

# 22. Analysis of ProjecTILs projections and SingleR annotations - 171 lines, 7,595 characters
RScript Scripts/Analysis/ProjecTILs/2.10x_annotation_analysis.R

# 23. Create reference atlases using 10x data - 140 lines, 5,277 characters
RScript Scripts/Analysis/ProjecTILs/3.Create_10x_reference_atlases.R

# 24. Perform DE Workup of the 10x reference atlases - 1,012 lines, 43,862 characters
RScript Scripts/Analysis/ProjecTILs/4.10x_reference_atlases_DE_workup.R

# 25. Perform DE Analysis of the 10x reference atlases - 570 lines, 28,968 characters
RScript Scripts/Analysis/ProjecTILs/5.10x_reference_atlases_analysis.R

# 26. Project Controls dataset onto custom 10x reference atlases - 159 lines, 8,040 characters
RScript Scripts/Analysis/ProjecTILs/6.Project_Controls_data_on_10x_atlases.R

# 27. Perform DE and analysis of projected Controls data - 924 lines, 45,782 characters
RScript Scripts/Analysis/ProjecTILs/7.Controls_vs._Reference_DE_testing.R

# 28. Output additional figures for Controls vs. Reference data - 742 lines, 36,484 characters
RScript Scripts/Analysis/ProjecTILs/8.Controls_projection_and_figures.R

# 29. Run scGSEA - 121 lines, 4,752 characters
Rscript Scripts/Analysis/GSEA/1.scGSEA_escape.R

# 30. Perform scGSEA analysis of 10x data - 157 lines, 6,625 characters
Rscript Scripts/Analysis/GSEA/2.scGSEA_escape_10x_Analysis.R

# 31. Perform scGSEA analysis of control data - 408 lines, 17,783 characters
Rscript Scripts/Analysis/GSEA/3.scGSEA_escape_Controls_Analysis.R

# 32. Run ReactomeGSA - 217 lines, 9,258 characters
Rscript Scripts/Analysis/GSEA/4.ReactomeGSA.R

# 33. Perform analysis of Reactome data - 255 lines, 10,805 characters
Rscript Scripts/Analysis/GSEA/5.ReactomeGSA_Analysis.R

# 34. Perform general V(D)J analysis - 410 lines, 19,920 characters
Rscript Scripts/Analysis/VDJ/1.General_Analysis.R

# 35. Perform comparative analysis - 213 lines, 9,356 characters
Rscript Scripts/Analysis/VDJ/2.Comparative_Analysis.R

# 36. Run ERGO-II webtool - 325 lines, 14,515 characters
#Rscript Scripts/Analysis/VDJ/3.ERGO_II_Workflow.R

# 36. Perform ERGO-II analysis - 252 lines, 12,148 characters
#Rscript Scripts/Analysis/VDJ/4.ERGO_II_Analysis.R

# 37. Perform Dominant Clonotype analysis - 409 lines, 17,456 characters
Rscript Scripts/Analysis/VDJ/5.Dominant_clonotypes_Analysis.R