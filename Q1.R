# Question 1:
# The series GSE135635 RNA-sequencing in plasmacytoid dendritic cells from patients with primary Sjogren's syndrome (pSS), non-Sjogren's sicca (nSS), and healthy donors (HC). Use the GEO data file GSE135635_VST_normalized_data_disc.csv.gz to find differentially expressed genes between
# A. pSS and nSS
# B. pSS and HC
# C. nSS and HC


library(edgeR)
library(GEOquery)

GSE <- "GSE135635"
