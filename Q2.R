# Question 2:
# The series GSE192829 is taken from a study aimed to explore the pathobiological markers of sarcoidosis in PBMCs by comparing the transcriptional signature of PBMCs from patients with pulmonary sarcoidosis and those of healthy controls by RNA sequencing. Use the GEO data file GSE192829_DK20193B_DK20141B_DK21066B_raw_counts.txt.gz to find differentially expressed between pulmonary sarcoidosis patients and healthy controls.


library(edgeR)
library(GEOquery)

GSE <- "GSE192829"
