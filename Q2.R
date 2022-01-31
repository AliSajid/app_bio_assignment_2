# Question 2:
# The series GSE192829 is taken from a study aimed to explore the pathobiological markers of sarcoidosis in PBMCs by comparing the transcriptional signature of PBMCs from patients with pulmonary sarcoidosis and those of healthy controls by RNA sequencing. Use the GEO data file GSE192829_DK20193B_DK20141B_DK21066B_raw_counts.txt.gz to find differentially expressed between pulmonary sarcoidosis patients and healthy controls.


library(edgeR)
library(GEOquery)
library(tidyverse)

GSE <- "GSE192829"

experiment <- getGEO(GSE)

sample_data <- pData(experiment[["GSE192829_series_matrix.txt.gz"]]) |>
  select(title, geo_accession, ends_with(":ch1")) |>
  rename_with(\(x) {str_to_lower(str_replace(str_remove(x, ":ch1"), " ", "_"))}) |>
  mutate(group = if_else(diagnosis == "Pulmonary sarcoidosis", "CASE", "CTRL"),
         rowname = title) |>
  as_tibble() |>
  column_to_rownames()

getGEOSuppFiles(GSE, baseDir = "datafiles", filter_regex = "raw")

count_data <- read_tsv("datafiles/GSE192829/GSE192829_DK20193B_DK20141B_DK21066B_raw_counts.txt.gz") |>
  rename(rowname = ...1) |>
  column_to_rownames()


dge <- DGEList(counts = count_data, samples = sample_data,
               genes = rownames(count_data))

keep <- filterByExpr(dge)

dge_subset <- dge[keep,,keep.lib.sizes=FALSE]


dge_subset <- calcNormFactors(dge_subset)

design <- with(sample_data, model.matrix(~ 0 + group))

dge_subset <- estimateDisp(dge_subset, design)

fit <- glmQLFit(dge_subset, design)

qlf <- glmQLFTest(fit, contrast = c(1, -1))

diff_genes <- topTags(qlf, n = 20, sort.by = "logFC") |>
  pluck("table") |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  write_csv("results/GSE192829-DGE.csv")
