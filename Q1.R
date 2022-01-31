# Question 1:
# The series GSE135635 RNA-sequencing in plasmacytoid dendritic cells from patients with primary Sjogren's syndrome (pSS), non-Sjogren's sicca (nSS), and healthy donors (HC). Use the GEO data file GSE135635_VST_normalized_data_disc.csv.gz to find differentially expressed genes between
# A. pSS and nSS
# B. pSS and HC
# C. nSS and HC


library(edgeR)
library(GEOquery)
library(tidyverse)
library(Homo.sapiens)

GSE <- "GSE135635"

experiment <- getGEO(GSE)

sample_data <-
  pData(experiment[["GSE135635-GPL18573_series_matrix.txt.gz"]]) |>
  dplyr::select(title, geo_accession, ends_with(":ch1")) |>
  rename_with(\(x) {
    str_to_lower(str_replace(str_remove(x, ":ch1"), " ", "_"))
  }) |>
  mutate(rowname = title) |>
  as_tibble() |>
  column_to_rownames()


getGEOSuppFiles(GSE, baseDir = "datafiles", filter_regex = "VST_normalized_data_disc")

count_data <-
  read_csv("datafiles/GSE135635/GSE135635_VST_normalized_data_disc.csv.gz") |>
  dplyr::rename(rowname = ...1) |>
  column_to_rownames()


gene_data <-
  select(
    Homo.sapiens,
    keys = rownames(count_data),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENENAME", "GENETYPE", "ENTREZID", "ENSEMBL"),
    multiVals = "first"
  ) |>
  as_tibble() |>
  group_by(ENSEMBL) |>
  slice_head(n = 1) |>
  ungroup() |>
  mutate(rowname = ENSEMBL) |>
  column_to_rownames()

dge <- DGEList(counts = count_data,
               samples = sample_data,
               genes = gene_data)

keep <- filterByExpr(dge)

dge_subset <- dge[keep, , keep.lib.sizes = FALSE]

dge_subset <- calcNormFactors(dge_subset)

design <- with(sample_data, model.matrix( ~ 0 + group))

dge_subset <- estimateDisp(dge_subset, design)

fit <- glmQLFit(dge_subset, design)

contrasts <- makeContrasts(
  A = grouppSS - groupnSS,
  B = grouppSS - groupHC,
  C = groupnSS - groupHC,
  levels = design
)

# Comparison A ------------------------------------------------------------

qlf.A <- glmQLFTest(fit, contrast = contrasts[, "A"])

diff_genes.A <- topTags(qlf.A, n = 10, sort.by = "logFC") |>
  pluck("table") |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  write_csv("results/GSE135635-A-DGE.csv")

# Comparison B ------------------------------------------------------------

qlf.B <- glmQLFTest(fit, contrast = contrasts[, "B"])

diff_genes.B <- topTags(qlf.B, n = 10, sort.by = "logFC") |>
  pluck("table") |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  write_csv("results/GSE135635-B-DGE.csv")

# Comparison C -------------------------------------------------------------

qlf.C <- glmQLFTest(fit, contrast = contrasts[, "C"])

diff_genes.C <- topTags(qlf.C, n = 10, sort.by = "logFC") |>
  pluck("table") |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  write_csv("results/GSE135635-C-DGE.csv")
