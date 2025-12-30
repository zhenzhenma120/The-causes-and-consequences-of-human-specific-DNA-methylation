library(edgeR)
library(dplyr)
library(readr)

input_dir  <- "merged_samples" #merged raw counts from all hybrid samples (hybC and hybH)
output_dir <- "outputs"

target_celltypes <- c("DA", "CNCC", "HP", "SKM", "IPSC")

dir.create(file.path(output_dir, "tmm_normalized"),
           showWarnings = FALSE, recursive = TRUE)

extract_celltype <- function(key) strsplit(key, "_") |> sapply(`[`, 1)
counts <- read_csv(
  file.path(input_dir, "all_samples_merged.csv"),
  show_col_types = FALSE
) |>
  column_to_rownames("gene")

sample_info <- read_tsv(
  file.path(output_dir, "sample_info.txt"),
  show_col_types = FALSE
)

sample_info$celltype <- extract_celltype(sample_info$key)

sample_info_filt <- sample_info |>
  filter(type == "Hybrid",
         celltype %in% target_celltypes)

samples_to_use <- intersect(sample_info_filt$filename,
                            colnames(counts))

counts_filt <- counts[, samples_to_use, drop = FALSE]
counts_filt <- counts_filt[rowSums(counts_filt) > 0, , drop = FALSE]

dge <- DGEList(counts = round(counts_filt))
dge <- calcNormFactors(dge, method = "TMM")

cpm_tmm <- cpm(dge, log = FALSE)

write.csv(
  cpm_tmm,
  file = file.path(output_dir,
                   "tmm_normalized",
                   "hybrid_samples_TMM_CPM.csv")
)

cat(sprintf(
  "Saved TMM-normalized CPM: %d genes × %d samples\n",
  nrow(cpm_tmm), ncol(cpm_tmm)
))
