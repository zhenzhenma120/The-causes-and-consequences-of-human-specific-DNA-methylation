#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(reticulate)
  library(dplyr)
  library(parallel)
})

# Source helper functions
source("01_enrichment_functions.R")
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript main_enrichment_pipeline.R <input_matrix.txt> <output_dir> [meth_p_thresh] [expr_p_thresh]")
}

input_file <- args[1]
output_dir <- args[2]
meth_p_threshold <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.05)
expr_p_threshold <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


cat(sprintf("[%s] Loading pathway annotations...\n", Sys.time()))
pickle_path <- "/oak/stanford/groups/hbfraser/Zhenzhen/backup/objs.pkl"
pathways <- load_pathways_from_pickle(pickle_path)
pathways <- pathways[[1]] #the second list is for promoters
cat(sprintf("Loaded %d pathway databases:\n", length(pathways)))
# for (db_name in names(pathways)) {
#   cat(sprintf("  - %s: %d genes\n", db_name, length(pathways[[db_name]])))
# }
saveRDS(pathways, file.path(output_dir, "01_raw_pathways.rds"))


cat(sprintf("\n[%s] Loading input data...\n", Sys.time()))
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- data[!duplicated(data[,1]), ]
rownames(data) <- data[,1]
data <- data[,-1]
cat(sprintf("Loaded data: %d genes x %d columns\n", nrow(data), ncol(data)))
cat(sprintf("Columns: %s\n", paste(colnames(data), collapse = ", ")))

cat(sprintf("\nInput data summary:\n"))
cat(sprintf("  - Genes with methylation data: %d\n", sum(!is.na(data$methylation))))
cat(sprintf("  - Genes with expression data: %d\n", sum(!is.na(data$expression))))
cat(sprintf("  - Genes with meth p < %.4f: %d\n", 
            meth_p_threshold, sum(data$methylation_p <= meth_p_threshold, na.rm = TRUE)))
cat(sprintf("  - Genes with expr p < %.4f: %d\n", 
            expr_p_threshold, sum(data$expression_p <= expr_p_threshold, na.rm = TRUE)))
saveRDS(data, file.path(output_dir, "02_input_data.rds"))

# ============================================================================
# filter pathways
# ============================================================================
cat(sprintf("\n[%s] Filtering pathways...\n", Sys.time()))

all_genes <- rownames(data)
min_pathway_size <- 20
max_pathway_size <- 200
sfari_all <- pathways$SFARI_all
filtered_pathways <- filter_pathways(
  pathways, 
  all_genes, 
  min_size = min_pathway_size,
  max_size = max_pathway_size
)

filtered_pathways$SFARI_all <- sfari_all # add sfari all genes $ sfari score == 1
saveRDS(filtered_pathways, file.path(output_dir, "03_filtered_pathways.rds"))

# ============================================================================
# build membership matrix and create sign matrices
# ============================================================================
membership_matrix <- build_sparse_membership(filtered_pathways, all_genes)
saveRDS(membership_matrix, file.path(output_dir, "04_membership_matrix.rds"))

# Create thresholded sign vectors
# NA = gene not in input data
# 0 = gene present but not significant (p >= threshold)
# 1/-1 = significant human/chimp bias (p < threshold)
meth_signs <- create_thresholded_signs(
  log2fc = data$methylation,
  pvals = data$methylation_p,
  p_threshold = meth_p_threshold,
  gene_names = all_genes
)
expr_signs <- create_thresholded_signs(
  log2fc = data$expression,
  pvals = data$expression_p,
  p_threshold = expr_p_threshold,
  gene_names = all_genes
)
discord_signs <- create_discordant_signs(meth_signs, expr_signs)

sign_data <- list(
  meth_signs = meth_signs,
  expr_signs = expr_signs,
  discord_signs = discord_signs,
  meth_p_threshold = meth_p_threshold,
  expr_p_threshold = expr_p_threshold
)

sign_matrices <- create_sign_matrices(
  membership_matrix, 
  meth_signs, 
  expr_signs, 
  discord_signs
)
saveRDS(sign_data, file.path(output_dir, "05_sign_vectors.rds"))
saveRDS(sign_matrices, file.path(output_dir, "05b_sign_matrices.rds"))

# ============================================================================
# test 1: methylation-only sign test
# ============================================================================
cat(sprintf("\n[%s] Running methylation-only enrichment test...\n", Sys.time()))

meth_results <- vectorized_binomial_test(
  membership_matrix,
  meth_signs,
  test_type = "greater"
)
cat("Top 20 meth results human biased\n")
meth_results %>% filter(!grepl("GO:", pathway)) %>% filter(n_genes_tested > 20) %>%filter(!grepl("Homo sapiens", pathway)) %>%arrange(p_human) %>% head(20) %>% print()
cat("\n\nTop 20 meth results chimp biased\n")
meth_results %>% filter(!grepl("GO:", pathway)) %>% filter(n_genes_tested > 20) %>%filter(!grepl("Homo sapiens", pathway)) %>%arrange(p_chimp) %>% head(20) %>% print()
cat(sprintf("Significant human-biased (p < 0.05): %d pathways\n", 
            sum(meth_results$p_human <= 0.05, na.rm = TRUE)))
cat(sprintf("Significant chimp-biased (p < 0.05): %d pathways\n", 
            sum(meth_results$p_chimp <= 0.05, na.rm = TRUE)))

saveRDS(meth_results, file.path(output_dir, "06_meth_enrichment_results.rds"))
write.table(meth_results, 
            file.path(output_dir, "06_meth_enrichment_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# test 2: repressive methylation (meth&expr) sign test
# ============================================================================
discord_results <- vectorized_binomial_test(
  membership_matrix,
  discord_signs,
  test_type = "greater"
)
cat("Top 20 meth results human biased\n")
discord_results %>% filter(!grepl("GO:", pathway))%>% filter(n_genes_tested > 20) %>%filter(!grepl("Homo sapiens", pathway)) %>%arrange(p_human) %>% head(50) %>% print()
cat("\n\nTop 20 meth results chimp biased\n")
discord_results %>% filter(!grepl("GO:", pathway))%>% filter(n_genes_tested > 20) %>%filter(!grepl("Homo sapiens", pathway)) %>%arrange(p_chimp) %>% head(20) %>% print()
cat(sprintf("Completed %d pathway tests\n", nrow(discord_results)))
cat(sprintf("Significant human-biased (p < 0.05): %d pathways\n", 
            sum(discord_results$p_human <= 0.05, na.rm = TRUE)))
cat(sprintf("Significant chimp-biased (p < 0.05): %d pathways\n", 
            sum(discord_results$p_chimp <= 0.05, na.rm = TRUE)))

saveRDS(discord_results, file.path(output_dir, "07_discord_enrichment_results.rds"))
write.table(discord_results, 
            file.path(output_dir, "07_discord_enrichment_results.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# summary stats
# ============================================================================
summary_stats <- list(
  n_genes_total = length(all_genes),
  n_genes_with_meth_data = sum(!is.na(data$meth_log2fc)),
  n_genes_with_expr_data = sum(!is.na(data$expr_log2fc)),
  meth_p_threshold = meth_p_threshold,
  expr_p_threshold = expr_p_threshold,
  n_genes_meth_human = sum(meth_signs == 1, na.rm = TRUE),
  n_genes_meth_chimp = sum(meth_signs == -1, na.rm = TRUE),
  n_genes_meth_nonsig = sum(meth_signs == 0, na.rm = TRUE),
  n_genes_expr_human = sum(expr_signs == 1, na.rm = TRUE),
  n_genes_expr_chimp = sum(expr_signs == -1, na.rm = TRUE),
  n_genes_expr_nonsig = sum(expr_signs == 0, na.rm = TRUE),
  n_genes_discord_human = sum(discord_signs == 1, na.rm = TRUE),
  n_genes_discord_chimp = sum(discord_signs == -1, na.rm = TRUE),
  n_genes_discord_nonsig = sum(discord_signs == 0, na.rm = TRUE),
  n_pathways_tested = nrow(meth_results),
  n_databases = length(filtered_pathways),
  pathway_size_range = c(min_pathway_size, max_pathway_size),
  meth_results_summary = list(
    human_sig_p05 = sum(meth_results$p_human <= 0.05, na.rm = TRUE),
    human_sig_p01 = sum(meth_results$p_human <= 0.01, na.rm = TRUE),
    chimp_sig_p05 = sum(meth_results$p_chimp <= 0.05, na.rm = TRUE),
    chimp_sig_p01 = sum(meth_results$p_chimp <= 0.01, na.rm = TRUE)
  ),
  discord_results_summary = list(
    human_sig_p05 = sum(discord_results$p_human <= 0.05, na.rm = TRUE),
    human_sig_p01 = sum(discord_results$p_human <= 0.01, na.rm = TRUE),
    chimp_sig_p05 = sum(discord_results$p_chimp <= 0.05, na.rm = TRUE),
    chimp_sig_p01 = sum(discord_results$p_chimp <= 0.01, na.rm = TRUE)
  )
)

saveRDS(summary_stats, file.path(output_dir, "08_summary_statistics.rds"))

save(
  pathways,
  filtered_pathways,
  data,
  membership_matrix,
  meth_signs,
  expr_signs,
  discord_signs,
  sign_matrices,
  meth_results,
  discord_results,
  summary_stats,
  file = file.path(output_dir, "all_results.RData")
)