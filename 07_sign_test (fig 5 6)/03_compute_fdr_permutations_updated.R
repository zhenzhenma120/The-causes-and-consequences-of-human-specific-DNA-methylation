#!/usr/bin/env Rscript
# Compute FDR via Permutation Testing - OPTIMIZED VERSION
# FDR = expected_false / observed_discoveries

suppressPackageStartupMessages({
  library(Matrix)
  library(parallel)
  library(dplyr)
})

source("01_enrichment_functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 4_compute_fdr_permutations_optimized.R <output_dir> <n_permutations>")
}
output_dir <- args[1]
n_permutations <- as.integer(args[2])


membership_matrix <- readRDS(file.path(output_dir, "04_membership_matrix.rds"))
sign_data <- readRDS(file.path(output_dir, "05_sign_vectors.rds"))
meth_results <- readRDS(file.path(output_dir, "06_meth_enrichment_results.rds"))
discord_results <- readRDS(file.path(output_dir, "07_discord_enrichment_results.rds"))
meth_signs <- sign_data$meth_signs
discord_signs <- sign_data$discord_signs
meth_results$gene_set_database <- get_pathway_database(meth_results$pathway)
discord_results$gene_set_database <- get_pathway_database(discord_results$pathway)
gene_set_summary <- table(meth_results$gene_set_database)

# Exclude GO terms from analysis
non_go_pathways <- meth_results$pathway[meth_results$gene_set_database != "GO"]
cat(sprintf("Pathways before GO exclusion: %d\n", nrow(meth_results)))
cat(sprintf("Pathways after GO exclusion: %d\n", length(non_go_pathways)))

meth_results <- meth_results %>% filter(gene_set_database != "GO")
discord_results <- discord_results %>% filter(gene_set_database != "GO")
membership_matrix <- membership_matrix[non_go_pathways, , drop = FALSE]

# Filter criteria:
# - More than 20 and less than 200 genes with methylation data
# - More than 10 and less than 200 genes with both meth and expr data (discord)

fdr_eligible_pathways <- filter_pathways_for_fdr(
  meth_results, discord_results,
  meth_min = 20, meth_max = 200,
  discord_min = 10, discord_max = 200
)

cat(sprintf("\nNum pathways for FDR: %d\n", length(fdr_eligible_pathways)))
meth_results$fdr_eligible <- meth_results$pathway %in% fdr_eligible_pathways
discord_results$fdr_eligible <- discord_results$pathway %in% fdr_eligible_pathways
saveRDS(fdr_eligible_pathways, file.path(output_dir, "08_fdr_eligible_pathways.rds"))

# ============================================================================
# run permutation
# ============================================================================
n_cores <- min(detectCores() - 1, 20)

membership_matrix_fdr <- membership_matrix[fdr_eligible_pathways, , drop = FALSE]
meth_null_results <- mclapply(1:n_permutations, function(i) {
  if (i %% 100 == 0) {
    cat(sprintf("  Completed %d/%d permutations\n", i, n_permutations))
  }

  permuted_signs <- permute_signs(meth_signs, NULL)
  perm_results <- vectorized_binomial_test(membership_matrix_fdr, permuted_signs)

  list(
    p_human = perm_results$p_human,
    p_chimp = perm_results$p_chimp
  )
}, mc.cores = n_cores)

meth_null_human <- do.call(cbind, lapply(meth_null_results, function(x) x$p_human))
meth_null_chimp <- do.call(cbind, lapply(meth_null_results, function(x) x$p_chimp))

rownames(meth_null_human) <- fdr_eligible_pathways
rownames(meth_null_chimp) <- fdr_eligible_pathways
saveRDS(meth_null_human, file.path(output_dir, "09_meth_null_human_optimized.rds"))
saveRDS(meth_null_chimp, file.path(output_dir, "09_meth_null_chimp_optimized.rds"))

# ============================================================================
# compute fdr
# ============================================================================
meth_results$fdr_human <- NA_real_
meth_results$fdr_chimp <- NA_real_

eligible_idx <- meth_results$fdr_eligible

if (sum(eligible_idx) > 0) {
  cat(sprintf("\nProcessing %d FDR-eligible pathways (KEGG+HPO combined)...\n", sum(eligible_idx)))

  # Get observed p-values for eligible pathways
  obs_pvals_human <- meth_results$p_human[eligible_idx]
  obs_pvals_chimp <- meth_results$p_chimp[eligible_idx]

  eligible_pathway_names <- meth_results$pathway[eligible_idx]

  local_idx <- 1
  for (i in which(eligible_idx)) {
    pathway_name <- meth_results$pathway[i]

    # --- FDR for HUMAN ---
    obs_p_human <- obs_pvals_human[local_idx]
    if (!is.na(obs_p_human)) {
      expected_false_human <- sum(meth_null_human <= obs_p_human, na.rm = TRUE) / n_permutations
      observed_discoveries_human <- sum(obs_pvals_human <= obs_p_human, na.rm = TRUE)
      if (observed_discoveries_human > 0) {
        meth_results$fdr_human[i] <- min(1, expected_false_human / observed_discoveries_human)
      }
    }

    # --- FDR for CHIMP ---
    obs_p_chimp <- obs_pvals_chimp[local_idx]
    if (!is.na(obs_p_chimp)) {
      expected_false_chimp <- sum(meth_null_chimp <= obs_p_chimp, na.rm = TRUE) / n_permutations
      observed_discoveries_chimp <- sum(obs_pvals_chimp <= obs_p_chimp, na.rm = TRUE)
      if (observed_discoveries_chimp > 0) {
        meth_results$fdr_chimp[i] <- min(1, expected_false_chimp / observed_discoveries_chimp)
      }
    }

    local_idx <- local_idx + 1
  }
}

cat(sprintf("\nMethylation enrichment FDR summary:\n"))
cat(sprintf("  Total pathways: %d\n", nrow(meth_results)))
cat(sprintf("  FDR-eligible pathways: %d\n", sum(meth_results$fdr_eligible)))
cat(sprintf("  HUMAN - FDR <= 0.05: %d | <= 0.10: %d | <= 0.20: %d\n",
            sum(meth_results$fdr_human <= 0.05, na.rm = TRUE),
            sum(meth_results$fdr_human <= 0.10, na.rm = TRUE),
            sum(meth_results$fdr_human <= 0.20, na.rm = TRUE)))
cat(sprintf("  CHIMP - FDR <= 0.05: %d | <= 0.10: %d | <= 0.20: %d\n",
            sum(meth_results$fdr_chimp <= 0.05, na.rm = TRUE),
            sum(meth_results$fdr_chimp <= 0.10, na.rm = TRUE),
            sum(meth_results$fdr_chimp <= 0.20, na.rm = TRUE)))

saveRDS(meth_results, file.path(output_dir, "10_meth_results_fdr_optimized.rds"))
write.table(meth_results,
            file.path(output_dir, "10_meth_results_fdr_optimized.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# run permutation for meth&expr sign test
# ============================================================================

discord_null_results <- mclapply(1:n_permutations, function(i) {
  if (i %% 100 == 0) {
    cat(sprintf("  Completed %d/%d permutations\n", i, n_permutations))
  }

  permuted_signs <- permute_signs(discord_signs, NULL)
  perm_results <- vectorized_binomial_test(membership_matrix_fdr, permuted_signs)

  list(
    p_human = perm_results$p_human,
    p_chimp = perm_results$p_chimp
  )
}, mc.cores = n_cores)

discord_null_human <- do.call(cbind, lapply(discord_null_results, function(x) x$p_human))
discord_null_chimp <- do.call(cbind, lapply(discord_null_results, function(x) x$p_chimp))

rownames(discord_null_human) <- fdr_eligible_pathways
rownames(discord_null_chimp) <- fdr_eligible_pathways

saveRDS(discord_null_human, file.path(output_dir, "11_discord_null_human_optimized.rds"))
saveRDS(discord_null_chimp, file.path(output_dir, "11_discord_null_chimp_optimized.rds"))

# ============================================================================
# compute fdr for meth&expr sign test
# ============================================================================

discord_results$fdr_human <- NA_real_
discord_results$fdr_chimp <- NA_real_

eligible_idx <- discord_results$fdr_eligible

if (sum(eligible_idx) > 0) {

  obs_pvals_human <- discord_results$p_human[eligible_idx]
  obs_pvals_chimp <- discord_results$p_chimp[eligible_idx]

  local_idx <- 1
  for (i in which(eligible_idx)) {

    # --- FDR for HUMAN ---
    obs_p_human <- obs_pvals_human[local_idx]
    if (!is.na(obs_p_human)) {
      expected_false_human <- sum(discord_null_human <= obs_p_human, na.rm = TRUE) / n_permutations
      observed_discoveries_human <- sum(obs_pvals_human <= obs_p_human, na.rm = TRUE)
      if (observed_discoveries_human > 0) {
        discord_results$fdr_human[i] <- min(1, expected_false_human / observed_discoveries_human)
      }
    }

    # --- FDR for CHIMP ---
    obs_p_chimp <- obs_pvals_chimp[local_idx]
    if (!is.na(obs_p_chimp)) {
      expected_false_chimp <- sum(discord_null_chimp <= obs_p_chimp, na.rm = TRUE) / n_permutations
      observed_discoveries_chimp <- sum(obs_pvals_chimp <= obs_p_chimp, na.rm = TRUE)
      if (observed_discoveries_chimp > 0) {
        discord_results$fdr_chimp[i] <- min(1, expected_false_chimp / observed_discoveries_chimp)
      }
    }

    local_idx <- local_idx + 1
  }

  cat(sprintf("  HUMAN FDR <= 0.05: %d | <= 0.10: %d | <= 0.20: %d\n",
              sum(discord_results$fdr_human <= 0.05, na.rm = TRUE),
              sum(discord_results$fdr_human <= 0.10, na.rm = TRUE),
              sum(discord_results$fdr_human <= 0.20, na.rm = TRUE)))
  cat(sprintf("  CHIMP FDR <= 0.05: %d | <= 0.10: %d | <= 0.20: %d\n",
              sum(discord_results$fdr_chimp <= 0.05, na.rm = TRUE),
              sum(discord_results$fdr_chimp <= 0.10, na.rm = TRUE),
              sum(discord_results$fdr_chimp <= 0.20, na.rm = TRUE)))
}

cat(sprintf("\nDiscordant enrichment FDR summary:\n"))
cat(sprintf("  Total pathways: %d\n", nrow(discord_results)))
cat(sprintf("  FDR-eligible pathways: %d\n", sum(discord_results$fdr_eligible)))
cat(sprintf("  HUMAN - FDR <= 0.05: %d | <= 0.10: %d | <= 0.20: %d\n",
            sum(discord_results$fdr_human <= 0.05, na.rm = TRUE),
            sum(discord_results$fdr_human <= 0.10, na.rm = TRUE),
            sum(discord_results$fdr_human <= 0.20, na.rm = TRUE)))
cat(sprintf("  CHIMP - FDR <= 0.05: %d | <= 0.10: %d | <= 0.20: %d\n",
            sum(discord_results$fdr_chimp <= 0.05, na.rm = TRUE),
            sum(discord_results$fdr_chimp <= 0.10, na.rm = TRUE),
            sum(discord_results$fdr_chimp <= 0.20, na.rm = TRUE)))

saveRDS(discord_results, file.path(output_dir, "12_discord_results_fdr_optimized.rds"))
write.table(discord_results,
            file.path(output_dir, "12_discord_results_fdr_optimized.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)