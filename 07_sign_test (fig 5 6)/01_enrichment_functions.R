# ============================================================================
# Load and filter pathways
# ============================================================================
load_pathways_from_pickle <- function(pickle_path) {
  tryCatch({
    # Load Python modules
    py <- import_builtins()
    pickle <- import("pickle")
    py_require("pandas")
    
    # Open and load pickle file
    pathway_dict <- py$open(pickle_path, "rb")
    pathways <- pickle$load(pathway_dict)
    return(pathways)
  }, error = function(e) {
    stop(sprintf("Error loading pickle file: %s\n%s", pickle_path, e$message))
  })
}

filter_pathways <- function(pathway_list, gene_universe, 
                            min_size = 5, max_size = 500) {
  filtered <- lapply(pathway_list, function(pathway_genes) {
    filtered_genes <- intersect(pathway_genes, gene_universe)
    pathway_size <- length(filtered_genes)
    if (pathway_size >= min_size & pathway_size <= max_size) {
      return(filtered_genes)
    } else {
      return(NULL)
    }
  })
  filtered <- Filter(Negate(is.null), filtered)
  return(filtered)
}

# ============================================================================
# Build membership matrix
# ============================================================================
build_sparse_membership <- function(pathway_list, gene_universe) {
  pathway_names <- names(pathway_list)
  gene_to_idx <- setNames(seq_along(gene_universe), gene_universe)
  triplet_list <- list()
  
  for (i in seq_along(pathway_list)) {
    genes <- unlist(pathway_list[[i]])  # Added unlist here
    gene_indices <- gene_to_idx[genes]
    gene_indices <- gene_indices[!is.na(gene_indices)]
    
    if (length(gene_indices) > 0) {
      triplet_list[[i]] <- data.frame(
        pathway = rep(i, length(gene_indices)),
        gene = gene_indices,
        value = rep(1, length(gene_indices))
      )
    }
  }
  
  triplet_list <- triplet_list[!sapply(triplet_list, is.null)]
  
  if (length(triplet_list) == 0) {
    stop("No genes matched between pathways and gene_universe")
  }
  
  triplets <- do.call(rbind, triplet_list)
  
  sparse_mat <- sparseMatrix(
    i = triplets$pathway,
    j = triplets$gene,
    x = triplets$value,
    dims = c(length(pathway_list), length(gene_universe)),
    dimnames = list(pathway_names, gene_universe)
  )
  
  return(sparse_mat)
}

# ============================================================================
# Create sign vectors
# ============================================================================
create_thresholded_signs <- function(log2fc, pvals, p_threshold, gene_names) {
  # Initialize with NA (missing data)
  signs <- rep(NA_real_, length(gene_names))
  names(signs) <- gene_names
  
  # For genes with data, apply p-value threshold
  has_data <- !is.na(log2fc) & !is.na(pvals)
  
  # Significant genes get their sign
  sig_genes <- has_data & (pvals <= p_threshold)
  signs[sig_genes] <- sign(log2fc[sig_genes])
  
  # Non-significant genes with data get 0
  nonsig_genes <- has_data & (pvals > p_threshold)
  signs[nonsig_genes] <- 0
  
  return(signs)
}

# ============================================================================
# Create discordant sign vector (opposite meth/expr, both pass threshold)
# ============================================================================
create_discordant_signs <- function(meth_signs, expr_signs) {
  # Initialize with NA
  discord <- rep(NA_real_, length(meth_signs))
  names(discord) <- names(meth_signs)
  
  # Only consider genes where BOTH meth and expr are significant (non-zero, non-NA)
  both_sig <- !is.na(meth_signs) & !is.na(expr_signs) & 
    (meth_signs != 0) & (expr_signs != 0)
  
  meth_vals <- meth_signs[both_sig]
  expr_vals <- expr_signs[both_sig]
  sig_idx <- which(both_sig)
  
  # Discordant: opposite signs (use meth sign for direction)
  is_discordant <- meth_vals * expr_vals < 0
  discord[sig_idx[is_discordant]] <- meth_vals[is_discordant]  # 1 or -1
  
  # Concordant: same signs (assign 0 to distinguish from discordant)
  is_concordant <- meth_vals * expr_vals > 0
  discord[sig_idx[is_concordant]] <- 0
  
  # Everything else stays NA
  return(discord)
}

create_sign_matrices <- function(membership_matrix, meth_signs, expr_signs, discord_signs) {
  # Get dimensions
  n_pathways <- nrow(membership_matrix)
  n_genes <- ncol(membership_matrix)
  pathway_names <- rownames(membership_matrix)
  gene_names <- colnames(membership_matrix)
  
  # Create sparse matrices by element-wise multiplication
  # membership_matrix has 1 where gene is in pathway, 0 otherwise
  # We want: NA where gene not in pathway, sign value where gene in pathway
  
  # For each sign vector, create a diagonal matrix and multiply
  meth_matrix <- create_pathway_sign_matrix(membership_matrix, meth_signs)
  expr_matrix <- create_pathway_sign_matrix(membership_matrix, expr_signs)
  discord_matrix <- create_pathway_sign_matrix(membership_matrix, discord_signs)
  
  return(list(
    meth_matrix = meth_matrix,
    expr_matrix = expr_matrix,
    discord_matrix = discord_matrix
  ))
}

create_pathway_sign_matrix <- function(membership_matrix, gene_signs) {
  # For each pathway, assign gene signs only to genes in that pathway
  # Genes not in pathway become NA
  
  n_pathways <- nrow(membership_matrix)
  n_genes <- ncol(membership_matrix)
  
  sign_matrix <- matrix(gene_signs, nrow = n_pathways, ncol = n_genes, byrow = TRUE)
  rownames(sign_matrix) <- rownames(membership_matrix)
  colnames(sign_matrix) <- colnames(membership_matrix)
  
  membership_dense <- as.matrix(membership_matrix)
  
  sign_matrix[membership_dense == 0] <- NA
  
  return(sign_matrix)
}

# ============================================================================
# Binomial sign test 
# ============================================================================
vectorized_binomial_test <- function(membership_mat, gene_signs, 
                                     test_type = "greater") {

  valid_idx <- !is.na(gene_signs) & (gene_signs != 0)
  signs_filtered <- gene_signs[valid_idx]
  membership_filtered <- membership_mat[, valid_idx, drop = FALSE]
  
  pathway_sizes <- Matrix::rowSums(membership_filtered)
  
  has_data_idx <- !is.na(gene_signs)
  membership_with_data <- membership_mat[, has_data_idx, drop = FALSE]
  pathway_sizes_with_data <- Matrix::rowSums(membership_with_data)
  
  human_idx <- signs_filtered == 1
  n_human_genes <- sum(human_idx)
  
  if (n_human_genes > 0) {
    n_human <- as.vector(membership_filtered[, human_idx, drop = FALSE] %*% 
                           rep(1, n_human_genes))
  } else {
    n_human <- rep(0, nrow(membership_filtered))
  }
  
  chimp_idx <- signs_filtered == -1
  n_chimp_genes <- sum(chimp_idx)
  
  if (n_chimp_genes > 0) {
    n_chimp <- as.vector(membership_filtered[, chimp_idx, drop = FALSE] %*% 
                           rep(1, n_chimp_genes))
  } else {
    n_chimp <- rep(0, nrow(membership_filtered))
  }
  
  total_genes <- length(signs_filtered)
  p_human_bg <- n_human_genes / total_genes
  p_chimp_bg <- n_chimp_genes / total_genes
  
  # P(X >= k) = 1 - P(X <= k-1) = pbinom(k-1, n, p, lower.tail=FALSE)
  p_human <- pbinom(n_human - 1, pathway_sizes, p_human_bg, lower.tail = FALSE)
  p_chimp <- pbinom(n_chimp - 1, pathway_sizes, p_chimp_bg, lower.tail = FALSE)
  p_human[pathway_sizes == 0] <- 1
  p_chimp[pathway_sizes == 0] <- 1
  
  results <- data.frame(
    pathway = rownames(membership_mat),
    n_genes_in_pathway = as.vector(Matrix::rowSums(membership_mat)),  # total in pathway
    n_genes_with_data = pathway_sizes_with_data,  # genes with meth/expr data
    n_genes_tested = pathway_sizes,  # genes with significant changes (non-zero)
    n_human = n_human,
    n_chimp = n_chimp,
    p_human = p_human,
    p_chimp = p_chimp,
    stringsAsFactors = FALSE
  )
  
  results$prop_human <- results$n_human / results$n_genes_tested
  results$prop_chimp <- results$n_chimp / results$n_genes_tested
  results$bg_prop_human <- p_human_bg
  results$bg_prop_chimp <- p_chimp_bg
  results$fe_human <- results$prop_human / p_human_bg
  results$fe_chimp <- results$prop_chimp / p_chimp_bg
  
  results$prop_human[results$n_genes_tested == 0] <- NA
  results$prop_chimp[results$n_genes_tested == 0] <- NA
  results$fe_human[results$n_genes_tested == 0] <- NA
  results$fe_chimp[results$n_genes_tested == 0] <- NA
  
  return(results)
}

# ============================================================================
# Permute gene signs
# ============================================================================
permute_signs <- function(gene_signs, permute_within = NULL) {
  permuted <- gene_signs
  
  # Only permute non-zero, non-NA values (keep 0s in place)
  nonzero_idx <- !is.na(gene_signs) & (gene_signs != 0)
  #nonzero_idx <- !is.na(gene_signs)
  if (is.null(permute_within)) {
    permuted[nonzero_idx] <- sample(gene_signs[nonzero_idx])
  } else {
    for (group in unique(permute_within[nonzero_idx])) {
      group_idx <- nonzero_idx & (permute_within == group)
      permuted[group_idx] <- sample(gene_signs[group_idx])
    }
  }
  
  names(permuted) <- names(gene_signs)
  return(permuted)
}

# ============================================================================
# Extract database name from pathway identifier (update: excludes GO terms)
# ============================================================================
extract_go_id <- function(pathway_name) {
  # Extract GO:xxxxxxx from pathway name
  # Example: "cell cycle (GO:0007049)" -> "GO:0007049"
  
  go_match <- regexpr("GO:\\d{7}", pathway_name)
  
  if (go_match > 0) {
    go_id <- regmatches(pathway_name, go_match)
    return(go_id)
  } else {
    return(NA_character_)
  }
}

get_pathway_database <- function(pathway_names, go_lookup = NULL, 
                                 categorize_go = FALSE) {
  
  categories <- sapply(pathway_names, function(pathway) {
    if (is.na(pathway) || pathway == "") return(NA_character_)
    
    if (grepl("\\(GO:", pathway)) {
      return("GO")  # Will be excluded
    }
    
    if (grepl(" - Homo sapiens \\(human\\)", pathway)) return("KEGG")
    
    return("HPO")
  }, USE.NAMES = FALSE)
  
  cat("done\n")
  return(categories)
}

# ============================================================================
# Filter pathways for FDR calculation based on coverage
# ============================================================================
filter_pathways_for_fdr <- function(meth_results, discord_results,
                                    meth_min = 20, meth_max = 200,
                                    discord_min = 10, discord_max = 200) {
  # Get pathways that pass methylation data filter
  meth_pass <- meth_results$pathway[
    meth_results$n_genes_tested > meth_min & 
      meth_results$n_genes_tested < meth_max
  ]
  
  # Get pathways that pass discord (meth+expr) data filter
  discord_pass <- discord_results$pathway[
    discord_results$n_genes_tested > discord_min & 
      discord_results$n_genes_tested < discord_max
  ]
  
  passing_pathways <- intersect(meth_pass, discord_pass)
  
  cat(sprintf("Pathways passing methylation data filter (%d < n < %d): %d\n", 
              meth_min, meth_max, length(meth_pass)))
  cat(sprintf("Pathways passing discord data filter (%d < n < %d): %d\n", 
              discord_min, discord_max, length(discord_pass)))
  cat(sprintf("Pathways passing BOTH filters: %d\n", length(passing_pathways)))
  
  return(passing_pathways)
}

# ============================================================================
# Extract pathway ID (for GO)
# ============================================================================
get_pathway_id <- function(pathway_names) {
  sapply(strsplit(pathway_names, "::"), function(x) {
    if (length(x) >= 2) x[2] else NA_character_
  })
}