#===============================================================================
# 05_CIS_TRANS_CORE.R
#===============================================================================
# DESCRIPTION:
#   Two-model framework for cis-trans regulatory classification of methylation
#   differences between parental and hybrid alleles. Treats biological replicates
#   as separate observations for increased statistical power.
#
# INPUT FORMATS:
#   BSseq object with samples named following patterns:
#     - parH* (parental H allele, e.g., parH_rep1, parH_rep2)
#     - parC* (parental C allele, e.g., parC_rep1)
#     - hybH* (hybrid H allele, e.g., hybH_rep1)
#     - hybC* (hybrid C allele, e.g., hybC_rep1, hybC_rep2)
#
# ADJUSTABLE PARAMETERS:
#   alpha        - Significance threshold for classification (default: 0.05)
#   min_coverage - Minimum read coverage per CpG site (default: 10)
#
# OUTPUTS:
#   cpg_results data frame with columns:
#     - chr, pos: Genomic coordinates
#     - classification: pure_cis, pure_trans, cis_plus_trans, cis_x_trans, conserved
#     - p_hybrid: P-value for hybrid H vs C effect
#     - p_interaction: P-value for allele x generation interaction
#     - parent_HC_diff, hybrid_HC_diff: Methylation differences
#
# CLASSIFICATION RULES:
#   - pure_cis: p_hybrid < alpha, p_interaction >= alpha
#   - pure_trans: p_hybrid >= alpha, p_interaction < alpha
#   - cis_plus_trans: both < alpha, effects same direction
#   - cis_x_trans: both < alpha, effects opposite direction
#   - conserved: both >= alpha
#
# ACKNOWLEDGMENTS:
#   Beta-binomial fitting functions adapted from DSS package
#   (DOI: 10.18129/B9.bioc.DSS)
#===============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

#===============================================================================
# BETA-BINOMIAL FITTING FUNCTIONS
# Adapted from DSS package (DOI: 10.18129/B9.bioc.DSS)
#===============================================================================

dbb <- function(size, x, mu, phi, log = TRUE) {
  # Beta-binomial density function (from DSS)
  size <- as.array(size)
  x <- as.array(x)
  
  tmp <- 1/phi - 1
  alpha <- mu * tmp
  beta <- tmp - alpha
  
  v <- lchoose(size, x) - lbeta(beta, alpha) + lbeta(size - x + beta, x + alpha)
  
  if(!log) return(exp(v))
  return(v)
}

BetaBinomialFit.oneCG <- function(Y, N, X, Z, n, p) {
  # Fit beta-binomial model for one CpG (adapted from DSS)
  c1 <- 0.001
  
  ix <- N > 0
  n_valid <- sum(ix)
  
  if(n_valid < p + 1) return(NULL)
  
  if(n_valid < length(ix)) {
    X <- X[ix, , drop = FALSE]
    Y <- Y[ix]
    N <- N[ix]
    Z <- Z[ix]
  }
  
  if(nrow(X) < ncol(X) + 1) return(NULL)
  if(any(abs(svd(X)$d) < 1e-8)) return(NULL)
  
  XTVinv <- t(X * N)
  tryCatch({
    beta0 <- solve(XTVinv %*% X) %*% (XTVinv %*% Z)
  }, error = function(e) return(NULL))
  
  if(is.null(beta0)) return(NULL)
  
  residuals <- (Z - X %*% beta0)^2 * N
  df <- n_valid - p
  if(df <= 0) {
    phiHat <- 0.1
  } else {
    phiHat <- (sum(residuals) - df) * n_valid / df / sum(N - 1)
  }
  phiHat <- min(max(c1, phiHat), 1 - c1)
  
  weights <- N / (1 + (N - 1) * phiHat)
  XTVinv <- t(X * weights)
  
  tryCatch({
    XTVinvX.inv <- solve(XTVinv %*% X)
    beta0 <- XTVinvX.inv %*% (XTVinv %*% Z)
    se.beta0 <- sqrt(diag(XTVinvX.inv))
  }, error = function(e) return(NULL))
  
  if(is.null(beta0)) return(NULL)
  
  return(list(
    beta0 = as.vector(beta0),
    se.beta0 = as.vector(se.beta0),
    var.beta0 = XTVinvX.inv,
    phi = phiHat
  ))
}

#===============================================================================
# SAMPLE IDENTIFICATION
#===============================================================================

IdentifyAndGroupSamples <- function(sample_names, verbose = TRUE) {
  # Identify sample groups from naming patterns
  parH_samples <- grep("parH|parent.*H|H.*parent", sample_names, value = TRUE, ignore.case = TRUE)
  parC_samples <- grep("parC|parent.*C|C.*parent", sample_names, value = TRUE, ignore.case = TRUE)
  hybH_samples <- grep("hybH|hybrid.*H|H.*hybrid", sample_names, value = TRUE, ignore.case = TRUE)
  hybC_samples <- grep("hybC|hybrid.*C|C.*hybrid", sample_names, value = TRUE, ignore.case = TRUE)
  
  missing_groups <- c()
  if(length(parH_samples) == 0) missing_groups <- c(missing_groups, "parH")
  if(length(parC_samples) == 0) missing_groups <- c(missing_groups, "parC")
  if(length(hybH_samples) == 0) missing_groups <- c(missing_groups, "hybH")
  if(length(hybC_samples) == 0) missing_groups <- c(missing_groups, "hybC")
  
  if(length(missing_groups) > 0) {
    stop("Missing sample groups: ", paste(missing_groups, collapse = ", "))
  }
  
  sample_groups <- list(
    parH = parH_samples,
    parC = parC_samples,
    hybH = hybH_samples,
    hybC = hybC_samples
  )
  
  if(verbose) {
    cat("Sample groups identified:\n")
    for(group in names(sample_groups)) {
      cat("  ", group, ":", paste(sample_groups[[group]], collapse = ", "), "\n")
    }
  }
  
  return(sample_groups)
}

ValidateSampleNamesForCisTrans <- function(sample_names, verbose = TRUE) {
  # Validate sample names contain required patterns
  tryCatch({
    IdentifyAndGroupSamples(sample_names, verbose = FALSE)
    return(TRUE)
  }, error = function(e) {
    if(verbose) cat("Validation failed:", e$message, "\n")
    return(FALSE)
  })
}

#===============================================================================
# MODEL FITTING FUNCTIONS
#===============================================================================

FitHybridOnlyModel <- function(BSobj, sample_groups = NULL) {
  # Model 1: Hybrid-only H vs C comparison
  if(is.null(sample_groups)) {
    sample_groups <- IdentifyAndGroupSamples(sampleNames(BSobj), verbose = FALSE)
  }
  
  hybrid_samples <- c(sample_groups$hybH, sample_groups$hybC)
  hybrid_indices <- which(sampleNames(BSobj) %in% hybrid_samples)
  BSobj_hybrid <- BSobj[, hybrid_indices]
  
  N0 <- getCoverage(BSobj_hybrid, type = "Cov")
  Y0 <- getCoverage(BSobj_hybrid, type = "M")
  
  hybrid_sample_names <- sampleNames(BSobj_hybrid)
  allele_assignments <- ifelse(hybrid_sample_names %in% sample_groups$hybH, "H", "C")
  
  design_hybrid <- data.frame(
    sample = hybrid_sample_names,
    allele = factor(allele_assignments, levels = c("C", "H")),
    stringsAsFactors = FALSE
  )
  
  X_hybrid <- model.matrix(~ allele, design_hybrid)
  
  c0 <- 0.1
  Z0_hybrid <- asin(2 * (Y0 + c0) / (N0 + 2 * c0) - 1)
  
  n_sites <- nrow(Y0)
  p <- ncol(X_hybrid)
  n <- nrow(X_hybrid)
  
  beta_hybrid <- matrix(NA, nrow = n_sites, ncol = p)
  var_beta_hybrid <- array(NA, dim = c(n_sites, p, p))
  phi_hybrid <- numeric(n_sites)
  
  for(i in 1:n_sites) {
    tmp <- BetaBinomialFit.oneCG(Y0[i, ], N0[i, ], X_hybrid, Z0_hybrid[i, ], n, p)
    if(is.null(tmp)) next
    beta_hybrid[i, ] <- tmp$beta0
    var_beta_hybrid[i, , ] <- tmp$var.beta0
    phi_hybrid[i] <- tmp$phi
  }
  
  return(list(beta = beta_hybrid, var_beta = var_beta_hybrid, phi = phi_hybrid,
              design = design_hybrid, X = X_hybrid, BSobj_used = BSobj_hybrid))
}

FitFullInteractionModel <- function(BSobj, sample_groups = NULL) {
  # Model 2: Full interaction model (allele + generation + interaction)
  if(is.null(sample_groups)) {
    sample_groups <- IdentifyAndGroupSamples(sampleNames(BSobj), verbose = FALSE)
  }
  
  all_sample_names <- sampleNames(BSobj)
  allele_assignments <- character(length(all_sample_names))
  generation_assignments <- character(length(all_sample_names))
  
  for(i in seq_along(all_sample_names)) {
    sample_name <- all_sample_names[i]
    if(sample_name %in% sample_groups$parH) {
      allele_assignments[i] <- "H"; generation_assignments[i] <- "parent"
    } else if(sample_name %in% sample_groups$parC) {
      allele_assignments[i] <- "C"; generation_assignments[i] <- "parent"
    } else if(sample_name %in% sample_groups$hybH) {
      allele_assignments[i] <- "H"; generation_assignments[i] <- "hybrid"
    } else if(sample_name %in% sample_groups$hybC) {
      allele_assignments[i] <- "C"; generation_assignments[i] <- "hybrid"
    }
  }
  
  design_full <- data.frame(
    sample = all_sample_names,
    allele = factor(allele_assignments, levels = c("C", "H")),
    generation = factor(generation_assignments, levels = c("parent", "hybrid")),
    stringsAsFactors = FALSE
  )
  
  X_full <- model.matrix(~ allele + generation + allele:generation, design_full)
  
  N0 <- getCoverage(BSobj, type = "Cov")
  Y0 <- getCoverage(BSobj, type = "M")
  
  c0 <- 0.1
  Z0_full <- asin(2 * (Y0 + c0) / (N0 + 2 * c0) - 1)
  
  n_sites <- nrow(Y0)
  p <- ncol(X_full)
  n <- nrow(X_full)
  
  beta_full <- matrix(NA, nrow = n_sites, ncol = p)
  var_beta_full <- array(NA, dim = c(n_sites, p, p))
  phi_full <- numeric(n_sites)
  
  for(i in 1:n_sites) {
    tmp <- BetaBinomialFit.oneCG(Y0[i, ], N0[i, ], X_full, Z0_full[i, ], n, p)
    if(is.null(tmp)) next
    beta_full[i, ] <- tmp$beta0
    var_beta_full[i, , ] <- tmp$var.beta0
    phi_full[i] <- tmp$phi
  }
  
  return(list(beta = beta_full, var_beta = var_beta_full, phi = phi_full,
              design = design_full, X = X_full, BSobj_used = BSobj))
}

#===============================================================================
# HYPOTHESIS TESTING
#===============================================================================

TestModelContrasts <- function(model_fit, contrast_matrix) {
  # Test linear contrasts from fitted model
  p <- ncol(model_fit$X)
  betas <- model_fit$beta
  
  if(nrow(contrast_matrix) != p) {
    stop("Contrast matrix dimension mismatch")
  }
  
  stat <- rep(NA, nrow(betas))
  
  for(i in 1:nrow(betas)) {
    if(any(is.na(betas[i, ]))) next
    Sigma <- model_fit$var_beta[i, , ]
    if(any(is.na(Sigma)) || any(eigen(Sigma)$values <= 0)) next
    
    tryCatch({
      contrast_beta <- t(contrast_matrix) %*% betas[i, ]
      contrast_var <- t(contrast_matrix) %*% Sigma %*% contrast_matrix
      if(any(eigen(contrast_var)$values <= 0)) next
      stat[i] <- as.numeric(t(contrast_beta) %*% solve(contrast_var) %*% contrast_beta)
    }, error = function(e) NULL)
  }
  
  df <- ncol(contrast_matrix)
  pvals <- 1 - pchisq(stat, df = df)
  
  if(ncol(contrast_matrix) == 1) {
    effect_size <- betas %*% contrast_matrix
  } else {
    effect_size <- rep(NA, length(pvals))
  }
  
  return(list(stat = stat, pvals = pvals, effect_size = effect_size))
}

#===============================================================================
# MAIN ANALYSIS FUNCTION
#===============================================================================

TwoModelCisTransAnalysis <- function(BSobj, alpha = 0.05, min_coverage = 10, verbose = TRUE) {
  # Main two-model cis-trans classification function
  if(!is(BSobj, "BSseq")) stop("BSobj must be a BSseq object")
  if(alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  
  sample_names <- sampleNames(BSobj)
  if(!ValidateSampleNamesForCisTrans(sample_names, verbose = verbose)) {
    stop("Invalid sample names")
  }
  
  sample_groups <- IdentifyAndGroupSamples(sample_names, verbose = verbose)
  
  # Filter by coverage
  coverage_mask <- rowSums(getCoverage(BSobj, type = "Cov") >= min_coverage) == ncol(BSobj)
  BSobj_filtered <- BSobj[coverage_mask, ]
  
  if(verbose) cat("CpGs with sufficient coverage:", sum(coverage_mask), "\n")
  if(sum(coverage_mask) == 0) stop("No CpGs with sufficient coverage")
  
  # Fit models
  model1_fit <- FitHybridOnlyModel(BSobj_filtered, sample_groups)
  hybrid_contrast <- matrix(c(0, 1), ncol = 1)
  model1_test <- TestModelContrasts(model1_fit, hybrid_contrast)
  
  model2_fit <- FitFullInteractionModel(BSobj_filtered, sample_groups)
  n_params <- ncol(model2_fit$X)
  interaction_contrast <- matrix(c(rep(0, n_params-1), 1), ncol = 1)
  model2_test <- TestModelContrasts(model2_fit, interaction_contrast)
  
  # Classification
  classification_results <- ClassifyTwoModelResults(
    BSobj_filtered, model1_test, model2_test, alpha, sample_groups, verbose
  )
  
  gr <- granges(BSobj_filtered)
  final_results <- data.frame(
    chr = as.character(seqnames(gr)),
    pos = start(gr),
    classification_results,
    stringsAsFactors = FALSE
  )
  
  summary_stats <- table(final_results$classification)
  
  return(list(
    cpg_results = final_results,
    model1_fit = model1_fit,
    model2_fit = model2_fit,
    model1_test = model1_test,
    model2_test = model2_test,
    summary = summary_stats,
    sample_groups = sample_groups
  ))
}

ClassifyTwoModelResults <- function(BSobj, model1_test, model2_test, alpha = 0.05, 
                                    sample_groups, verbose = TRUE) {
  # Classify CpGs based on two-model results
  p_hybrid <- model1_test$pvals
  p_interaction <- model2_test$pvals
  p_hybrid_adj <- p.adjust(p_hybrid, method = "BH")
  p_interaction_adj <- p.adjust(p_interaction, method = "BH")
  
  hybrid_effect <- model1_test$effect_size
  interaction_effect <- model2_test$effect_size
  
  # Calculate biological differences
  M_matrix <- getCoverage(BSobj, type = "M")
  Cov_matrix <- getCoverage(BSobj, type = "Cov")
  meth_proportions <- M_matrix / Cov_matrix
  meth_proportions[is.na(meth_proportions)] <- 0
  
  sample_names <- sampleNames(BSobj)
  parH_idx <- which(sample_names %in% sample_groups$parH)
  parC_idx <- which(sample_names %in% sample_groups$parC)
  hybH_idx <- which(sample_names %in% sample_groups$hybH)
  hybC_idx <- which(sample_names %in% sample_groups$hybC)
  
  parent_HC_diff <- rowMeans(meth_proportions[, parH_idx, drop = FALSE]) - 
    rowMeans(meth_proportions[, parC_idx, drop = FALSE])
  hybrid_HC_diff <- rowMeans(meth_proportions[, hybH_idx, drop = FALSE]) - 
    rowMeans(meth_proportions[, hybC_idx, drop = FALSE])
  
  decision_metric <- (parent_HC_diff - hybrid_HC_diff) * hybrid_HC_diff
  
  # Classification
  classification <- rep("unclassified", length(p_hybrid))
  sig_hybrid <- p_hybrid < alpha & !is.na(p_hybrid)
  sig_interaction <- p_interaction < alpha & !is.na(p_interaction)
  
  classification[sig_hybrid & !sig_interaction] <- "pure_cis"
  classification[!sig_hybrid & !sig_interaction] <- "conserved"
  classification[!sig_hybrid & sig_interaction] <- "pure_trans"
  
  both_sig <- sig_hybrid & sig_interaction
  substantial_effect <- abs(hybrid_HC_diff) > 0.01
  
  classification[both_sig & substantial_effect & decision_metric > 0] <- "cis_plus_trans"
  classification[both_sig & substantial_effect & decision_metric < 0] <- "cis_x_trans"
  
  return(data.frame(
    classification = classification,
    p_hybrid = p_hybrid,
    p_interaction = p_interaction,
    p_hybrid_adj = p_hybrid_adj,
    p_interaction_adj = p_interaction_adj,
    hybrid_effect = hybrid_effect,
    interaction_effect = interaction_effect,
    parent_HC_diff = parent_HC_diff,
    hybrid_HC_diff = hybrid_HC_diff,
    decision_metric = decision_metric,
    sig_hybrid = sig_hybrid,
    sig_interaction = sig_interaction,
    stringsAsFactors = FALSE
  ))
}

#===============================================================================
# PLOTTING FUNCTIONS
#===============================================================================

PlotTwoModelResults <- function(results, title = "Two-Model Cis-Trans Analysis") {
  # Scatter plot of classification results
  plot_data <- results$cpg_results[results$cpg_results$classification != "unclassified", ]
  n_unclassified <- sum(results$cpg_results$classification == "unclassified")
  
  if(nrow(plot_data) == 0) return(NULL)
  
  colors <- c("conserved" = "#999999", "pure_cis" = "#1F78B4", "pure_trans" = "#33A02C",
              "cis_plus_trans" = "#FF7F00", "cis_x_trans" = "#E31A1C")
  
  p <- ggplot(plot_data, aes(x = parent_HC_diff, y = hybrid_HC_diff, color = classification)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors, name = "Regulation Type") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right") +
    labs(title = title, subtitle = paste0(n_unclassified, " unclassified excluded"),
         x = "H-C Difference (Parents)", y = "H-C Difference (Hybrids)") +
    coord_fixed(ratio = 1)
  
  return(list(plot = p, data = plot_data, n_unclassified = n_unclassified))
}

PlotTwoModelSummary <- function(results) {
  # Summary bar plot
  classified_data <- results$cpg_results[results$cpg_results$classification != "unclassified", ]
  if(nrow(classified_data) == 0) return(NULL)
  
  summary_counts <- table(classified_data$classification)
  plot_data <- data.frame(
    Classification = names(summary_counts),
    Count = as.numeric(summary_counts),
    Proportion = as.numeric(summary_counts) / sum(summary_counts)
  )
  
  p <- ggplot(plot_data, aes(x = Classification, y = Proportion, fill = Classification)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "Regulatory Classification Summary",
         subtitle = paste0(nrow(classified_data), " classified sites"),
         x = "Regulatory Pattern", y = "Proportion") +
    scale_fill_brewer(type = "qual", palette = "Set3")
  
  return(p)
}

PlotTwoModelResultsLog2FC <- function(results, title = "Two-Model Analysis (Log2FC)", 
                                      pseudocount = 0.01) {
  # Log2FC scatter plot
  plot_data <- results$cpg_results[results$cpg_results$classification != "unclassified", ]
  n_unclassified <- sum(results$cpg_results$classification == "unclassified")
  
  if(nrow(plot_data) == 0) return(NULL)
  
  baseline_meth <- 0.5
  plot_data$parent_log2FC <- log2(
    pmax(pseudocount, pmin(1-pseudocount, baseline_meth + plot_data$parent_HC_diff/2)) /
      pmax(pseudocount, pmin(1-pseudocount, baseline_meth - plot_data$parent_HC_diff/2))
  )
  plot_data$hybrid_log2FC <- log2(
    pmax(pseudocount, pmin(1-pseudocount, baseline_meth + plot_data$hybrid_HC_diff/2)) /
      pmax(pseudocount, pmin(1-pseudocount, baseline_meth - plot_data$hybrid_HC_diff/2))
  )
  
  colors <- c("conserved" = "#999999", "pure_cis" = "#1F78B4", "pure_trans" = "#33A02C",
              "cis_plus_trans" = "#FF7F00", "cis_x_trans" = "#E31A1C")
  
  p <- ggplot(plot_data, aes(x = parent_log2FC, y = hybrid_log2FC, color = classification)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors, name = "Regulation Type") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right") +
    labs(title = title, subtitle = paste0(n_unclassified, " unclassified excluded"),
         x = "Log2(H/C) Parents", y = "Log2(H/C) Hybrids") +
    coord_fixed(ratio = 1)
  
  return(list(plot = p, data = plot_data, n_unclassified = n_unclassified))
}

#===============================================================================
# TEST DATA GENERATION
#===============================================================================

GenerateTestDataWithReplicates <- function(n_cpgs = 2000, seed = 123) {
  # Generate test data with replicate structure
  set.seed(seed)
  
  chr_names <- paste0("chr", rep(1:3, each = n_cpgs/3))
  positions <- rep(seq(1000, by = 100, length.out = n_cpgs/3), 3)
  
  gr <- GRanges(seqnames = chr_names, ranges = IRanges(start = positions, width = 1))
  
  sample_names <- c("parH_rep1", "parH_rep2", "parC", "hybH", "hybC_rep1", "hybC_rep2")
  n_samples <- length(sample_names)
  
  coverage_matrix <- matrix(rpois(n_cpgs * n_samples, lambda = 40), 
                            nrow = n_cpgs, ncol = n_samples, dimnames = list(NULL, sample_names))
  
  meth_matrix <- matrix(0, nrow = n_cpgs, ncol = n_samples, dimnames = list(NULL, sample_names))
  
  for(i in 1:n_cpgs) {
    base_meth <- rbeta(1, 2, 2)
    cis_effect <- rnorm(1, 0, 0.15)
    trans_effect <- rnorm(1, 0, 0.1)
    
    sample_effects <- list(
      parH_rep1 = base_meth + cis_effect + rnorm(1, 0, 0.03),
      parH_rep2 = base_meth + cis_effect + rnorm(1, 0, 0.03),
      parC = base_meth + rnorm(1, 0, 0.03),
      hybH = base_meth + cis_effect + trans_effect + rnorm(1, 0, 0.03),
      hybC_rep1 = base_meth + rnorm(1, 0, 0.03),
      hybC_rep2 = base_meth + rnorm(1, 0, 0.03)
    )
    
    for(sample_name in sample_names) {
      sample_meth <- max(0.01, min(0.99, sample_effects[[sample_name]]))
      meth_matrix[i, sample_name] <- rbinom(1, coverage_matrix[i, sample_name], sample_meth)
    }
  }
  
  BSobj <- BSseq(gr = gr, M = meth_matrix, Cov = coverage_matrix, sampleNames = sample_names)
  
  return(list(BSobj = BSobj, coverage = coverage_matrix, methylation = meth_matrix))
}