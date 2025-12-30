#===============================================================================
# 02_regional_analysis.R
#===============================================================================
# DESCRIPTION:
#   Functions for aggregating CpG-level data into genomic regions and 
#   performing regional cis-trans classification analysis.
#
# INPUT FORMATS:
#   BSseq object - Filtered methylation data from 01_data_processing.R
#   GRanges object - Genomic regions from BED file
#
# ADJUSTABLE PARAMETERS:
#   min_cpgs_per_region  - Minimum CpGs required per region (default: 3)
#   min_total_coverage   - Minimum total coverage per region (default: 50)
#   alpha                - Significance threshold (default: 0.05)
#   analyze_heterogeneity - Calculate within-region heterogeneity (default: TRUE)
#
# OUTPUTS:
#   regional_results - Data frame with classification for each region
#   heterogeneity_results - Within-region pattern consistency metrics (support fraction: percentage of individual CpG sites out of all sites with matching classification as pooled regional data)
#   confidence_results - Site-level support for regional classifications (support fraction: percentage of individual CpG sites out of all sites with matching classification as pooled regional data)
#
# DEPENDENCIES:
#   Requires: 05_CIS_TRANS_CORE.R (TwoModelCisTransAnalysis function)
#===============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(data.table)
})

#-------------------------------------------------------------------------------
# REGIONAL DATA AGGREGATION
#-------------------------------------------------------------------------------

AggregateRegionalData <- function(BSobj, regions_gr, min_cpgs_per_region = 3, 
                                  min_total_coverage = 50) {
  # Aggregate CpG sites within regions
  cpg_gr <- granges(BSobj)
  overlaps <- findOverlaps(cpg_gr, regions_gr)
  
  cov_matrix <- getCoverage(BSobj, type = "Cov")
  meth_matrix <- getCoverage(BSobj, type = "M")
  
  n_regions <- length(regions_gr)
  n_samples <- ncol(BSobj)
  sample_names <- sampleNames(BSobj)
  
  regional_cov <- matrix(0, nrow = n_regions, ncol = n_samples)
  regional_meth <- matrix(0, nrow = n_regions, ncol = n_samples)
  regional_cpg_count <- integer(n_regions)
  regional_valid <- logical(n_regions)
  
  colnames(regional_cov) <- sample_names
  colnames(regional_meth) <- sample_names
  
  for(i in 1:n_regions) {
    region_cpgs <- queryHits(overlaps)[subjectHits(overlaps) == i]
    regional_cpg_count[i] <- length(region_cpgs)
    
    if(length(region_cpgs) == 0) next
    
    if(length(region_cpgs) == 1) {
      regional_meth[i, ] <- meth_matrix[region_cpgs, ]
      regional_cov[i, ] <- cov_matrix[region_cpgs, ]
    } else {
      regional_meth[i, ] <- colSums(meth_matrix[region_cpgs, ], na.rm = TRUE)
      regional_cov[i, ] <- colSums(cov_matrix[region_cpgs, ], na.rm = TRUE)
    }
    
    regional_valid[i] <- (regional_cpg_count[i] >= min_cpgs_per_region) && 
      (sum(regional_cov[i, ]) >= min_total_coverage)
  }
  
  valid_regions <- which(regional_valid)
  if(length(valid_regions) == 0) stop("No regions meet minimum requirements")
  
  regional_BSobj <- BSseq(
    gr = regions_gr[valid_regions],
    M = regional_meth[valid_regions, , drop = FALSE],
    Cov = regional_cov[valid_regions, , drop = FALSE],
    sampleNames = sample_names
  )
  
  region_metadata <- data.frame(
    region_index = valid_regions,
    region_name = if("name" %in% colnames(mcols(regions_gr))) {
      mcols(regions_gr)$name[valid_regions]
    } else paste0("region_", valid_regions),
    chr = as.character(seqnames(regions_gr[valid_regions])),
    start = start(regions_gr[valid_regions]),
    end = end(regions_gr[valid_regions]),
    width = width(regions_gr[valid_regions]),
    n_cpgs = regional_cpg_count[valid_regions],
    total_coverage = rowSums(regional_cov[valid_regions, ]),
    stringsAsFactors = FALSE
  )
  
  return(list(
    regional_BSobj = regional_BSobj,
    metadata = region_metadata,
    all_regions_valid = regional_valid,
    all_regions_cpg_count = regional_cpg_count,
    cpg_region_mapping = overlaps
  ))
}

#-------------------------------------------------------------------------------
# HETEROGENEITY ANALYSIS
#-------------------------------------------------------------------------------

AnalyzeRegionalHeterogeneity <- function(BSobj, regions_gr, site_level_results = NULL) {
  # Analyze heterogeneity of cis-trans patterns within regions
  cpg_gr <- granges(BSobj)
  overlaps <- findOverlaps(cpg_gr, regions_gr)
  
  cov_matrix <- getCoverage(BSobj, type = "Cov")
  meth_matrix <- getCoverage(BSobj, type = "M")
  meth_prop <- meth_matrix / cov_matrix
  meth_prop[is.na(meth_prop)] <- 0
  
  sample_names <- sampleNames(BSobj)
  parH_samples <- grep("parH", sample_names, value = TRUE, ignore.case = TRUE)
  parC_samples <- grep("parC", sample_names, value = TRUE, ignore.case = TRUE)
  hybH_samples <- grep("hybH", sample_names, value = TRUE, ignore.case = TRUE)
  hybC_samples <- grep("hybC", sample_names, value = TRUE, ignore.case = TRUE)
  
  if(length(parH_samples) == 0 || length(parC_samples) == 0 || 
     length(hybH_samples) == 0 || length(hybC_samples) == 0) {
    return(NULL)
  }
  
  parent_HC_diff <- rowMeans(meth_prop[, parH_samples, drop = FALSE]) - 
    rowMeans(meth_prop[, parC_samples, drop = FALSE])
  hybrid_HC_diff <- rowMeans(meth_prop[, hybH_samples, drop = FALSE]) - 
    rowMeans(meth_prop[, hybC_samples, drop = FALSE])
  
  n_regions <- length(regions_gr)
  heterogeneity_results <- data.frame(
    region_index = 1:n_regions,
    n_cpgs_in_region = integer(n_regions),
    parent_diff_mean = numeric(n_regions),
    parent_diff_sd = numeric(n_regions),
    parent_diff_cv = numeric(n_regions),
    hybrid_diff_mean = numeric(n_regions),
    hybrid_diff_sd = numeric(n_regions),
    hybrid_diff_cv = numeric(n_regions),
    pattern_consistency = character(n_regions),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:n_regions) {
    region_cpgs <- queryHits(overlaps)[subjectHits(overlaps) == i]
    heterogeneity_results$n_cpgs_in_region[i] <- length(region_cpgs)
    
    if(length(region_cpgs) == 0) {
      heterogeneity_results$pattern_consistency[i] <- "no_cpgs"
      next
    }
    
    if(length(region_cpgs) == 1) {
      heterogeneity_results$pattern_consistency[i] <- "single_cpg"
      heterogeneity_results$parent_diff_mean[i] <- parent_HC_diff[region_cpgs]
      heterogeneity_results$hybrid_diff_mean[i] <- hybrid_HC_diff[region_cpgs]
      next
    }
    
    region_parent_diffs <- parent_HC_diff[region_cpgs]
    region_hybrid_diffs <- hybrid_HC_diff[region_cpgs]
    region_parent_diffs <- region_parent_diffs[!is.na(region_parent_diffs)]
    region_hybrid_diffs <- region_hybrid_diffs[!is.na(region_hybrid_diffs)]
    
    if(length(region_parent_diffs) == 0 || length(region_hybrid_diffs) == 0) {
      heterogeneity_results$pattern_consistency[i] <- "insufficient_data"
      next
    }
    
    heterogeneity_results$parent_diff_mean[i] <- mean(region_parent_diffs)
    heterogeneity_results$parent_diff_sd[i] <- sd(region_parent_diffs)
    heterogeneity_results$hybrid_diff_mean[i] <- mean(region_hybrid_diffs)
    heterogeneity_results$hybrid_diff_sd[i] <- sd(region_hybrid_diffs)
    
    heterogeneity_results$parent_diff_cv[i] <- ifelse(
      abs(heterogeneity_results$parent_diff_mean[i]) > 0.01,
      abs(heterogeneity_results$parent_diff_sd[i] / heterogeneity_results$parent_diff_mean[i]), NA)
    heterogeneity_results$hybrid_diff_cv[i] <- ifelse(
      abs(heterogeneity_results$hybrid_diff_mean[i]) > 0.01,
      abs(heterogeneity_results$hybrid_diff_sd[i] / heterogeneity_results$hybrid_diff_mean[i]), NA)
    
    avg_cv <- mean(c(heterogeneity_results$parent_diff_cv[i], 
                     heterogeneity_results$hybrid_diff_cv[i]), na.rm = TRUE)
    
    if(is.na(avg_cv)) {
      heterogeneity_results$pattern_consistency[i] <- "low_signal"
    } else if(avg_cv < 0.5) {
      heterogeneity_results$pattern_consistency[i] <- "consistent"
    } else if(avg_cv < 1.0) {
      heterogeneity_results$pattern_consistency[i] <- "moderately_heterogeneous"
    } else {
      heterogeneity_results$pattern_consistency[i] <- "highly_heterogeneous"
    }
  }
  
  return(heterogeneity_results)
}

#-------------------------------------------------------------------------------
# MAIN REGIONAL ANALYSIS FUNCTION
#-------------------------------------------------------------------------------

RegionalCisTransAnalysis <- function(BSobj, regions_gr, 
                                     min_cpgs_per_region = 3,
                                     min_total_coverage = 50,
                                     alpha = 0.05,
                                     analyze_heterogeneity = TRUE) {
  # Main function for regional cis-trans analysis
  regional_data <- AggregateRegionalData(BSobj, regions_gr, min_cpgs_per_region, 
                                         min_total_coverage)
  
  heterogeneity_results <- NULL
  if(analyze_heterogeneity) {
    heterogeneity_results <- AnalyzeRegionalHeterogeneity(BSobj, regions_gr)
  }
  
  regional_sample_names <- sampleNames(regional_data$regional_BSobj)
  if(!ValidateSampleNamesForCisTrans(regional_sample_names, verbose = FALSE)) {
    stop("Regional BSseq object has invalid sample names")
  }
  
  regional_cistrans_results <- TwoModelCisTransAnalysis(
    BSobj = regional_data$regional_BSobj,
    alpha = alpha,
    min_coverage = 1,
    verbose = FALSE
  )
  
  final_results <- CombineRegionalResults(
    cistrans_results = regional_cistrans_results,
    regional_metadata = regional_data$metadata,
    heterogeneity_results = heterogeneity_results,
    all_regions_gr = regions_gr,
    all_regions_valid = regional_data$all_regions_valid,
    min_cpgs_per_region = min_cpgs_per_region,
    min_total_coverage = min_total_coverage
  )
  
  return(final_results)
}

#-------------------------------------------------------------------------------
# COMBINE REGIONAL RESULTS
#-------------------------------------------------------------------------------

CombineRegionalResults <- function(cistrans_results, regional_metadata,
                                   heterogeneity_results, all_regions_gr,
                                   all_regions_valid, min_cpgs_per_region,
                                   min_total_coverage) {
  # Combine all regional results into single data frame
  n_total_regions <- length(all_regions_gr)
  
  region_results <- data.frame(
    region_id = 1:n_total_regions,
    region_name = if("name" %in% colnames(mcols(all_regions_gr))) {
      mcols(all_regions_gr)$name
    } else paste0("region_", 1:n_total_regions),
    chr = as.character(seqnames(all_regions_gr)),
    start = start(all_regions_gr),
    end = end(all_regions_gr),
    width = width(all_regions_gr),
    meets_requirements = all_regions_valid,
    failure_reason = character(n_total_regions),
    classification = "unclassified",
    p_hybrid = NA, p_interaction = NA,
    hybrid_effect = NA, interaction_effect = NA,
    parent_HC_diff = NA, hybrid_HC_diff = NA,
    n_cpgs = 0, total_coverage = 0,
    stringsAsFactors = FALSE
  )
  
  if(!is.null(heterogeneity_results)) {
    region_results$pattern_consistency <- heterogeneity_results$pattern_consistency
    region_results$parent_diff_cv <- heterogeneity_results$parent_diff_cv
    region_results$hybrid_diff_cv <- heterogeneity_results$hybrid_diff_cv
    region_results$n_cpgs <- heterogeneity_results$n_cpgs_in_region
    
    no_cpg_regions <- heterogeneity_results$n_cpgs_in_region == 0
    insufficient_cpg_regions <- heterogeneity_results$n_cpgs_in_region > 0 & 
      heterogeneity_results$n_cpgs_in_region < min_cpgs_per_region
    
    region_results$failure_reason[no_cpg_regions] <- "no_cpgs"
    region_results$failure_reason[insufficient_cpg_regions] <- 
      paste0("too_few_cpgs_(<", min_cpgs_per_region, ")")
  }
  
  valid_region_indices <- which(all_regions_valid)
  
  if(length(valid_region_indices) > 0) {
    for(i in seq_along(valid_region_indices)) {
      orig_idx <- valid_region_indices[i]
      if(i <= nrow(cistrans_results$cpg_results)) {
        cistrans_row <- cistrans_results$cpg_results[i, ]
        region_results$classification[orig_idx] <- cistrans_row$classification
        region_results$p_hybrid[orig_idx] <- cistrans_row$p_hybrid
        region_results$p_interaction[orig_idx] <- cistrans_row$p_interaction
        region_results$hybrid_effect[orig_idx] <- cistrans_row$hybrid_effect
        region_results$interaction_effect[orig_idx] <- cistrans_row$interaction_effect
        region_results$parent_HC_diff[orig_idx] <- cistrans_row$parent_HC_diff
        region_results$hybrid_HC_diff[orig_idx] <- cistrans_row$hybrid_HC_diff
      }
      if(i <= nrow(regional_metadata)) {
        region_results$n_cpgs[orig_idx] <- regional_metadata$n_cpgs[i]
        region_results$total_coverage[orig_idx] <- regional_metadata$total_coverage[i]
      }
    }
  }
  
  region_results$failure_reason[region_results$meets_requirements] <- "none"
  
  return(list(
    region_results = region_results,
    cistrans_results = cistrans_results,
    regional_metadata = regional_metadata,
    heterogeneity_results = heterogeneity_results,
    summary_stats = list(
      total_regions = n_total_regions,
      regions_with_cpgs = sum(region_results$n_cpgs > 0),
      regions_meeting_requirements = sum(region_results$meets_requirements),
      classified_regions = sum(region_results$meets_requirements & 
                                 region_results$classification != "unclassified")
    )
  ))
}

#-------------------------------------------------------------------------------
# CONFIDENCE METRICS
#-------------------------------------------------------------------------------

CalculateRegionalConfidence <- function(BSobj, regions_gr, regional_results, 
                                        site_level_results = NULL) {
  # Calculate confidence metrics for regional classifications
  if(is.null(site_level_results)) {
    site_level_results <- TwoModelCisTransAnalysis(BSobj, alpha = 0.05, verbose = FALSE)
  }
  
  cpg_gr <- granges(BSobj)
  overlaps <- findOverlaps(cpg_gr, regions_gr)
  
  confidence_metrics <- data.frame(
    region_id = 1:length(regions_gr),
    n_supporting_sites = integer(length(regions_gr)),
    n_opposing_sites = integer(length(regions_gr)),
    n_unclassified_sites = integer(length(regions_gr)),
    support_fraction = numeric(length(regions_gr)),
    confidence_level = character(length(regions_gr)),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:length(regions_gr)) {
    region_cpgs <- queryHits(overlaps)[subjectHits(overlaps) == i]
    
    if(length(region_cpgs) == 0) {
      confidence_metrics$confidence_level[i] <- "no_sites"
      next
    }
    
    site_classifications <- site_level_results$cpg_results$classification[region_cpgs]
    regional_classification <- regional_results$region_results$classification[i]
    
    if(regional_classification == "unclassified") {
      confidence_metrics$confidence_level[i] <- "region_unclassified"
      next
    }
    
    supporting <- sum(site_classifications == regional_classification, na.rm = TRUE)
    unclassified <- sum(site_classifications == "unclassified", na.rm = TRUE)
    opposing <- length(site_classifications) - supporting - unclassified
    
    confidence_metrics$n_supporting_sites[i] <- supporting
    confidence_metrics$n_opposing_sites[i] <- opposing
    confidence_metrics$n_unclassified_sites[i] <- unclassified
    
    classified_sites <- supporting + opposing
    confidence_metrics$support_fraction[i] <- if(classified_sites > 0) {
      supporting / classified_sites
    } else 0
    
    # change confidence metrics here
    if(classified_sites == 0) {
      confidence_metrics$confidence_level[i] <- "no_classified_sites"
    } else if(confidence_metrics$support_fraction[i] >= 0.66) {
      confidence_metrics$confidence_level[i] <- "high_confidence"
    } else if(confidence_metrics$support_fraction[i] >= 0.33) {
      confidence_metrics$confidence_level[i] <- "medium_confidence"
    } else {
      confidence_metrics$confidence_level[i] <- "low_confidence"
    }
  }
  
  return(confidence_metrics)
}