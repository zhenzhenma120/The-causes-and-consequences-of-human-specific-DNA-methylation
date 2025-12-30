#!/usr/bin/env Rscript
# ==============================================================================
# DMR DETECTION - binary scoring with CROPS penalty optimization
# ==============================================================================
#
# DESCRIPTION:
#   Identifies DMRs using changepoint detection with CROPS-based penalty
#   optimization and BINARY scoring. Finds optimal penalty for each genomic segment.
#
#   Scoring: Binary (1 for target class, 0 for all others)
#
# INPUT FORMAT:
#   CSV file with columns:
#     - chr: Chromosome number (numeric)
#     - pos: Genomic position (integer)
#     - classification: Target class label
#     - parent_HC_diff: Methylation difference (-1 to 1)
#
# OUTPUT FORMAT:
#   CSV/RDS with columns:
#     chr, start, end, n_cpgs, target_class, target_proportion,
#     mean_effect, median_effect, effect_sd, dominant_direction,
#     directional_consistency, optimal_penalty
#
# USAGE:
#   Rscript 02_dmr_binary_crops.R [maxgap] [celltype]
#
#     Rscript 02_dmr_binary_crops.R 500
#     Rscript 02_dmr_binary_crops.R 300 DA
#     Rscript 02_dmr_binary_crops.R 1000 DA,CNCC,IPSC
#
# PARAMETERS:
#   maxgap   - Max gap between CpGs (default: 500)
#   celltype - Cell type(s) to process (default: all)
#   
# CONFIGURABLE:
#   pen_value_range - CROPS search range (default: 0.01 to 2.0)
#   min_sites       - Min CpGs per DMR (default: 5)
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(changepoint)
})

CONFIG <- list(
  basedir = "",
  output_dir = "",
  default_celltypes = c("DA", "CNCC", "IPSC", "SKM", "HEP"),
  target_classes = c("pure_cis", "pure_trans", "cis_plus_trans", "cis_x_trans"),
  min_sites = 5,
  target_proportion_threshold = 0.5,
  min_effect_threshold = 0.05,
  conserved_effect_range = 0.05,
  
  # CROPS ptimization range
  pen_value_range = c(0.01, 2.0)
)

create_binary_scores <- function(data, target_class) {
  as.numeric(data$classification == target_class)
}

trim_binary_edges <- function(scores) {
  ones <- which(scores == 1)
  if (length(ones) == 0) return(seq_along(scores))
  min(ones):max(ones)
}

find_elbow_distance <- function(penalties, metric_values) {
  if (length(penalties) < 3) return(median(penalties, na.rm = TRUE))
  
  valid <- is.finite(penalties) & is.finite(metric_values)
  pen <- penalties[valid]
  met <- metric_values[valid]
  
  if (length(pen) < 3) return(median(penalties, na.rm = TRUE))
  
  ord <- order(pen)
  pen <- pen[ord]
  met <- met[ord]
  
  tryCatch({
    A <- met[length(met)] - met[1]
    B <- pen[1] - pen[length(pen)]
    C <- pen[length(pen)] * met[1] - pen[1] * met[length(met)]

    dists <- abs(A * pen + B * met + C) / sqrt(A^2 + B^2)
    pen[which.max(dists)]
  }, error = function(e) median(pen))
}
evaluate_crops <- function(scores, pen_range, min_seg) {
  tryCatch({
    crops <- cpt.mean(scores, method = "PELT", minseglen = min_seg,
                      penalty = "CROPS", pen.value = pen_range)
    
    pen_vals <- crops@pen.value.full
    cpts_mat <- crops@cpts.full
    
    if (is.null(pen_vals) || is.null(cpts_mat)) return(NULL)
    
    results <- lapply(seq_along(pen_vals), function(i) {
      cpts <- if (i <= nrow(cpts_mat)) {
        row <- cpts_mat[i, ]
        row[!is.na(row) & row > 0 & row <= length(scores)]
      } else numeric(0)
      
      data.frame(penalty = pen_vals[i], n_cpts = length(cpts))
    })
    
    do.call(rbind, results)
  }, error = function(e) NULL)
}
get_optimal_penalty_crops <- function(scores, pen_range, min_sites, target_threshold) {
  crops_perf <- evaluate_crops(scores, pen_range, min_sites)
  
  if (!is.null(crops_perf) && nrow(crops_perf) >= 3) {
    return(find_elbow_distance(crops_perf$penalty, crops_perf$n_cpts))
  }
  test_pens <- seq(pen_range[1], pen_range[2], length.out = 10)
  best <- test_pens[1]
  
  for (pen in test_pens) {
    tryCatch({
      cpt <- cpt.mean(scores, method = "PELT", minseglen = min_sites,
                      penalty = "Manual", pen.value = pen)
      cpts <- cpts(cpt)
      
      if (length(cpts) == 0) {
        if (mean(scores) > target_threshold) {
          best <- pen
          break
        }
      } else {
        breaks <- c(0, cpts, length(scores))
        valid <- 0
        for (j in 1:(length(breaks) - 1)) {
          seg <- scores[(breaks[j] + 1):breaks[j + 1]]
          if (length(seg) >= min_sites && mean(seg) > target_threshold) {
            valid <- valid + 1
          }
        }
        if (valid > 0) best <- pen
      }
    }, error = function(e) {})
  }
  
  best
}

args <- commandArgs(trailingOnly = TRUE)
maxgap <- if (length(args) >= 1) as.numeric(args[1]) else 500
celltypes <- if (length(args) >= 2) unlist(strsplit(args[2], ",")) else CONFIG$default_celltypes

cat("\n=== DMR Detection (Binary + CROPS) ===\n")
cat(sprintf("Max gap: %d | Cell types: %s\n", maxgap, paste(celltypes, collapse = ", ")))
cat(sprintf("Penalty range: [%.2f, %.2f]\n\n", CONFIG$pen_value_range[1], CONFIG$pen_value_range[2]))

detect_dmrs_binary_crops <- function(data, target_class, min_sites, max_gap,
                                     target_threshold, min_effect, conserved_range, pen_range) {
  
  scores <- create_binary_scores(data, target_class)
  
  ord <- order(data$chr, data$pos)
  data <- data[ord, ]
  scores <- scores[ord]
  
  results <- list()
  
  for (chrom in unique(data$chr)) {
    chr_idx <- data$chr == chrom
    chr_data <- data[chr_idx, ]
    chr_scores <- scores[chr_idx]
    
    if (length(chr_scores) < min_sites * 2) next
    
    # Split by gaps
    gaps <- diff(chr_data$pos)
    breaks <- which(gaps > max_gap)
    
    seg_starts <- c(1, breaks + 1)
    seg_ends <- c(breaks, nrow(chr_data))
    
    for (s in seq_along(seg_starts)) {
      si <- seg_starts[s]
      ei <- seg_ends[s]
      
      if (ei - si + 1 < min_sites * 2) next
      
      seg_scores <- chr_scores[si:ei]
      seg_data <- chr_data[si:ei, ]
      
      # Get optimal penalty for this segment
      opt_pen <- get_optimal_penalty_crops(seg_scores, pen_range, min_sites, target_threshold)
      
      # Changepoint detection
      tryCatch({
        cpt <- cpt.mean(seg_scores, method = "PELT", minseglen = min_sites,
                        penalty = "Manual", pen.value = opt_pen)
        cpts <- cpts(cpt)
        
        all_breaks <- c(0, cpts, length(seg_scores))
        
        for (i in 1:(length(all_breaks) - 1)) {
          sub_start <- all_breaks[i] + 1
          sub_end <- all_breaks[i + 1]
          
          if (sub_end - sub_start + 1 < min_sites) next
          
          sub_scores <- seg_scores[sub_start:sub_end]
          sub_data <- seg_data[sub_start:sub_end, ]
          
          target_prop <- mean(sub_scores)
          if (target_prop <= target_threshold) next
          
          # Edge trim
          trim_idx <- trim_binary_edges(sub_scores)
          if (length(trim_idx) < min_sites) next
          
          final_scores <- sub_scores[trim_idx]
          final_data <- sub_data[trim_idx, ]
          effects <- final_data$parent_HC_diff
          target_prop <- mean(final_scores)
          
          # Filters
          passes <- FALSE
          dom_dir <- NA
          dir_cons <- NA
          
          if (target_class == "conserved") {
            small <- abs(effects) <= conserved_range
            prop_small <- mean(small)
            if (prop_small >= target_threshold) {
              passes <- TRUE
              dom_dir <- "stable"
              dir_cons <- prop_small
            }
          } else {
            strong <- effects[abs(effects) > min_effect]
            if (length(strong) > 0) {
              pos_n <- sum(strong > 0)
              neg_n <- sum(strong < 0)
              dir_cons <- max(pos_n, neg_n) / length(strong)
              dom_dir <- ifelse(pos_n > neg_n, "hyper", "hypo")
              passes <- TRUE
            }
          }
          
          if (passes) {
            results[[length(results) + 1]] <- data.frame(
              chr = paste0("chr", chrom),
              start = min(final_data$pos),
              end = max(final_data$pos),
              n_cpgs = nrow(final_data),
              target_class = target_class,
              target_proportion = target_prop,
              mean_effect = mean(abs(effects)),
              median_effect = median(effects),
              effect_sd = sd(effects),
              dominant_direction = dom_dir,
              directional_consistency = dir_cons,
              optimal_penalty = opt_pen,
              stringsAsFactors = FALSE
            )
          }
        }
      }, error = function(e) {})
    }
  }
  
  if (length(results) > 0) {
    combined <- do.call(rbind, results)
    if (target_class == "conserved") {
      combined <- combined[order(-combined$target_proportion, combined$effect_sd), ]
    } else {
      combined <- combined[order(-combined$target_proportion, -combined$mean_effect), ]
    }
    return(combined)
  }
  data.frame()
}

if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

dmr_list <- list()

for (celltype in celltypes) {
  cat(celltype)
  
  data_file <- sprintf("%scistrans_results_%s_5bp_coverage/01_individual_cpg_analysis/%s_individual_cpg_results.csv",
                       CONFIG$basedir, celltype, celltype)
  
  data <- read.csv(data_file)
  
  n_cores <- min(detectCores() - 1, length(CONFIG$target_classes), 4)
  n_cores <- max(n_cores, 1)
  
  
  all_dmrs <- if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    clusterEvalQ(cl, { library(changepoint); library(dplyr) })
    clusterExport(cl, c("detect_dmrs_binary_crops", "create_binary_scores", 
                        "trim_binary_edges", "get_optimal_penalty_crops",
                        "evaluate_crops", "find_elbow_distance",
                        "data", "maxgap", "CONFIG"), envir = environment())
    
    result <- tryCatch({
      foreach(target = CONFIG$target_classes, .combine = list, .multicombine = TRUE) %dopar% {
        dmrs <- detect_dmrs_binary_crops(data, target, CONFIG$min_sites, maxgap,
                                         CONFIG$target_proportion_threshold,
                                         CONFIG$min_effect_threshold,
                                         CONFIG$conserved_effect_range,
                                         CONFIG$pen_value_range)
        list(class = target, dmrs = dmrs)
      }
    }, finally = stopCluster(cl))
    result
  } else {
    lapply(CONFIG$target_classes, function(target) {
      cat(sprintf("  %s...\n", target))
      dmrs <- detect_dmrs_binary_crops(data, target, CONFIG$min_sites, maxgap,
                                       CONFIG$target_proportion_threshold,
                                       CONFIG$min_effect_threshold,
                                       CONFIG$conserved_effect_range,
                                       CONFIG$pen_value_range)
      list(class = target, dmrs = dmrs)
    })
  }
  
  dmr_results <- list()
  for (r in all_dmrs) {
    if (!is.null(r) && nrow(r$dmrs) > 0) {
      dmr_results[[r$class]] <- r$dmrs
      cat(sprintf("  %s: %d DMRs\n", r$class, nrow(r$dmrs)))
    }
  }
  
  if (length(dmr_results) > 0) {
    combined <- do.call(rbind, dmr_results)
    combined <- combined[order(-abs(combined$mean_effect)), ]
    
    cat(sprintf("\nTotal: %d DMRs\n", nrow(combined)))
    print(table(combined$target_class))
    
    csv_file <- sprintf("%s%s_dmrs_binary_%d.csv", CONFIG$output_dir, celltype, maxgap)
    write.csv(combined, csv_file, row.names = FALSE)
    
    gr <- GRanges(
      seqnames = combined$chr,
      ranges = IRanges(start = combined$start, end = combined$end),
      target_class = combined$target_class,
      n_cpgs = combined$n_cpgs,
      mean_effect = combined$mean_effect,
      dominant_direction = combined$dominant_direction,
      optimal_penalty = combined$optimal_penalty
    )
    saveRDS(gr, sprintf("%s%s_dmrs_binary_%d.rds", CONFIG$output_dir, celltype, maxgap))
    
    dmr_list[[celltype]] <- combined
  }
}

saveRDS(dmr_list, sprintf("%sdmr_list_binary_%d.rds", CONFIG$output_dir, maxgap))
