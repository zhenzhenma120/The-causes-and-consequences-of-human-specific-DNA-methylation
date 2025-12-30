#===============================================================================
# 03_visualization.R
#===============================================================================
# DESCRIPTION:
#   Plotting functions for visualizing individual CpG and regional cis-trans
#   analysis results.
#
# INPUT FORMATS:
#   results object - Output from TwoModelCisTransAnalysis()
#   regional_results object - Output from RegionalCisTransAnalysis()
#   confidence_results - Output from CalculateRegionalConfidence()
#
# ADJUSTABLE PARAMETERS:
#   color_by_confidence - Use shapes for confidence levels (default: FALSE)
#   show_unclassified   - Include unclassified regions (default: FALSE)
#   pseudocount         - For log2FC calculations (default: 0.01)
#
# OUTPUTS:
#     - Scatter plots (H-C difference and log2FC views)
#     - Summary bar plots
#     - Heterogeneity plots
#     - Confidence metrics plots
#     - Genomic distribution plots
#===============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(GenomicRanges)
})

# Classification color palette (consistent across all plots)
CLASSIFICATION_COLORS <- c(
  "conserved" = "#999999",
  "pure_cis" = "#1F78B4",
  "pure_trans" = "#33A02C",
  "cis_plus_trans" = "#FF7F00",
  "cis_x_trans" = "#E31A1C",
  "unclassified" = "#DDDDDD"
)

CONFIDENCE_SHAPES <- c(
  "high_confidence" = 16, "medium_confidence" = 17, "low_confidence" = 15,
  "no_classified_sites" = 4, "region_unclassified" = 3, "no_sites" = 2
)

#-------------------------------------------------------------------------------
# INDIVIDUAL CpG PLOTS
#-------------------------------------------------------------------------------

PlotTwoModelResults <- function(results, title = "Two-Model Cis-Trans Analysis") {
  # Scatter plot of individual CpG results
  plot_data <- results$cpg_results[results$cpg_results$classification != "unclassified", ]
  n_unclassified <- sum(results$cpg_results$classification == "unclassified")
  
  if(nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot(plot_data, aes(x = parent_HC_diff, y = hybrid_HC_diff, color = classification)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = CLASSIFICATION_COLORS, name = "Regulation Type") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right") +
    labs(title = title,
         subtitle = paste0(n_unclassified, " unclassified sites excluded"),
         x = "H-C Difference in Parents", y = "H-C Difference in Hybrids",
         caption = "Diagonal = cis-only expectation") +
    coord_fixed(ratio = 1)
  
  return(list(plot = p, data = plot_data, n_unclassified = n_unclassified))
}

PlotTwoModelResultsLog2FC <- function(results, title = "Two-Model Analysis (Log2FC)", 
                                      pseudocount = 0.01) {
  # Log2FC scatter plot
  plot_data <- results$cpg_results[results$cpg_results$classification != "unclassified", ]
  n_unclassified <- sum(results$cpg_results$classification == "unclassified")
  
  if(nrow(plot_data) == 0) return(NULL)
  
  baseline_meth <- 0.5
  plot_data$parent_H_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth + plot_data$parent_HC_diff/2))
  plot_data$parent_C_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth - plot_data$parent_HC_diff/2))
  plot_data$hybrid_H_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth + plot_data$hybrid_HC_diff/2))
  plot_data$hybrid_C_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth - plot_data$hybrid_HC_diff/2))
  
  plot_data$parent_log2FC <- log2(plot_data$parent_H_est / plot_data$parent_C_est)
  plot_data$hybrid_log2FC <- log2(plot_data$hybrid_H_est / plot_data$hybrid_C_est)
  
  p <- ggplot(plot_data, aes(x = parent_log2FC, y = hybrid_log2FC, color = classification)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = CLASSIFICATION_COLORS, name = "Regulation Type") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right") +
    labs(title = title,
         subtitle = paste0(n_unclassified, " unclassified sites excluded"),
         x = "Log2(H/C) in Parents", y = "Log2(H/C) in Hybrids",
         caption = "Diagonal = cis-only expectation") +
    coord_fixed(ratio = 1)
  
  return(list(plot = p, data = plot_data, n_unclassified = n_unclassified))
}

PlotTwoModelSummary <- function(results) {
  # Summary bar plot of classifications
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

#-------------------------------------------------------------------------------
# REGIONAL PLOTS
#-------------------------------------------------------------------------------

PlotRegionalClassification <- function(regional_results, 
                                       title = "Regional Cis-Trans Classification",
                                       show_unclassified = FALSE) {
  # Regional scatter plot
  plot_data <- regional_results$region_results
  if(!show_unclassified) {
    plot_data <- plot_data[plot_data$meets_requirements & 
                             plot_data$classification != "unclassified", ]
  }
  if(nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot(plot_data, aes(x = parent_HC_diff, y = hybrid_HC_diff)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.7) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right") +
    labs(title = title, x = "Parent H-C Difference (Regional)", 
         y = "Hybrid H-C Difference (Regional)",
         caption = "Diagonal = cis-only expectation") +
    coord_fixed(ratio = 1)
  
  if("confidence_level" %in% colnames(plot_data)) {
    p <- p + geom_point(aes(color = classification, shape = confidence_level), 
                        size = 3, alpha = 0.8) +
      scale_color_manual(values = CLASSIFICATION_COLORS, name = "Classification") +
      scale_shape_manual(values = CONFIDENCE_SHAPES, name = "Confidence")
  } else {
    p <- p + geom_point(aes(color = classification), size = 3, alpha = 0.8) +
      scale_color_manual(values = CLASSIFICATION_COLORS, name = "Classification")
  }
  
  return(p)
}

PlotRegionalClassificationLog2FC <- function(regional_results, 
                                             title = "Regional Classification (Log2FC)",
                                             pseudocount = 0.01) {
  # Regional log2FC scatter plot
  plot_data <- regional_results$region_results
  plot_data <- plot_data[plot_data$meets_requirements & 
                           plot_data$classification != "unclassified", ]
  if(nrow(plot_data) == 0) return(NULL)
  
  baseline_meth <- 0.5
  plot_data$parent_H_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth + plot_data$parent_HC_diff/2))
  plot_data$parent_C_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth - plot_data$parent_HC_diff/2))
  plot_data$hybrid_H_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth + plot_data$hybrid_HC_diff/2))
  plot_data$hybrid_C_est <- pmax(pseudocount, pmin(1-pseudocount, 
                                                    baseline_meth - plot_data$hybrid_HC_diff/2))
  
  plot_data$parent_log2FC <- log2(plot_data$parent_H_est / plot_data$parent_C_est)
  plot_data$hybrid_log2FC <- log2(plot_data$hybrid_H_est / plot_data$hybrid_C_est)
  
  p <- ggplot(plot_data, aes(x = parent_log2FC, y = hybrid_log2FC)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(title = title, x = "Log2(H/C) in Parents", y = "Log2(H/C) in Hybrids") +
    coord_fixed(ratio = 1)
  
  if("confidence_level" %in% colnames(plot_data)) {
    p <- p + geom_point(aes(color = classification, shape = confidence_level), 
                        size = 3, alpha = 0.8) +
      scale_color_manual(values = CLASSIFICATION_COLORS, name = "Classification") +
      scale_shape_manual(values = CONFIDENCE_SHAPES, name = "Confidence")
  } else {
    p <- p + geom_point(aes(color = classification), size = 3, alpha = 0.8) +
      scale_color_manual(values = CLASSIFICATION_COLORS, name = "Classification")
  }
  
  return(p)
}

PlotRegionalSummary <- function(regional_results, by_confidence = FALSE) {
  # Regional summary bar plot
  classified_data <- regional_results$region_results[
    regional_results$region_results$meets_requirements & 
      regional_results$region_results$classification != "unclassified", ]
  if(nrow(classified_data) == 0) return(NULL)
  
  if(by_confidence && "confidence_level" %in% colnames(classified_data)) {
    summary_data <- as.data.frame(table(classified_data$classification, 
                                        classified_data$confidence_level))
    colnames(summary_data) <- c("Classification", "Confidence", "Count")
    
    p <- ggplot(summary_data, aes(x = Classification, y = Count, fill = Confidence)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Regional Classification Summary",
           subtitle = paste0(nrow(classified_data), " classified regions"),
           x = "Regulatory Pattern", y = "Number of Regions") +
      scale_fill_brewer(type = "qual", palette = "Set2")
  } else {
    summary_counts <- table(classified_data$classification)
    plot_data <- data.frame(
      Classification = names(summary_counts),
      Count = as.numeric(summary_counts)
    )
    
    p <- ggplot(plot_data, aes(x = Classification, y = Count, fill = Classification)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(title = "Regional Classification Summary",
           subtitle = paste0(nrow(classified_data), " classified regions"),
           x = "Regulatory Pattern", y = "Number of Regions") +
      scale_fill_brewer(type = "qual", palette = "Set3")
  }
  
  return(p)
}

#-------------------------------------------------------------------------------
# HETEROGENEITY AND CONFIDENCE PLOTS
#-------------------------------------------------------------------------------

PlotRegionalHeterogeneity <- function(regional_results, metric = "cv") {
  # Plot regional heterogeneity metrics
  if(is.null(regional_results$heterogeneity_results)) return(NULL)
  
  het_data <- regional_results$heterogeneity_results
  region_data <- regional_results$region_results
  plot_data <- merge(het_data, region_data, by.x = "region_index", by.y = "region_id")
  plot_data <- plot_data[plot_data$meets_requirements & 
                           plot_data$classification != "unclassified" &
                           plot_data$n_cpgs_in_region > 1, ]
  
  if(nrow(plot_data) == 0) return(NULL)
  
  if(metric == "cv") {
    valid_cv <- !is.na(plot_data$parent_diff_cv) & is.finite(plot_data$parent_diff_cv)
    if(sum(valid_cv) == 0) return(NULL)
    
    p <- ggplot(plot_data[valid_cv, ], aes(x = classification, y = parent_diff_cv, 
                                           fill = classification)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(title = "Regional Heterogeneity", x = "Classification", 
           y = "Coefficient of Variation") +
      scale_fill_brewer(type = "qual", palette = "Set3")
    return(p)
  }
  
  if(metric == "consistency" && "pattern_consistency" %in% colnames(plot_data)) {
    consistency_summary <- as.data.frame(table(plot_data$classification, 
                                               plot_data$pattern_consistency))
    colnames(consistency_summary) <- c("Classification", "Consistency", "Count")
    
    p <- ggplot(consistency_summary, aes(x = Classification, y = Count, fill = Consistency)) +
      geom_bar(stat = "identity", position = "fill") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Pattern Consistency by Classification", x = "Classification", 
           y = "Proportion") +
      scale_fill_brewer(type = "qual", palette = "Set2")
    return(p)
  }
  
  return(NULL)
}

PlotConfidenceMetrics <- function(regional_results, confidence_results = NULL) {
  # Plot confidence metrics
  if(is.null(confidence_results)) {
    if("support_fraction" %in% colnames(regional_results$region_results)) {
      plot_data <- regional_results$region_results
    } else return(NULL)
  } else {
    plot_data <- merge(regional_results$region_results, confidence_results, by = "region_id")
  }
  
  plot_data <- plot_data[plot_data$meets_requirements & 
                           plot_data$classification != "unclassified", ]
  if(nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot(plot_data, aes(x = classification, y = support_fraction, fill = classification)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_hline(yintercept = c(0.6, 0.8), linetype = "dashed", color = "red", alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "Site-Level Support for Regional Classifications",
         x = "Classification", y = "Fraction of Supporting Sites") +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    ylim(0, 1)
  
  return(p)
}

#-------------------------------------------------------------------------------
# GENOMIC DISTRIBUTION PLOTS
#-------------------------------------------------------------------------------

PlotGenomicDistribution <- function(regional_results, plot_type = "chromosome") {
  # Plot genomic distribution of classifications
  plot_data <- regional_results$region_results[
    regional_results$region_results$meets_requirements & 
      regional_results$region_results$classification != "unclassified", ]
  if(nrow(plot_data) == 0) return(NULL)
  
  if(plot_type == "chromosome") {
    chr_summary <- as.data.frame(table(plot_data$chr, plot_data$classification))
    colnames(chr_summary) <- c("Chromosome", "Classification", "Count")
    
    p <- ggplot(chr_summary, aes(x = Chromosome, y = Count, fill = Classification)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Classifications by Chromosome", x = "Chromosome", 
           y = "Number of Regions") +
      scale_fill_brewer(type = "qual", palette = "Set3")
    return(p)
  }
  
  return(NULL)
}

#-------------------------------------------------------------------------------
# COMPREHENSIVE PLOT CREATION
#-------------------------------------------------------------------------------

CreateAllRegionalPlots <- function(regional_results, confidence_results = NULL,
                                   BSobj = NULL, site_level_results = NULL,
                                   output_dir = NULL) {
  # Create all regional analysis plots
  if(!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots <- list()
  
  tryCatch({plots$classification <- PlotRegionalClassification(regional_results)}, 
           error = function(e) NULL)
  tryCatch({plots$summary <- PlotRegionalSummary(regional_results)}, 
           error = function(e) NULL)
  tryCatch({plots$summary_confidence <- PlotRegionalSummary(regional_results, by_confidence = TRUE)}, 
           error = function(e) NULL)
  tryCatch({plots$heterogeneity_cv <- PlotRegionalHeterogeneity(regional_results, metric = "cv")}, 
           error = function(e) NULL)
  tryCatch({plots$heterogeneity_consistency <- PlotRegionalHeterogeneity(regional_results, metric = "consistency")}, 
           error = function(e) NULL)
  tryCatch({plots$confidence <- PlotConfidenceMetrics(regional_results, confidence_results)}, 
           error = function(e) NULL)
  tryCatch({plots$chr_distribution <- PlotGenomicDistribution(regional_results, plot_type = "chromosome")}, 
           error = function(e) NULL)
  tryCatch({plots$classification_log2fc <- PlotRegionalClassificationLog2FC(regional_results)}, 
           error = function(e) NULL)
  
  if(!is.null(output_dir)) {
    SaveRegionalPlots(plots, output_dir)
  }
  
  plots <- plots[!sapply(plots, is.null)]
  return(plots)
}

SaveRegionalPlots <- function(plots, output_dir) {
  # Save all plots to files
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for(plot_name in names(plots)) {
    if(!is.null(plots[[plot_name]])) {
      filename <- file.path(output_dir, paste0("regional_", plot_name, ".png"))
      tryCatch({
        ggsave(filename, plots[[plot_name]], width = 12, height = 8, dpi = 300)
      }, error = function(e) NULL)
    }
  }
}