#===============================================================================
# 04_main_multibed_pipeline.R - for classifying regional methylation data for predefined regions
#===============================================================================
# DESCRIPTION:
#   Integrated pipeline for individual CpG analysis and regional analysis across
#   multiple BED files. Orchestrates data loading, analysis, and output generation.
#   This code will run individual CpG classification first in order to compute heterogeneity metrics
#
# INPUT FORMATS:
#   methylation_files - Vector of paths to DSS-format methylation files
#   sample_names      - Vector of sample names (must contain parH, parC, hybH, hybC)
#   bed_files         - Named vector/list of BED file paths
#
# ADJUSTABLE PARAMETERS:
#   alpha                    - Significance threshold (default: 0.05)
#   min_coverage_per_sample  - Minimum coverage per sample (default: 5)
#   min_samples_with_coverage - Samples required with coverage (default: all)
#   min_cpgs_per_region      - Minimum CpGs per region (default: 3)
#   min_total_coverage       - Minimum total coverage per region (default: 50)
#   save_all_intermediates   - Save RData files for replotting (default: TRUE)
#   create_all_plots         - Generate all visualization plots (default: TRUE)
#   calculate_confidence     - Calculate site-level support metrics (default: TRUE)
#
# OUTPUTS:
#   Directory structure:
#     main_output_dir/
#     ├── 01_individual_cpg_analysis/
#     │   ├── individual_cpg_results.csv      (per-CpG classifications)
#     │   ├── individual_cpg_scatter_plot.png
#     │   ├── individual_cpg_summary_plot.png
#     │   └── individual_analysis_complete.RData
#     ├── 02_regional_[bed_name]/             (one per BED file)
#     │   ├── regional_results.csv            (per-region classifications)
#     │   ├── confidence_metrics.csv
#     │   ├── heterogeneity_analysis.csv
#     │   └── plots/
#     ├── 03_cross_bed_comparison/            (if multiple BED files)
#     ├── cross_bed_summary.csv
#     └── report.txt
#
# DEPENDENCIES:
#   Requires: 05_CIS_TRANS_CORE.R, 01_data_processing.R, 02_regional_analysis.R,
#             03_visualization.R
#===============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

#-------------------------------------------------------------------------------
# MAIN PIPELINE FUNCTION
#-------------------------------------------------------------------------------

RunCompleteMultiBEDPipeline <- function(
    methylation_files,
    sample_names,
    bed_files,
    main_output_dir = "complete_analysis",
    alpha = 0.05,
    min_coverage_per_sample = 5,
    min_samples_with_coverage = NULL,
    min_cpgs_per_region = 3,
    min_total_coverage = 50,
    save_all_intermediates = TRUE,
    create_all_plots = TRUE,
    calculate_confidence = TRUE,
    verbose = TRUE
) {
  
  if(verbose) cat("Starting cis-trans analysis pipeline at", as.character(Sys.time()), "\n")
  
  # Input validation
  if(length(methylation_files) != length(sample_names)) {
    stop("Number of methylation files must match number of sample names")
  }
  if(length(bed_files) == 0) stop("At least one BED file required")
  
  missing_meth <- methylation_files[!file.exists(methylation_files)]
  if(length(missing_meth) > 0) stop("Missing files: ", paste(missing_meth, collapse = ", "))
  
  missing_bed <- bed_files[!file.exists(bed_files)]
  if(length(missing_bed) > 0) stop("Missing BED files: ", paste(missing_bed, collapse = ", "))
  
  if(!dir.exists(main_output_dir)) dir.create(main_output_dir, recursive = TRUE)
  
  pipeline_params <- list(
    methylation_files = methylation_files,
    sample_names = sample_names,
    bed_files = bed_files,
    alpha = alpha,
    min_coverage_per_sample = min_coverage_per_sample,
    min_cpgs_per_region = min_cpgs_per_region,
    min_total_coverage = min_total_coverage,
    analysis_start_time = Sys.time()
  )
  
  if(save_all_intermediates) {
    save(pipeline_params, file = file.path(main_output_dir, "pipeline_parameters.RData"))
  }
  
  #--- STEP 1: Load and prepare data ---
  if(verbose) cat("Loading methylation data...\n")
  BSobj <- LoadAndCreateBSseq(methylation_files, sample_names)
  
  if(is.null(min_samples_with_coverage)) {
    min_samples_with_coverage <- length(sample_names)
  }
  
  BSobj_filtered <- FilterBSseqObject(BSobj, min_coverage_per_sample, 
                                      min_samples_with_coverage)
  
  if(save_all_intermediates) {
    save(BSobj, BSobj_filtered, file = file.path(main_output_dir, "BSseq_objects.RData"))
  }
  
  if(verbose) cat("Filtered:", nrow(BSobj_filtered), "CpGs\n")
  
  #--- STEP 2: Individual CpG analysis ---
  if(verbose) cat("Running individual CpG analysis...\n")
  
  individual_output_dir <- file.path(main_output_dir, "01_individual_cpg_analysis")
  if(!dir.exists(individual_output_dir)) dir.create(individual_output_dir, recursive = TRUE)
  
  individual_results <- TwoModelCisTransAnalysis(BSobj_filtered, alpha, 
                                                 min_coverage_per_sample, verbose = FALSE)
  
  write.csv(individual_results$cpg_results, 
            file.path(individual_output_dir, "individual_cpg_results.csv"), 
            row.names = FALSE)
  
  if(save_all_intermediates) {
    save(individual_results, 
         file = file.path(individual_output_dir, "individual_analysis_complete.RData"))
  }
  
  if(create_all_plots) {
    individual_scatter <- PlotTwoModelResults(individual_results)
    individual_summary <- PlotTwoModelSummary(individual_results)
    individual_log2fc <- PlotTwoModelResultsLog2FC(individual_results)
    
    if(!is.null(individual_scatter$plot)) {
      ggsave(file.path(individual_output_dir, "individual_cpg_scatter_plot.png"), 
             individual_scatter$plot, width = 12, height = 8, dpi = 300)
    }
    if(!is.null(individual_log2fc$plot)) {
      ggsave(file.path(individual_output_dir, "individual_cpg_log2fc_plot.png"), 
             individual_log2fc$plot, width = 12, height = 8, dpi = 300)
    }
    if(!is.null(individual_summary)) {
      ggsave(file.path(individual_output_dir, "individual_cpg_summary_plot.png"), 
             individual_summary, width = 10, height = 6, dpi = 300)
    }
  }
  
  #--- STEP 3: Regional analysis for each BED file ---
  if(verbose) cat("Processing BED files...\n")
  
  if(is.null(names(bed_files))) {
    bed_names <- paste0("bedfile_", seq_along(bed_files))
    names(bed_files) <- bed_names
  } else {
    bed_names <- names(bed_files)
  }
  
  all_regional_results <- list()
  regional_summary_data <- data.frame()
  
  for(i in seq_along(bed_files)) {
    bed_name <- bed_names[i]
    bed_file_path <- bed_files[i]
    
    if(verbose) cat("  Processing:", bed_name, "\n")
    
    bed_output_dir <- file.path(main_output_dir, sprintf("02_regional_%s", bed_name))
    if(!dir.exists(bed_output_dir)) dir.create(bed_output_dir, recursive = TRUE)
    
    tryCatch({
      regions_gr <- LoadRegionsFromBED(bed_file_path, region_name = bed_name)
      
      regional_results <- RegionalCisTransAnalysis(
        BSobj_filtered, regions_gr, min_cpgs_per_region, 
        min_total_coverage, alpha, analyze_heterogeneity = TRUE
      )
      
      confidence_results <- NULL
      if(calculate_confidence) {
        confidence_results <- CalculateRegionalConfidence(
          BSobj_filtered, regions_gr, regional_results, individual_results
        )
        regional_results$region_results <- merge(
          regional_results$region_results,
          confidence_results[, c("region_id", "n_supporting_sites", "n_opposing_sites", 
                                 "support_fraction", "confidence_level")],
          by.x = "region_id", by.y = "region_id", all.x = TRUE
        )
      }
      
      write.csv(regional_results$region_results, 
                file.path(bed_output_dir, "regional_results.csv"), row.names = FALSE)
      
      if(!is.null(confidence_results)) {
        write.csv(confidence_results, 
                  file.path(bed_output_dir, "confidence_metrics.csv"), row.names = FALSE)
      }
      
      if(!is.null(regional_results$heterogeneity_results)) {
        write.csv(regional_results$heterogeneity_results,
                  file.path(bed_output_dir, "heterogeneity_analysis.csv"), row.names = FALSE)
      }
      
      if(save_all_intermediates) {
        save(regional_results, confidence_results, regions_gr,
             file = file.path(bed_output_dir, "regional_analysis_complete.RData"))
      }
      
      if(create_all_plots) {
        plots_dir <- file.path(bed_output_dir, "plots")
        CreateAllRegionalPlots(regional_results, confidence_results, 
                               BSobj_filtered, individual_results, plots_dir)
      }
      
      all_regional_results[[bed_name]] <- list(
        bed_file = bed_file_path,
        regional_results = regional_results,
        confidence_results = confidence_results,
        regions_gr = regions_gr,
        status = "success"
      )
      
      summary_row <- data.frame(
        bed_name = bed_name,
        status = "success",
        total_regions = regional_results$summary_stats$total_regions,
        meeting_requirements = regional_results$summary_stats$regions_meeting_requirements,
        classified = regional_results$summary_stats$classified_regions,
        success_rate = round(100 * regional_results$summary_stats$classified_regions / 
                               regional_results$summary_stats$total_regions, 1),
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      all_regional_results[[bed_name]] <- list(status = "failed", error = e$message)
      summary_row <- data.frame(
        bed_name = bed_name, status = "failed",
        total_regions = NA, meeting_requirements = NA, 
        classified = NA, success_rate = NA,
        stringsAsFactors = FALSE
      )
    })
    
    regional_summary_data <- rbind(regional_summary_data, summary_row)
  }
  
  #--- STEP 4: Cross-BED summary ---
  write.csv(regional_summary_data, 
            file.path(main_output_dir, "cross_bed_summary.csv"), row.names = FALSE)
  
  successful_beds <- regional_summary_data[regional_summary_data$status == "success", ]
  
  if(nrow(successful_beds) > 1 && create_all_plots) {
    comparison_dir <- file.path(main_output_dir, "03_cross_bed_comparison")
    if(!dir.exists(comparison_dir)) dir.create(comparison_dir)
    
    comparison_plots <- CreateCrossBEDComparisonPlots(successful_beds, all_regional_results)
    if(!is.null(comparison_plots)) {
      for(plot_name in names(comparison_plots)) {
        if(!is.null(comparison_plots[[plot_name]])) {
          ggsave(file.path(comparison_dir, paste0(plot_name, ".png")), 
                 comparison_plots[[plot_name]], width = 12, height = 8, dpi = 300)
        }
      }
    }
  }
  
  #--- STEP 5: Generate reports ---
  GenerateCompleteAnalysisReport(
    individual_results, all_regional_results, regional_summary_data, pipeline_params,
    file.path(main_output_dir, "COMPLETE_ANALYSIS_REPORT.txt")
  )
    
  runtime_minutes <- as.numeric(difftime(Sys.time(), pipeline_params$analysis_start_time, units = "mins"))
  
  if(verbose) {
    cat("\nAnalysis completed in", round(runtime_minutes, 1), "minutes\n")
    cat("Results saved to:", main_output_dir, "\n")
  }
  
  return(list(
    pipeline_params = pipeline_params,
    individual_results = individual_results,
    all_regional_results = all_regional_results,
    regional_summary_data = regional_summary_data,
    BSobj_filtered = BSobj_filtered,
    main_output_dir = main_output_dir,
    runtime_minutes = runtime_minutes
  ))
}

#-------------------------------------------------------------------------------
# CROSS-BED COMPARISON PLOTS
#-------------------------------------------------------------------------------

CreateCrossBEDComparisonPlots <- function(successful_beds, all_regional_results) {
  # Create comparison plots across BED files
  if(nrow(successful_beds) < 2) return(NULL)
  
  plots <- list()
  
  success_plot_data <- successful_beds[order(successful_beds$success_rate, decreasing = TRUE), ]
  plots$success_rates <- ggplot(success_plot_data, 
                                aes(x = reorder(bed_name, success_rate), y = success_rate)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste0(success_rate, "%")), hjust = -0.1) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Classification Success Rates", x = "BED File", y = "Success Rate (%)") +
    ylim(0, max(success_plot_data$success_rate) * 1.1)
  
  # Classification comparison
  classification_data <- data.frame()
  for(bed_name in successful_beds$bed_name) {
    result <- all_regional_results[[bed_name]]
    if(!is.null(result$regional_results)) {
      classified <- result$regional_results$region_results[
        result$regional_results$region_results$meets_requirements & 
          result$regional_results$region_results$classification != "unclassified", ]
      if(nrow(classified) > 0) {
        class_counts <- table(classified$classification)
        for(class in names(class_counts)) {
          classification_data <- rbind(classification_data,
                                       data.frame(bed_name = bed_name, classification = class, 
                                                  count = as.numeric(class_counts[class])))
        }
      }
    }
  }
  
  if(nrow(classification_data) > 0) {
    plots$classification_comparison <- ggplot(classification_data, 
                                              aes(x = bed_name, y = count, fill = classification)) +
      geom_bar(stat = "identity", position = "fill") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Classification Distribution", x = "BED File", y = "Proportion") +
      scale_fill_brewer(type = "qual", palette = "Set3")
  }
  
  return(plots)
}

#-------------------------------------------------------------------------------
# REPORT GENERATION
#-------------------------------------------------------------------------------

GenerateCompleteAnalysisReport <- function(individual_results, all_regional_results, 
                                           regional_summary_data, pipeline_params, output_file) {
  # Generate comprehensive text report
  sink(output_file)
  
  cat("COMPLETE CIS-TRANS METHYLATION ANALYSIS REPORT\n")
  cat(strrep("=", 60), "\n\n")
  cat("Analysis completed:", as.character(Sys.time()), "\n")
  cat("Runtime:", round(as.numeric(difftime(Sys.time(), pipeline_params$analysis_start_time, 
                                            units = "mins")), 1), "minutes\n\n")
  
  cat("PARAMETERS\n")
  cat(strrep("-", 40), "\n")
  cat("Alpha:", pipeline_params$alpha, "\n")
  cat("Min coverage per sample:", pipeline_params$min_coverage_per_sample, "\n")
  cat("Min CpGs per region:", pipeline_params$min_cpgs_per_region, "\n")
  cat("Min total coverage:", pipeline_params$min_total_coverage, "\n\n")
  
  cat("INDIVIDUAL CpG RESULTS\n")
  cat(strrep("-", 40), "\n")
  cat("Total CpGs:", nrow(individual_results$cpg_results), "\n")
  for(class in names(individual_results$summary)) {
    count <- individual_results$summary[class]
    percent <- round(100 * count / sum(individual_results$summary), 1)
    cat("  ", class, ":", count, "(", percent, "%)\n")
  }
  
  cat("\nREGIONAL RESULTS\n")
  cat(strrep("-", 40), "\n")
  for(bed_name in names(all_regional_results)) {
    result <- all_regional_results[[bed_name]]
    cat("\n", bed_name, ":\n")
    if(result$status == "success") {
      stats <- result$regional_results$summary_stats
      cat("  Total regions:", stats$total_regions, "\n")
      cat("  Classified:", stats$classified_regions, "\n")
    } else {
      cat("  Status: FAILED\n")
    }
  }
  
  sink()
}