#!/usr/bin/env Rscript
#===============================================================================
# 06_RUN_ANALYSIS.R
#===============================================================================
# DESCRIPTION:
#   Runs the complete cis-trans methylation analysis pipeline for a specific
#   cell type. Reads metadata to identify samples and processes all BED files.
#
# USAGE:
#   Rscript 06_RUN_ANALYSIS.R <CELLTYPE>
#   Example: Rscript 06_RUN_ANALYSIS.R DA
#
# INPUT FORMATS:
#   Metadata CSV file with columns:
#     - cell.type: Cell type identifier (DA, CNCC, IPSC, SKM, HEP)
#     - filename: Original BAM filename (converted to DSS filename)
#     - species: H or C (for Human/Chimp allele)
#     - system: Par or Hyb (for Parental/Hybrid)
#     - replicate: Replicate number
#
#   DSS methylation files (tab-separated):
#     chr    pos    N    X
#
#   BED files (standard BED format):
#     chr    start    end    [name]
#
# ADJUSTABLE PARAMETERS (modify in STEP 6 below):
#   alpha                    - Significance threshold (default: 0.05)
#   min_coverage_per_sample  - Minimum coverage per sample (default: 5)
#   min_cpgs_per_region      - Minimum CpGs per region (default: 3)
#   min_total_coverage       - Minimum total coverage per region (default: 15)
#
# OUTPUTS:
#   Directory: cistrans_results_<CELLTYPE>/
#     - Individual CpG analysis results and plots
#     - Regional analysis results for each BED file
#     - Summary reports and navigation index
#
# CELL TYPES used in this paper:
#   DA, CNCC, IPSC, SKM, HEP
#===============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(DSS)
})

#===============================================================================
# STEP 1: SET WORKING DIRECTORY AND LOAD MODULES
#===============================================================================

#setwd('')

source("05_CIS_TRANS_CORE.R")
source("01_data_processing.R")
source("02_regional_analysis.R")
source("03_visualization.R")
source("04_main_multibed_pipeline.R")

#===============================================================================
# STEP 2: GET CELL TYPE TO RUN
#===============================================================================

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  stop("Usage: Rscript 06_RUN_ANALYSIS.R <CELLTYPE>\nValid types: DA, CNCC, IPSC, SKM, HEP")
}

CELLTYPE <- args[1]
valid_celltypes <- c("DA", "CNCC", "IPSC", "SKM", "HEP")
if(!CELLTYPE %in% valid_celltypes) {
  stop("Invalid cell type. Must be: ", paste(valid_celltypes, collapse = ", "))
}

cat("Analyzing cell type:", CELLTYPE, "\n")

#===============================================================================
# STEP 3: READ METADATA AND BUILD FILE PATHS
# MODIFY METADATA HERE
#===============================================================================

metadata_file <- ""
cat("metadata file: ", metadata_file)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
celltype_metadata <- metadata[metadata$cell.type == CELLTYPE, ]

if(nrow(celltype_metadata) == 0) stop("No samples found for: ", CELLTYPE)

# Build DSS file paths
#celltype_metadata$dss_filename <- gsub("\\.cov\\.gz$", "_DSS.txt", celltype_metadata$filename)
celltype_metadata$dss_filepath <- file.path(celltype_metadata$dss_filename)

# Create sample names (parH, parC, hybH, hybC pattern)
celltype_metadata$sample_name <- ifelse(
  celltype_metadata$system == "Hyb",
  paste0("hyb", toupper(celltype_metadata$species), "_rep", celltype_metadata$replicate),
  paste0("par", toupper(celltype_metadata$species), "_rep", celltype_metadata$replicate)
)

cat("Samples:\n")
for(i in 1:nrow(celltype_metadata)) {
  cat("  ", celltype_metadata$sample_name[i], ":", celltype_metadata$dss_filename[i], "\n")
}

#===============================================================================
# STEP 4: VALIDATE FILES
#===============================================================================

missing_files <- celltype_metadata$dss_filepath[!file.exists(celltype_metadata$dss_filepath)]
if(length(missing_files) > 0) {
  cat("ERROR: Missing DSS files:\n")
  for(file in missing_files) cat("  ", file, "\n")
  stop("Missing files")
}

#===============================================================================
# STEP 5: DEFINE BED FILES
# MODIFY PATHS HERE FOR YOUR BED FILES
#===============================================================================

bed_base_dir <- ""
cat("bed_base_dir: ", metadata_file)

bed_files <- character(0)

# Promoters (modify path as needed)
#promoter_file <- ""
#if(file.exists(promoter_file)) bed_files["promoters"] <- promoter_file
#"promoters" will be the name of regional subfolder

# Add additional BED files here as needed:
# enhancer_file <- file.path(bed_base_dir, paste0(CELLTYPE, "_enhancers.bed"))
# if(file.exists(enhancer_file)) bed_files["enhancers"] <- enhancer_file
#"enhancers" will be the name of regional subfolder

cat("\nBED files:\n")
for(bed_name in names(bed_files)) {
  cat("  ", bed_name, ":", bed_files[bed_name], "\n")
}

if(length(bed_files) == 0) stop("No BED files found")

missing_bed <- bed_files[!file.exists(unlist(bed_files))]
if(length(missing_bed) > 0) stop("Missing BED files: ", paste(missing_bed, collapse = ", "))

#===============================================================================
# STEP 6: RUN ANALYSIS
# ADJUST PARAMETERS HERE
#===============================================================================

output_directory <- paste0("cistrans_results_", CELLTYPE)

results <- RunCompleteMultiBEDPipeline(
  methylation_files = celltype_metadata$dss_filepath,
  sample_names = celltype_metadata$sample_name,
  bed_files = bed_files,
  main_output_dir = output_directory,
  
  # ADJUSTABLE PARAMETERS
  alpha = 0.05,                     # Significance threshold
  min_coverage_per_sample = 5,      # Min coverage per sample (increase for WGBS)
  min_samples_with_coverage = NULL, # NULL = all samples required
  min_cpgs_per_region = 3,          # Min CpGs per region
  min_total_coverage = 15,          # Min total coverage per region
  
  save_all_intermediates = TRUE,
  create_all_plots = TRUE,
  calculate_confidence = TRUE,
  verbose = TRUE
)

#===============================================================================
# STEP 7: SAVE RESULTS AND SUMMARY
#===============================================================================

results_file <- file.path(results$main_output_dir, paste0("complete_results_", CELLTYPE, ".RData"))
save(results, celltype_metadata, file = results_file)
summary_file <- file.path(results$main_output_dir, paste0("SUMMARY_", CELLTYPE, ".txt"))
sink(summary_file)
cat("ANALYSIS SUMMARY FOR", CELLTYPE, "\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("Runtime:", round(results$runtime_minutes, 1), "minutes\n")
cat("Individual CpGs:", nrow(results$individual_results$cpg_results), "\n\n")

cat("CLASSIFICATION RESULTS:\n")
individual_summary <- results$individual_results$summary
for(class in names(individual_summary)) {
  count <- individual_summary[class]
  percent <- round(100 * count / sum(individual_summary), 1)
  cat("  ", class, ":", count, "(", percent, "%)\n")
}
sink()
cat("\nAnalysis complete for", CELLTYPE, "\n")
cat("Results:", results$main_output_dir, "\n")