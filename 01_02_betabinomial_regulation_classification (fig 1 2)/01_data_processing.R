#===============================================================================
# 01_data_processing.R
#===============================================================================
# DESCRIPTION:
#   Functions for loading methylation data, creating BSseq objects, and 
#   processing BED files for regional analysis.
#
# INPUT FORMATS:
#   Methylation files (DSS format, tab-separated):
#     chr    pos    N    X
#     1      1000   50   25
#     (chr=chromosome, pos=position, N=total reads, X=methylated reads)
#
#   BED files (standard BED format):
#     chr    start    end    [name]
#     1      1000     2000   region_1
#     (0-based coordinates, name column optional)
#
# ADJUSTABLE PARAMETERS:
#   min_coverage_per_sample  - Minimum read coverage per sample (default: 5)
#   min_samples_with_coverage - Minimum samples with coverage (default: all)
#   remove_sex_chromosomes   - Remove X/Y chromosomes (default: TRUE)
#   remove_mitochondrial     - Remove MT chromosome (default: TRUE)
#
# OUTPUTS:
#   BSseq objects containing filtered methylation data
#   GRanges objects containing genomic regions
#
# DEPENDENCIES:
#   bsseq, GenomicRanges, data.table
#===============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(data.table)
})

#-------------------------------------------------------------------------------
# DATA LOADING FUNCTIONS
#-------------------------------------------------------------------------------

LoadMethylationFiles <- function(file_paths, sample_names) {
  # Load methylation data from DSS-formatted files
  if(length(file_paths) != length(sample_names)) {
    stop("Number of files must match number of sample names")
  }
  
  data_list <- list()
  for(i in seq_along(file_paths)) {
    if(!file.exists(file_paths[i])) stop("File not found: ", file_paths[i])
    
    dat <- read.table(file_paths[i], header = TRUE, stringsAsFactors = FALSE)
    required_cols <- c("chr", "pos", "N", "X")
    if(!all(required_cols %in% colnames(dat))) {
      stop("File must have columns: ", paste(required_cols, collapse = ", "))
    }
    
    dat$chr <- gsub("^chr", "", dat$chr)
    if(any(dat$X > dat$N)) stop("Invalid data: X > N in ", file_paths[i])
    
    data_list[[i]] <- dat
  }
  names(data_list) <- sample_names
  return(data_list)
}

CreateBSseqObject <- function(data_list, sample_names = NULL) {
  # Create BSseq object from loaded data
  if(is.null(sample_names)) sample_names <- names(data_list)
  
  # Validate sample names for cis-trans analysis
  has_parH <- any(grepl("parH", sample_names, ignore.case = TRUE))
  has_parC <- any(grepl("parC", sample_names, ignore.case = TRUE))
  has_hybH <- any(grepl("hybH", sample_names, ignore.case = TRUE))
  has_hybC <- any(grepl("hybC", sample_names, ignore.case = TRUE))
  
  if(!(has_parH && has_parC && has_hybH && has_hybC)) {
    warning("Sample names should include parH, parC, hybH, hybC patterns")
  }
  
  BSobj <- makeBSseqData(data_list, sample_names)
  return(BSobj)
}

LoadAndCreateBSseq <- function(file_paths, sample_names, output_file = NULL) {
  # Combined loading and BSseq creation
  data_list <- LoadMethylationFiles(file_paths, sample_names)
  BSobj <- CreateBSseqObject(data_list, sample_names)
  if(!is.null(output_file)) save(BSobj, file = output_file)
  return(BSobj)
}

FilterBSseqObject <- function(BSobj, 
                              min_coverage_per_sample = 5,
                              min_samples_with_coverage = NULL,
                              remove_sex_chromosomes = TRUE,
                              remove_mitochondrial = TRUE) {
  # Filter BSseq object by coverage and chromosome
  if(is.null(min_samples_with_coverage)) {
    min_samples_with_coverage <- ncol(BSobj)
  }
  
  cov_matrix <- getCoverage(BSobj, type = "Cov")
  adequate_coverage <- rowSums(cov_matrix >= min_coverage_per_sample) >= min_samples_with_coverage
  
  gr <- granges(BSobj)
  keep_chr <- rep(TRUE, length(gr))
  
  if(remove_sex_chromosomes) {
    keep_chr <- keep_chr & !seqnames(gr) %in% c("X", "Y", "chrX", "chrY")
  }
  if(remove_mitochondrial) {
    keep_chr <- keep_chr & !seqnames(gr) %in% c("M", "MT", "chrM", "chrMT")
  }
  
  BSobj_filtered <- BSobj[adequate_coverage & keep_chr, ]
  return(BSobj_filtered)
}

#-------------------------------------------------------------------------------
# BED FILE FUNCTIONS
#-------------------------------------------------------------------------------

LoadRegionsFromBED <- function(bed_file, region_name = "custom") {
  # Load genomic regions from BED file
  if(!file.exists(bed_file)) stop("BED file not found: ", bed_file)
  
  # Detect header
  first_line <- readLines(bed_file, n = 1)
  first_cols <- strsplit(first_line, "\t")[[1]]
  has_header <- length(first_cols) >= 3 && 
    (is.na(suppressWarnings(as.numeric(first_cols[2]))) ||
     is.na(suppressWarnings(as.numeric(first_cols[3]))))
  
  bed_data <- read.table(bed_file, header = has_header, stringsAsFactors = FALSE, sep = "\t")
  colnames(bed_data)[1:3] <- c("chr", "start", "end")
  bed_data$chr <- gsub("^chr", "", bed_data$chr)
  
  if(ncol(bed_data) < 4) {
    bed_data$name <- paste0(region_name, "_", 1:nrow(bed_data))
  } else {
    colnames(bed_data)[4] <- "name"
  }
  
  regions_gr <- GRanges(
    seqnames = bed_data$chr,
    ranges = IRanges(start = bed_data$start + 1, end = bed_data$end),
    name = bed_data$name
  )
  return(regions_gr)
}

CreateGenomicTiles <- function(BSobj, tile_width = 100) {
  # Create genomic tiles for tiling analysis
  gr <- granges(BSobj)
  chr_info <- split(gr, seqnames(gr))
  all_tiles <- GRanges()
  
  for(chr_name in names(chr_info)) {
    chr_sites <- chr_info[[chr_name]]
    chr_start <- min(start(chr_sites))
    chr_end <- max(start(chr_sites))
    
    tile_starts <- seq(chr_start, chr_end, by = tile_width)
    tile_ends <- pmin(tile_starts + tile_width - 1, chr_end)
    valid_tiles <- tile_starts <= chr_end
    
    chr_tiles <- GRanges(
      seqnames = chr_name,
      ranges = IRanges(start = tile_starts[valid_tiles], end = tile_ends[valid_tiles]),
      name = paste0(chr_name, "_", tile_starts[valid_tiles], "_", tile_ends[valid_tiles])
    )
    all_tiles <- c(all_tiles, chr_tiles)
  }
  return(all_tiles)
}

#-------------------------------------------------------------------------------
# SAMPLE VALIDATION
#-------------------------------------------------------------------------------

ValidateSampleNamesForCisTrans <- function(sample_names, verbose = TRUE) {
  # Validate sample names contain required patterns
  has_parH <- any(grepl("parH", sample_names, ignore.case = TRUE))
  has_parC <- any(grepl("parC", sample_names, ignore.case = TRUE))
  has_hybH <- any(grepl("hybH", sample_names, ignore.case = TRUE))
  has_hybC <- any(grepl("hybC", sample_names, ignore.case = TRUE))
  
  if(!(has_parH && has_parC && has_hybH && has_hybC)) {
    if(verbose) cat("Sample validation failed. Required: parH, parC, hybH, hybC\n")
    return(FALSE)
  }
  return(TRUE)
}

#-------------------------------------------------------------------------------
# TEST DATA GENERATION
#-------------------------------------------------------------------------------

GenerateQuadrantScatteredTestData <- function(output_dir = "quadrant_test_data",
                                              n_cpgs = 20000,
                                              n_regions = 200,
                                              seed = 456) {
  # Generate test data with CpGs scattered across quadrants
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  set.seed(seed)
  
  chromosomes <- c("1", "2", "3", "4", "5")
  positions_per_chr <- n_cpgs %/% length(chromosomes)
  
  all_positions <- data.frame()
  for(chr in chromosomes) {
    chr_positions <- data.frame(
      chr = paste0("chr", chr),
      pos = sort(sample(5000:500000, positions_per_chr, replace = FALSE))
    )
    all_positions <- rbind(all_positions, chr_positions)
  }
  
  sample_files <- c("parH_rep1.txt", "parH_rep2.txt", "parC_rep1.txt", 
                    "hybH_rep1.txt", "hybC_rep1.txt", "hybC_rep2.txt")
  
  n_total_cpgs <- nrow(all_positions)
  cpg_quadrants <- sample(1:4, n_total_cpgs, replace = TRUE)
  
  cpg_effects <- data.frame(
    cpg_index = 1:n_total_cpgs,
    quadrant = cpg_quadrants,
    parent_effect = numeric(n_total_cpgs),
    hybrid_effect = numeric(n_total_cpgs)
  )
  
  for(i in 1:n_total_cpgs) {
    quadrant <- cpg_quadrants[i]
    if(quadrant == 1) {
      base_effect <- rnorm(1, 0.15, 0.05)
      cpg_effects$parent_effect[i] <- base_effect
      cpg_effects$hybrid_effect[i] <- base_effect + rnorm(1, 0, 0.02)
    } else if(quadrant == 2) {
      cpg_effects$parent_effect[i] <- rnorm(1, -0.1, 0.05)
      cpg_effects$hybrid_effect[i] <- rnorm(1, 0.12, 0.05)
    } else if(quadrant == 3) {
      base_effect <- rnorm(1, -0.15, 0.05)
      cpg_effects$parent_effect[i] <- base_effect
      cpg_effects$hybrid_effect[i] <- base_effect + rnorm(1, 0, 0.02)
    } else {
      cpg_effects$parent_effect[i] <- rnorm(1, 0.15, 0.05)
      cpg_effects$hybrid_effect[i] <- rnorm(1, -0.1, 0.05)
    }
  }
  
  for(i in seq_along(sample_files)) {
    sample_base <- gsub("_rep[12].txt", "", sample_files[i])
    coverage <- pmax(rpois(n_total_cpgs, 50), 10)
    meth_prob <- rep(0.4, n_total_cpgs)
    
    for(j in 1:n_total_cpgs) {
      if(grepl("H", sample_base)) {
        if(grepl("par", sample_base)) {
          meth_prob[j] <- meth_prob[j] + cpg_effects$parent_effect[j]
        } else {
          meth_prob[j] <- meth_prob[j] + cpg_effects$hybrid_effect[j]
        }
      } else {
        if(grepl("par", sample_base)) {
          meth_prob[j] <- meth_prob[j] + rnorm(1, -0.02, 0.03)
        } else {
          meth_prob[j] <- meth_prob[j] + rnorm(1, -0.01, 0.03)
        }
      }
    }
    
    meth_prob <- pmax(0.01, pmin(0.99, meth_prob))
    methylated <- rbinom(n_total_cpgs, coverage, meth_prob)
    
    output_data <- data.frame(
      chr = all_positions$chr,
      pos = all_positions$pos,
      N = coverage,
      X = methylated
    )
    write.table(output_data, file.path(output_dir, sample_files[i]), 
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Create regions
  regions_per_chr <- n_regions %/% length(chromosomes)
  all_regions <- data.frame()
  
  for(chr in chromosomes) {
    chr_positions <- all_positions[all_positions$chr == paste0("chr", chr), ]
    chr_start <- min(chr_positions$pos)
    chr_end <- max(chr_positions$pos)
    region_width <- (chr_end - chr_start) / (regions_per_chr + 1)
    
    region_starts <- seq(chr_start, chr_end - region_width, length.out = regions_per_chr)
    region_ends <- pmin(region_starts + region_width, chr_end)
    
    chr_regions <- data.frame(
      chr = paste0("chr", chr),
      start = round(region_starts),
      end = round(region_ends),
      name = paste0("chr", chr, "_region_", 1:regions_per_chr)
    )
    all_regions <- rbind(all_regions, chr_regions)
  }
  
  bed_file <- file.path(output_dir, "dense_cpg_regions.bed")
  write.table(all_regions[, 1:4], bed_file, sep = "\t", row.names = FALSE, 
              col.names = FALSE, quote = FALSE)
  
  return(list(
    data_files = file.path(output_dir, sample_files),
    bed_file = bed_file,
    sample_names = gsub(".txt", "", sample_files)
  ))
}