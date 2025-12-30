#!/usr/bin/env Rscript
# ==============================================================================
# DMR GENOMIC FEATURE ENRICHMENT AND REGULATION ABUNDANCE
# ==============================================================================
#
# DESCRIPTION:
#   1. GENOMIC ENRICHMENT ANALYSIS
#      - Calculates observed vs expected DMR overlap with genomic features
#      - Features include: gene elements (UTRs, exons, introns, upstream/downstream)
#        and putative cCREs (regulatory elements exerting cis effects on gene expression, e.g. enhancers, CTCF binding sites)
#      - Outputs log2 fold-change enrichment heatmaps
#
#   2. DMR COMPOSITION VISUALIZATION
#      - Creates stacked bar plots showing regulation type proportions
#      - Separates hypermethylated (Hu>Ch) and hypomethylated (Hu<Ch) DMRs
#      - Shows relative contributions of each regulation mechanisms on methylation divergence across cell types
#
# INPUT:
#   DMR files (from changepoint detection)
#   Genomic feature BED files:
#     - gene.bed, exon.bed, cds.bed, intron.bed, utr5.bed, utr3.bed
#     - GRCh38-ELS.bed (enhancers), GRCh38-CTCF.bed
#
# OUTPUT:
#     - unified_features_log2fc_{maxgap}.csv (enrichment matrix)
#     - enrichment_details_{celltype}_{regulation}_{maxgap}.csv (per-combo details)
#     - unified_genomic_heatmap_{maxgap}.pdf/jpg (heatmap visualization)
#     - bar plots of Hu>Ch DMRs and Hu<Ch DMRs
#
# USAGE:
#   Rscript dmr_genomic_analysis.R [maxgap]
#     Rscript dmr_genomic_analysis.R 5000
#     Rscript dmr_genomic_analysis.R 500
#
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(cowplot)
  library(RColorBrewer)
})

CONFIG <- list(
  input_dir = "changepoint/",
  gene_features_dir = "gene_features/",
  ccre_dir = "ccre/",
  celltypes = c("CNCC", "DA", "IPSC", "SKM", "HEP"),
  regulation_classes = c("pure_cis", "pure_trans", "cis_plus_trans", "cis_x_trans"),
  
  celltype_colors = c(
    "CNCC" = "#7C9C3D",
    "DA" = "#C4D600",
    "HEP" = "#FF8C00",
    "IPSC" = "#87CEEB",
    "SKM" = "#DC143C"
  ),
  
  regulation_colors = c(
    "conserved" = "lightgray",
    "pure_cis" = "#FF4500",
    "pure_trans" = "#4169E1",
    "cis_x_trans" = "#228B22",
    "cis_plus_trans" = "#87CEEB"
  ),
  
  gene_bed_files = list(
    "gene.bed" = "Gene_body",
    "exon.bed" = "Exon",
    "cds.bed" = "CDS",
    "intron.bed" = "Intron",
    "utr5.bed" = "5UTR",
    "utr3.bed" = "3UTR"
  ),
  
  ccre_bed_files = list(
    "GRCh38-ELS.bed" = "Enhancers",
    "GRCh38-CTCF.bed" = "CTCF"
  ),
  
  feature_display = c(
    "Upstream" = "Upstream (2kb)",
    "5UTR" = "5' UTR",
    "3UTR" = "3' UTR",
    "Exon" = "Exon",
    "Intron" = "Intron",
    "Downstream" = "Downstream (3kb)"
  ),
  
  regulation_display = c(
    "pure_cis" = "pure cis",
    "pure_trans" = "pure trans",
    "cis_plus_trans" = "cis + trans",
    "cis_x_trans" = "cis × trans"
  )
)

read_bed_to_granges <- function(bed_file) {
  if (!file.exists(bed_file)) {
    warning("BED file not found: ", bed_file)
    return(GRanges())
  }
  
  bed <- read.table(bed_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  if (ncol(bed) < 3) stop("BED file must have at least 3 columns")
  
  colnames(bed)[1:3] <- c("chr", "start", "end")
  
  gr <- GRanges(
    seqnames = bed$chr,
    ranges = IRanges(start = bed$start + 1, end = bed$end)  # Convert to 1-based
  )
  
  if (ncol(bed) >= 4) mcols(gr)$name <- bed[[4]]
  if (ncol(bed) >= 5) mcols(gr)$score <- bed[[5]]
  
  return(gr)
}


filter_standard_chromosomes <- function(gr) {
  standard_chrs <- c(paste0("chr", c(1:22, "X", "Y")), c(1:22, "X", "Y"))
  existing <- intersect(as.character(seqnames(gr)), standard_chrs)
  
  if (length(existing) == 0) {
    warning("No standard chromosomes found")
    return(GRanges())
  }
  
  gr_filtered <- gr[seqnames(gr) %in% existing]
  keepSeqlevels(gr_filtered, existing, pruning.mode = "coarse")
}
load_all_dmrs <- function(input_dir, celltypes, maxgap, format = "rds") {
  all_dmrs <- list()
  
  for (celltype in celltypes) {
    if (format == "rds") {
      file <- file.path(input_dir, sprintf("%s_changepoint_dmrs_simple_%s_.rds", celltype, maxgap))
      if (file.exists(file)) {
        dmrs <- readRDS(file)
        if (length(dmrs) > 0) {
          dmrs$celltype <- celltype
          all_dmrs[[celltype]] <- dmrs
          cat(sprintf("Loaded %d DMRs for %s\n", length(dmrs), celltype))
        }
      }
    } else {
      file <- file.path(input_dir, sprintf("%s_changepoint_dmrs_simple_%s_.csv", celltype, maxgap))
      if (file.exists(file)) {
        dmrs <- read.csv(file, stringsAsFactors = FALSE)
        if (nrow(dmrs) > 0) {
          dmrs$celltype <- celltype
          all_dmrs[[celltype]] <- dmrs
          cat(sprintf("Loaded %d DMRs for %s\n", nrow(dmrs), celltype))
        }
      }
    }
  }
  
  if (length(all_dmrs) == 0) {
    stop("No DMR data found. Check file paths.")
  }
  
  return(all_dmrs)
}

# ==============================================================================
# ENRICHMENT BY GENOMIC FEATURES
# ==============================================================================

load_genomic_features <- function(config) {
  cat("\n=== Loading Genomic Features ===\n")
  
  gene_features <- list()
  for (bed_file in names(config$gene_bed_files)) {
    feature_name <- config$gene_bed_files[[bed_file]]
    bed_path <- file.path(config$gene_features_dir, bed_file)
    
    gr <- read_bed_to_granges(bed_path)
    gr <- filter_standard_chromosomes(gr)
    gene_features[[feature_name]] <- gr
    cat(sprintf("  %s: %d regions\n", feature_name, length(gr)))
  }
  
  if (length(gene_features[["Gene_body"]]) > 0) {
    gene_bodies <- gene_features[["Gene_body"]]
    
    upstream <- promoters(gene_bodies, upstream = 2000, downstream = 0)
    gene_features[["Upstream"]] <- filter_standard_chromosomes(upstream)
    cat(sprintf("  Upstream: %d regions\n", length(gene_features[["Upstream"]])))
    
    downstream <- flank(gene_bodies, width = 3000, start = FALSE)
    gene_features[["Downstream"]] <- filter_standard_chromosomes(downstream)
    cat(sprintf("  Downstream: %d regions\n", length(gene_features[["Downstream"]])))
  }
  
  ccre_features <- list()
  for (bed_file in names(config$ccre_bed_files)) {
    feature_name <- config$ccre_bed_files[[bed_file]]
    bed_path <- file.path(config$ccre_dir, bed_file)
    
    gr <- read_bed_to_granges(bed_path)
    gr <- filter_standard_chromosomes(gr)
    ccre_features[[feature_name]] <- gr
    cat(sprintf("  %s: %d regions\n", feature_name, length(gr)))
  }
  
  list(gene = gene_features, ccre = ccre_features)
}


create_hierarchical_features <- function(gene_features) {
  
  priority <- c("Upstream", "5UTR", "CDS", "3UTR", "Exon", "Intron", "Gene_body", "Downstream")
  available <- intersect(priority, names(gene_features))
  
  hierarchical <- list()
  hierarchical[[available[1]]] <- gene_features[[available[1]]]
  
  if (length(available) > 1) {
    for (i in 2:length(available)) {
      current_name <- available[i]
      current_regions <- gene_features[[current_name]]
      
      if (length(current_regions) > 0) {
        # Combine all previous features
        prev_features <- do.call(c, hierarchical[1:(i-1)])
        
        # Remove overlapping regions
        overlaps <- findOverlaps(current_regions, prev_features)
        if (length(overlaps) > 0) {
          current_exclusive <- current_regions[-unique(queryHits(overlaps))]
        } else {
          current_exclusive <- current_regions
        }
        
        hierarchical[[current_name]] <- current_exclusive
        cat(sprintf("  %s: %d regions (exclusive)\n", current_name, length(current_exclusive)))
      } else {
        hierarchical[[current_name]] <- GRanges()
      }
    }
  }
  
  return(hierarchical)
}

#' Calculate enrichment for DMRs in genomic features
calculate_enrichment <- function(all_dmrs, all_features) {
  cat("\n=== Calculating Enrichment ===\n")
  
  all_combined <- do.call(c, all_features)
  universe <- reduce(all_combined)
  universe_length <- sum(width(universe))
  cat(sprintf("Universe size: %.2f Mb\n", universe_length / 1e6))
  
  reduced_features <- lapply(all_features, reduce)
  feature_lengths <- sapply(reduced_features, function(gr) sum(width(gr)))
  
  log2fc_results <- list()
  enrichment_details <- list()
  
  for (celltype in names(all_dmrs)) {
    dmrs <- all_dmrs[[celltype]]
    regulations <- unique(dmrs$target_class)
    
    for (regulation in regulations) {
      reg_dmrs <- dmrs[dmrs$target_class == regulation]
      if (length(reg_dmrs) == 0) next
      
      combo_name <- paste(celltype, regulation, sep = "_")
      total_dmrs <- length(reg_dmrs)
      
      results <- sapply(names(all_features), function(fname) {
        feature_gr <- reduced_features[[fname]]
        feature_len <- feature_lengths[fname]
        
        # Observed overlaps
        overlaps <- findOverlaps(reg_dmrs, feature_gr)
        observed <- length(unique(queryHits(overlaps)))
        
        # Expected overlaps
        expected <- total_dmrs * (feature_len / universe_length)
        
        # Log2 fold enrichment
        if (observed == 0 && expected > 0) {
          log2(0.1 / expected)  # Pseudocount for zero observed
        } else if (observed > 0 && expected > 0) {
          log2(observed / expected)
        } else {
          NA
        }
      })
      
      log2fc_results[[combo_name]] <- round(results, 3)
      cat(sprintf("  %s: %d DMRs\n", combo_name, total_dmrs))
    }
  }
  log2fc_matrix <- do.call(rbind, log2fc_results)
  
  list(
    matrix = log2fc_matrix,
    universe_length = universe_length,
    feature_lengths = feature_lengths
  )
}

create_enrichment_heatmap <- function(enrichment, config, output_dir, maxgap) {
  cat("\n=== Creating Enrichment Heatmap ===\n")
  
  mat <- enrichment$matrix
  
  row_meta <- data.frame(
    combo = rownames(mat),
    celltype = sapply(strsplit(rownames(mat), "_"), `[`, 1),
    regulation = sapply(strsplit(rownames(mat), "_"), function(x) paste(x[-1], collapse = "_")),
    stringsAsFactors = FALSE
  )
  row_meta <- row_meta[row_meta$regulation != "conserved", ]
  mat <- mat[row_meta$combo, , drop = FALSE]
  
  row_meta <- row_meta[order(row_meta$regulation, row_meta$celltype), ]
  mat <- mat[row_meta$combo, , drop = FALSE]
  
  mat_t <- t(mat)
  
  gene_names <- c("Upstream", "5UTR", "Intron", "Exon", "3UTR", "Downstream")
  ccre_names <- c("Enhancers", "CTCF")
  
  gene_mat <- mat_t[intersect(gene_names, rownames(mat_t)), , drop = FALSE]
  ccre_mat <- mat_t[intersect(ccre_names, rownames(mat_t)), , drop = FALSE]
  
  rownames(gene_mat) <- sapply(rownames(gene_mat), function(x) {
    if (x %in% names(config$feature_display)) config$feature_display[x] else x
  })
  
  colnames(gene_mat) <- gsub("_", " ", colnames(gene_mat))
  colnames(ccre_mat) <- gsub("_", " ", colnames(ccre_mat))
  
  col_anno <- columnAnnotation(
    `Cell Type` = row_meta$celltype,
    Regulation = gsub("_", " ", row_meta$regulation),
    col = list(
      `Cell Type` = config$celltype_colors[row_meta$celltype],
      Regulation = setNames(
        config$regulation_colors[row_meta$regulation],
        gsub("_", " ", row_meta$regulation)
      )
    ),
    annotation_height = unit(c(1, 1), "cm"),
    annotation_name_gp = gpar(fontsize = 12)
  )
  col_fun <- colorRamp2(c(-4, 0, 4), c("#053061", "white", "#67001F"))
  
  ht_gene <- Heatmap(
    gene_mat,
    name = "Log2FC",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 13),
    rect_gp = gpar(col = "gray90", lwd = 0.5),
    top_annotation = col_anno,
    heatmap_legend_param = list(
      title = "Log₂(Obs/Exp)",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      direction = "horizontal"
    ),
    width = unit(12, "cm"),
    height = unit(5, "cm"),
    row_title = "Gene Features",
    row_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  ht_ccre <- Heatmap(
    ccre_mat,
    name = "Log2FC_ccre",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 13),
    rect_gp = gpar(col = "gray90", lwd = 0.5),
    width = unit(12, "cm"),
    height = unit(2, "cm"),
    row_title = "cCREs",
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    show_heatmap_legend = FALSE
  )
  
  ht_list <- ht_gene %v% ht_ccre
  
  pdf(file.path(output_dir, sprintf("enrichment_heatmap_%s.pdf", maxgap)), 
      width = 14, height = 10)
  draw(ht_list, heatmap_legend_side = "right")
  dev.off()
  
  jpeg(file.path(output_dir, sprintf("enrichment_heatmap_%s.jpg", maxgap)),
       width = 12, height = 6, units = "in", res = 600)
  draw(ht_list, heatmap_legend_side = "right")
  dev.off()
  
  cat("Heatmap saved\n")
}

create_composition_barplots <- function(all_dmrs, config, output_dir, maxgap) {
  cat("\n=== Creating Composition Bar Plots ===\n")
  
  dmr_data <- do.call(rbind, lapply(names(all_dmrs), function(ct) {
    df <- as.data.frame(all_dmrs[[ct]])
    df$celltype <- ct
    df
  }))
  
  dmr_data <- dmr_data[dmr_data$target_class != "conserved", ]
  
  dmr_counts <- dmr_data %>%
    group_by(celltype, target_class, dominant_direction) %>%
    summarise(
      dmr_count = n(),
      sites_in_dmrs = sum(n_cpgs),
      .groups = "drop"
    ) %>%
    rename(regulation_class = target_class, direction = dominant_direction)
  
  complete_data <- dmr_counts %>%
    complete(
      celltype = config$celltypes,
      regulation_class = config$regulation_classes,
      direction = c("hyper", "hypo"),
      fill = list(dmr_count = 0, sites_in_dmrs = 0)
    ) %>%
    group_by(celltype, direction) %>%
    mutate(
      total_dmrs = sum(dmr_count),
      dmr_proportion = ifelse(total_dmrs == 0, 0, dmr_count / total_dmrs)
    ) %>%
    ungroup()
  
  complete_data$celltype <- factor(complete_data$celltype, levels = config$celltypes)
  complete_data$regulation_class <- factor(complete_data$regulation_class, 
                                           levels = config$regulation_classes)
  complete_data$direction <- factor(complete_data$direction, levels = c("hyper", "hypo"))
  
  bar_theme <- theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      plot.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  plot_hyper <- ggplot(
    complete_data %>% filter(direction == "hyper"),
    aes(x = celltype, y = dmr_proportion, fill = regulation_class)
  ) +
    geom_bar(stat = "identity", position = "stack", width = 0.7,
             color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = config$regulation_colors,
      name = "Regulation Type",
      labels = config$regulation_display
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), 
                       expand = c(0, 0), limits = c(0, 1)) +
    bar_theme +
    labs(title = "Hu > Ch DMRs", x = "Cell Type", y = "Proportion of DMRs")
  
  plot_hypo <- ggplot(
    complete_data %>% filter(direction == "hypo"),
    aes(x = celltype, y = dmr_proportion, fill = regulation_class)
  ) +
    geom_bar(stat = "identity", position = "stack", width = 0.7,
             color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = config$regulation_colors,
      name = "Regulation Type",
      labels = config$regulation_display
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = c(0, 0), limits = c(0, 1)) +
    bar_theme +
    labs(title = "Hu < Ch DMRs", x = "Cell Type", y = "Proportion of DMRs")
  
  combined <- plot_grid(plot_hyper, plot_hypo, nrow = 1, align = "h", rel_widths = c(1, 1))
  
  ggsave(file.path(output_dir, sprintf("dmr_stacked_bars_%s.pdf", maxgap)),
         combined, width = 14, height = 5)
  ggsave(file.path(output_dir, sprintf("dmr_stacked_bars_%s.jpg", maxgap)),
         combined, width = 14, height = 5, dpi = 600)
  ggsave(file.path(output_dir, sprintf("dmr_hyper_only_%s.png", maxgap)),
         plot_hyper, width = 6, height = 5, dpi = 600)
  ggsave(file.path(output_dir, sprintf("dmr_hypo_only_%s.png", maxgap)),
         plot_hypo, width = 6, height = 5, dpi = 600) 
  return(complete_data)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  maxgap <- if (length(args) >= 1) args[1] else "5000"
    
  output_dir <- file.path(CONFIG$input_dir, sprintf("combined_analysis_%s/", maxgap))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  all_dmrs_gr <- load_all_dmrs(CONFIG$input_dir, CONFIG$celltypes, maxgap, format = "rds")
  all_dmrs_df <- load_all_dmrs(CONFIG$input_dir, CONFIG$celltypes, maxgap, format = "csv")
  
  features <- load_genomic_features(CONFIG)
  hierarchical_gene <- create_hierarchical_features(features$gene)
  # select features for analysis (exclude Gene_body and CDS) and combine gene features with cCRE features
  gene_for_analysis <- hierarchical_gene[c("Upstream", "5UTR", "Intron", "Exon", "3UTR", "Downstream")]
  gene_for_analysis <- gene_for_analysis[sapply(gene_for_analysis, length) > 0]
  all_features <- c(gene_for_analysis, features$ccre)
  enrichment <- calculate_enrichment(all_dmrs_gr, all_features)
  
  write.csv(enrichment$matrix,
            file.path(output_dir, sprintf("enrichment_log2fc_%s.csv", maxgap)))
  create_enrichment_heatmap(enrichment, CONFIG, output_dir, maxgap)
  
  composition_data <- create_composition_barplots(all_dmrs_df, CONFIG, output_dir, maxgap)
  write.csv(composition_data,
            file.path(output_dir, sprintf("dmr_composition_summary_%s.csv", maxgap)),
            row.names = FALSE)
  
  cat(sprintf("\n=== Analysis Complete ===\nResults saved to: %s\n", output_dir))
}
main()
