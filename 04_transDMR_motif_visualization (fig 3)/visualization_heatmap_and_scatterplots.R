.libPaths("~/R/x86_64-pc-linux-gnu-library/4.3")

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(grid)
library(ggplot2)

load_motif_data <- function(base_dir, maxgap,
                            cell_types = c("DA", "CNCC", "IPSC", "HEP", "SKM"),
                            directions = c("hyper", "hypo"),
                            regulation_groups = c("pure_cis", "pure_trans")) {
  
  tables_dir <- file.path(paste0(base_dir, "_motifs"), 
                          paste0("motif_tables_", maxgap, "bp_relaxed/pure_trans_relaxed_out_tables"))
  pvalue_file <- file.path(tables_dir, "table1b_neg_log_p_values.csv")
  
  if (!file.exists(pvalue_file)) stop("P-value file not found: ", pvalue_file)
  
  pvalue_matrix <- as.matrix(read.csv(pvalue_file, row.names = 1, check.names = FALSE))
  mode(pvalue_matrix) <- "numeric"
  
  all_cols <- colnames(pvalue_matrix)
  keep_cols <- c()
  
  for (cell_type in cell_types) {
    for (regulation_group in regulation_groups) {
      for (direction in directions) {
        pattern <- paste0("^", cell_type, '_', regulation_group, "_", direction, "_parent")
        keep_cols <- c(keep_cols, all_cols[grepl(pattern, all_cols)])
      }
    }
  }
  
  keep_cols <- unique(keep_cols)
  if (length(keep_cols) == 0) stop("No matching columns found")
  
  pvalue_matrix[, keep_cols, drop = FALSE]
}

load_expression_data <- function(expression_dir, cell_types) {
  file_patterns <- c(DA = "DA_Parental.txt", CNCC = "CNCC_Parental.txt", 
                     IPSC = "IPSC_Parental.txt", HEP = "HEP_Parental_gilad.txt",
                     SKM = "SKM_Parental_gilad.txt")
  
  expression_data <- list()
  for (cell_type in cell_types) {
    if (cell_type %in% names(file_patterns)) {
      file_path <- file.path(expression_dir, file_patterns[cell_type])
      if (file.exists(file_path)) {
        data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        expression_data[[cell_type]] <- data.frame(
          gene = data$Gene, log2fc = data$log2FoldChange,
          pvalue = data$pvalue, significant = !is.na(data$pvalue) & data$pvalue < 0.05,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  expression_data
}

extract_tf_names <- function(motif_names) {
  tf_names <- gsub("\\(.*", "", motif_names)
  tf_names <- gsub("/.*", "", tf_names)
  trimws(tf_names)
}

create_tf_expression_matrix <- function(motif_matrix, expression_data, selected_cell_types) {
  motif_names <- rownames(motif_matrix)
  tf_names <- extract_tf_names(motif_names)
  available_cell_types <- intersect(selected_cell_types, names(expression_data))
  
  if (length(available_cell_types) == 0) stop("No expression data for selected cell types")
  
  expr_matrix <- matrix(NA, nrow = length(motif_names), ncol = length(available_cell_types))
  sig_matrix <- matrix(FALSE, nrow = length(motif_names), ncol = length(available_cell_types))
  rownames(expr_matrix) <- rownames(sig_matrix) <- motif_names
  colnames(expr_matrix) <- colnames(sig_matrix) <- available_cell_types
  
  for (i in seq_along(motif_names)) {
    tf_name <- tf_names[i]
    for (cell_type in available_cell_types) {
      expr_data <- expression_data[[cell_type]]
      match_idx <- which(toupper(expr_data$gene) == toupper(tf_name))[1]
      if (!is.na(match_idx)) {
        expr_matrix[i, cell_type] <- expr_data$log2fc[match_idx]
        sig_matrix[i, cell_type] <- expr_data$significant[match_idx]
      }
    }
  }
  
  has_expression <- rowSums(!is.na(expr_matrix)) > 0
  list(expression = expr_matrix[has_expression, , drop = FALSE],
       significance = sig_matrix[has_expression, , drop = FALSE],
       kept_indices = which(has_expression))
}

filter_significant_motifs <- function(pvalue_matrix, threshold = 0.05, min_sig = 2) {
  neg_log_threshold <- -log10(threshold)
  sig_counts <- rowSums(pvalue_matrix > neg_log_threshold, na.rm = TRUE)
  keep_motifs <- sig_counts >= min_sig
  sig_counts_filtered <- sig_counts[keep_motifs]
  names(sig_counts)[keep_motifs][order(-sig_counts_filtered, names(sig_counts_filtered))]
}

# ============================================================================
# heatmap fig 3
# ============================================================================

create_cell_fun <- function(sig_matrix = NULL) {
  function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "white", lwd = 0.5, fill = NA))
    if (!is.null(sig_matrix) && sig_matrix[i, j]) {
      grid.rect(x, y, width, height, gp = gpar(col = "black", lwd = 1, fill = NA))
    }
  }
}

add_line_breaks <- function(text_vector, max_chars = 15) {
  sapply(text_vector, function(text) {
    if (nchar(text) <= max_chars) return(text)
    mid_point <- ceiling(nchar(text) / 2)
    spaces <- gregexpr(" ", text)[[1]]
    if (length(spaces) > 1 && spaces[1] != -1) {
      closest_space <- spaces[which.min(abs(spaces - mid_point))]
      paste0(substr(text, 1, closest_space - 1), "\n", substr(text, closest_space + 1, nchar(text)))
    } else {
      paste0(substr(text, 1, mid_point), "\n", substr(text, mid_point + 1, nchar(text)))
    }
  })
}

extract_tf_family_info <- function(names) {
  names_no_slash <- sub("/.*", "", names)
  tf_families <- gsub(".*\\(([^)]+)\\).*", "\\1", names_no_slash)
  tf_families[!grepl("\\(", names_no_slash)] <- "Unknown"
  cleaned_names <- trimws(gsub("\\([^)]*\\)", "", names_no_slash))
  list(families = tf_families, cleaned_names = cleaned_names)
}

create_dual_heatmap <- function(pvalue_matrix, expr_matrix, motif_sig_matrix, expr_sig_matrix) {
  
  motif_plot_matrix <- pvalue_matrix
  motif_plot_matrix[is.na(motif_plot_matrix)] <- 0
  motif_plot_matrix[motif_plot_matrix < 0] <- 0
  max_val <- quantile(motif_plot_matrix[motif_plot_matrix > 0], 0.99, na.rm = TRUE)
  if (is.finite(max_val) && max_val > 10) motif_plot_matrix[motif_plot_matrix > max_val] <- max_val
  
  reorder_by_direction <- function(mat) {
    cols <- colnames(mat)
    hyper <- sort(cols[grepl("hyper", cols, ignore.case = TRUE)])
    hypo <- sort(cols[grepl("hypo", cols, ignore.case = TRUE)])
    other <- sort(cols[!grepl("hyper|hypo", cols, ignore.case = TRUE)])
    mat[, c(hyper, hypo, other), drop = FALSE]
  }
  
  motif_plot_matrix <- reorder_by_direction(motif_plot_matrix)
  motif_sig_matrix <- reorder_by_direction(motif_sig_matrix)
  original_colnames <- colnames(motif_plot_matrix)
  
  family_info <- extract_tf_family_info(rownames(motif_plot_matrix))
  rownames(motif_plot_matrix) <- rownames(motif_sig_matrix) <- family_info$cleaned_names
  rownames(expr_matrix) <- rownames(expr_sig_matrix) <- family_info$cleaned_names
  
  motif_colors <- colorRamp2(seq(0, 10, length.out = 9),
    c("#4D9FD6", "#5FAEDE", "#E3F4FE", "#FFFEE8", "#FAE669", "#F1C232", "#E09A37", "#E76A41", "#9C2007"))
  expr_colors <- colorRamp2(c(-5, 0, 5), c("#4E64B5", "white", "red3"))
  
  category_colors <- c(bZIP = "#B85450", DM = "#376F98", Forkhead = "#CCCCCC", HMG = "#D8F0D8",
    Homeobox = "#90C990", IRF = "#FFFF93", MADS = "black", NR = "#CD9F98", Stat = "#00CED1",
    TEA = "#2D5F2E", TEAD = "#8B7020", Zf = "#2E8B57", Other = "gray60")
  
  simplified_family <- sapply(family_info$families, function(f) {
    if (grepl("bZIP", f, ignore.case = TRUE)) "bZIP"
    else if (grepl("^DM$", f, ignore.case = TRUE)) "DM"
    else if (grepl("Forkhead", f, ignore.case = TRUE)) "Forkhead"
    else if (grepl("HMG", f, ignore.case = TRUE)) "HMG"
    else if (grepl("Homeobox", f, ignore.case = TRUE)) "Homeobox"
    else if (grepl("IRF", f, ignore.case = TRUE)) "IRF"
    else if (grepl("MADS", f, ignore.case = TRUE)) "MADS"
    else if (grepl("^NR$|Nuclear receptor", f, ignore.case = TRUE)) "NR"
    else if (grepl("Stat", f, ignore.case = TRUE)) "Stat"
    else if (grepl("^TEA$", f, ignore.case = TRUE)) "TEA"
    else if (grepl("TEAD", f, ignore.case = TRUE)) "TEAD"
    else if (grepl("^Zf$|zinc finger", f, ignore.case = TRUE)) "Zf"
    else "Other"
  })
  
  simplified_colors <- category_colors[names(category_colors) %in% unique(simplified_family)]
  
  ht_initial <- Heatmap(motif_plot_matrix, name = "init", col = motif_colors,
    cluster_rows = TRUE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
    show_row_names = FALSE, show_column_names = FALSE)
  ht_drawn <- draw(ht_initial, show_heatmap_legend = FALSE)
  row_order <- row_order(ht_drawn)
  row_dend <- row_dend(ht_drawn)
  grid.newpage()
  
  motif_plot_ordered <- motif_plot_matrix[row_order, , drop = FALSE]
  motif_sig_ordered <- motif_sig_matrix[row_order, , drop = FALSE]
  colnames(motif_plot_ordered) <- colnames(motif_sig_ordered) <- 
    sapply(colnames(motif_plot_ordered), function(n) strsplit(n, "_")[[1]][1])
  
  desired_order <- c("CNCC", "DA", "HEP", "IPSC", "SKM")
  available_cols <- intersect(desired_order, colnames(expr_matrix))
  expr_ordered <- expr_matrix[row_order, available_cols, drop = FALSE]
  expr_sig_ordered <- expr_sig_matrix[row_order, available_cols, drop = FALSE]
  simplified_family_ordered <- simplified_family[row_order]
  
  col_annotation <- HeatmapAnnotation(
    `Methylation Bias` = ifelse(grepl("hyper", original_colnames), "Hu>Ch",
      ifelse(grepl("hypo", original_colnames), "Ch>Hu", "Other")),
    col = list(`Methylation Bias` = c("Hu>Ch" = "#B55E5E", "Ch>Hu" = "#3E5270", "Other" = "grey80")),
    annotation_height = unit(3, "mm"), show_annotation_name = FALSE, show_legend = FALSE)
  
  ht_motif <- Heatmap(motif_plot_ordered, name = "motif", col = motif_colors,
    cluster_rows = row_dend, cluster_columns = FALSE,
    column_split = ifelse(grepl("hyper", original_colnames), "Hu>Ch",
      ifelse(grepl("hypo", original_colnames), "Ch>Hu", "Other")),
    bottom_annotation = col_annotation,
    column_title = "Motif Enrichment in trans-DMRs",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_gap = unit(2, "mm"), show_row_names = FALSE, show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8), column_names_rot = 90,
    cell_fun = create_cell_fun(motif_sig_ordered),
    heatmap_legend_param = list(title = "-log10(p-value)", title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8), legend_direction = "horizontal"),
    width = unit(5.5, "cm"),
    left_annotation = rowAnnotation(Family = simplified_family_ordered,
      col = list(Family = simplified_colors), width = unit(5, "mm"),
      annotation_legend_param = list(Family = list(title = "TF Family", ncol = 2))))
  
  row_names_anno <- rowAnnotation(
    names = anno_text(add_line_breaks(family_info$cleaned_names[row_order], 18),
      gp = gpar(fontsize = 7), just = "left", width = unit(2, "cm")),
    width = unit(2, "cm"), show_annotation_name = FALSE)
  
  ht_expr <- Heatmap(expr_ordered, name = "expression", col = expr_colors,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    show_column_names = TRUE, na_col = "grey80", column_names_gp = gpar(fontsize = 8),
    cell_fun = create_cell_fun(NULL),
    heatmap_legend_param = list(title = "log2FC(H/C)", legend_direction = "horizontal"),
    width = unit(2.5, "cm"), column_title = "TF Expression",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"))
  
  list(heatmap = ht_motif + row_names_anno + ht_expr, row_order = row_order,
       cleaned_names = family_info$cleaned_names)
}

# ============================================================================
# tf scatterplot fig 3
# ============================================================================

plot_tf_gene <- function(gene_name, expression_dir, motif_file, output_dir) {
  
  cell_colors <- c(CNCC = "#7C9C3D", DA = "#C4D600", HEP = "#FF8C00", 
                   IPSC = "#87CEEB", SKM = "#DC143C")
  
  parental_files <- list.files(expression_dir, pattern = "*Parental*", full.names = TRUE)
  expression <- list()
  for (file in parental_files) {
    cell_type <- strsplit(basename(file), "_")[[1]][1]
    df <- read.table(file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    gene_match <- grep(paste0("^", gene_name, "$"), rownames(df), ignore.case = TRUE)
    if (length(gene_match) > 0) expression[[cell_type]] <- df[gene_match[1], "log2FoldChange"]
  }
  
  motif_df <- read.csv(motif_file, row.names = 1, stringsAsFactors = FALSE)
  clean_names <- sapply(rownames(motif_df), function(x) strsplit(x, "[(/]")[[1]][1])
  gene_match <- grep(paste0("^", gene_name, "$"), clean_names, ignore.case = TRUE)
  if (length(gene_match) == 0) return(NULL)
  gene_row <- motif_df[gene_match[1], ]
  
  plot_data <- data.frame()
  for (cell in c("CNCC", "DA", "HEP", "IPSC", "SKM")) {
    if (cell %in% names(expression)) {
      for (cond in c("hypo", "hyper")) {
        col_name <- paste0(cell, "_pure_trans_", cond, "_parent")
        if (col_name %in% colnames(gene_row)) {
          plot_data <- rbind(plot_data, data.frame(
            cell_type = cell,
            condition = ifelse(cond == "hyper", "Hu>Ch trans-DMRs", "Ch>Hu trans-DMRs"),
            expression = expression[[cell]],
            motif_fc = gene_row[[col_name]]))
        }
      }
    }
  }
  
  if (nrow(plot_data) == 0) return(NULL)
  
  p <- ggplot(plot_data, aes(x = expression, y = motif_fc)) +
    geom_point(aes(color = cell_type), size = 1.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~condition, ncol = 2, scales = "free") +
    scale_color_manual(values = cell_colors) +
    labs(x = "Expression (log2FC Hu/Ch)",
         y = "Fold-Enrichment of motif\nin Hu>Ch or Ch>Hu trans-DMRs",
         title = gene_name, color = "Cell Type") +
    theme_bw(base_size = 12) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(face = "bold"),
          legend.position = "right")
  
  ggsave(file.path(output_dir, paste0(gene_name, "_correlation.jpg")), 
         p, width = 5, height = 2.5, dpi = 600)
  p
}


run_analysis <- function(base_dir, maxgap, expression_dir,
                         cell_types = c("DA", "CNCC", "IPSC", "HEP", "SKM"),
                         directions = c("hyper", "hypo"),
                         regulation_groups = c("pure_trans"),
                         significance_threshold = 0.01,
                         min_significant = 2,
                         width = 9.5, height = 10) {
  
  pvalue_matrix <- load_motif_data(base_dir, maxgap, cell_types, directions, regulation_groups)
  significant_motifs <- filter_significant_motifs(pvalue_matrix, significance_threshold, min_significant)
  if (length(significant_motifs) == 0) stop("No significant motifs found")
  
  pvalue_filtered <- pvalue_matrix[significant_motifs, , drop = FALSE]
  neg_log_threshold <- -log10(0.001)
  motif_sig_matrix <- !is.na(pvalue_filtered) & pvalue_filtered > neg_log_threshold
  
  expression_data <- load_expression_data(expression_dir, cell_types)
  if (length(expression_data) == 0) stop("No expression data loaded")
  
  tf_data <- create_tf_expression_matrix(pvalue_filtered, expression_data, cell_types)
  pvalue_filtered <- pvalue_filtered[tf_data$kept_indices, , drop = FALSE]
  motif_sig_matrix <- motif_sig_matrix[tf_data$kept_indices, , drop = FALSE]
  
  heatmap_result <- create_dual_heatmap(pvalue_filtered, tf_data$expression, 
                                        motif_sig_matrix, tf_data$significance)
  
  output_dir <- file.path(paste0(base_dir, "_motifs"), paste0("complex_dual_heatmaps_", maxgap, "bp"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  file_base <- paste0("dual_heatmap_", maxgap, "bp_", paste(regulation_groups, collapse = "_"))
  output_pdf <- file.path(output_dir, paste0(file_base, ".pdf"))
  output_jpg <- file.path(output_dir, paste0(file_base, ".jpg"))
  
  pdf(output_pdf, width = width, height = height)
  draw(heatmap_result$heatmap, heatmap_legend_side = "bottom", merge_legend = FALSE)
  decorate_heatmap_body("expression", {grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.5))}, slice = 1)
  dev.off()
  
  jpeg(output_jpg, width = width, height = height, units = "in", res = 600)
  draw(heatmap_result$heatmap, heatmap_legend_side = "bottom", merge_legend = FALSE)
  decorate_heatmap_body("expression", {grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.5))}, slice = 1)
  dev.off()
  
  motif_file <- file.path(paste0(base_dir, "_motifs"),
    paste0("motif_tables_", maxgap, "bp_relaxed/pure_trans_relaxed_out_tables/table2_enrichment_fold_changes.csv"))
  
  for (gene in heatmap_result$cleaned_names) {
    tryCatch(plot_tf_gene(gene, expression_dir, motif_file, output_dir),
             error = function(e) message("Skipped ", gene, ": ", e$message))
  }
  
  list(cleaned_names = heatmap_result$cleaned_names, output_dir = output_dir,
       heatmap_pdf = output_pdf, heatmap_jpg = output_jpg)
}
BASE_DIR <- "Changepoint/"
EXPRESSION_DIR <- "TPM_expression/"

results <- run_analysis(
  base_dir = BASE_DIR,
  maxgap = "5000",
  expression_dir = EXPRESSION_DIR,
  cell_types = c("DA", "CNCC", "IPSC", "SKM", "HEP"),
  directions = c("hyper", "hypo"),
  regulation_groups = c("pure_trans"),
  significance_threshold = 0.01,
  min_significant = 2
)