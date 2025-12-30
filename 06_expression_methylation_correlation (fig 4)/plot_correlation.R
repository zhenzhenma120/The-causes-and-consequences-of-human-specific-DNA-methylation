library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
file1 <- args[1]
file2 <- args[2]
suffix1 <- args[3]
suffix2 <- args[4]
cell_type <- args[5]
output_prefix <- args[6]

read_correlation_file <- function(filepath) {
  data <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(data)
}

extract_correlations <- function(data, species_system, cell_type_filter) {
  subset_data <- data[data$Species_System == species_system & data$Cell_Type == cell_type_filter, ]
  subset_data <- subset_data[order(subset_data$Bin), ]
  return(list(
    correlations = subset_data$Spearman_Correlation,
    bins = subset_data$Bin,
    regions = subset_data$Region_Type
  ))
}

plot_conditions_bar <- function(cell_type, conditions, correlations_list, 
                                title_suffix, bar_colors) {
  vec_length <- length(correlations_list[[1]])
  
  tss_bin <- round(vec_length / 2)
  
  cat(sprintf("Processing %s %s: total bins = %d, TSS at bin %d\n", 
              cell_type, title_suffix, vec_length, tss_bin))
  
  data_matrix <- matrix(NA, nrow = length(conditions), ncol = vec_length)
  
  for (i in 1:length(conditions)) {
    if (!is.null(correlations_list[[i]])) {
      data_matrix[i, ] <- correlations_list[[i]]
    }
  }
  
  y_min <- min(data_matrix, na.rm = TRUE) - 0.1
  y_max <- max(data_matrix, na.rm = TRUE) + 0.1
  #y_min = -0.4
  #y_max = 0.4
  y_range <- y_max - y_min
  
  ylim <- c(y_min, y_max)
  
  bp <- barplot(data_matrix, beside = TRUE, col = bar_colors, 
                ylim = ylim, main = paste(cell_type, title_suffix),
                xlab = "", ylab = "Spearman Correlation", 
                space = c(0, 0.7), las = 1, border = "black", xaxt = "n", yaxt = "n")
  y_ticks <- seq(floor(y_min * 5) / 5, ceiling(y_max * 5) / 5, by = 0.2)
  axis(2, at = y_ticks, las = 1, cex.axis = 1)
  abline(h = 0, col = "grey50", lty = "dashed", lwd = 1.5)
  
  tss_pos <- mean(bp[, tss_bin])
  
  segments(x0 = tss_pos+1.2, y0 = y_min, 
           x1 = tss_pos+1.2, y1 = y_max, 
           col = "gray90", lty = "dashed", lwd = 2)
  
  upstream_mid <- tss_bin / 2
  downstream_mid <- tss_bin + (vec_length - tss_bin) / 2
  
  # text_y_bottom <- y_min + 0.05 * y_range
  # text(mean(bp[, upstream_mid])+2.3, text_y_bottom, "Upstream", cex = 1.0, font = 2, col = "blue")
  # text(mean(bp[, downstream_mid])+2.3, text_y_bottom, "Downstream", cex = 1.0, font = 2, col = "purple")
  
  start_pos <- mean(bp[, 1])
  tss_label_pos <- tss_pos+1
  end_pos <- mean(bp[, ncol(bp)])
  
  axis(1, at = colMeans(bp), labels = FALSE, tcl = -0.3)
  
  axis(1, at = c(start_pos, tss_label_pos+0.4, end_pos), 
       labels = c("-1kb", "TSS", "+1kb"), 
       tick = FALSE, cex.axis = 1.3, padj = -0.2)
}

cat("Reading input files...\n")
data1 <- read_correlation_file(file1)
data2 <- read_correlation_file(file2)

species_systems <- c("Parent_Human", "Parent_Chimp", "Hybrid_Human", "Hybrid_Chimp")

correlations_file1 <- list()
correlations_file2 <- list()

for (species in species_systems) {
  result1 <- extract_correlations(data1, species, cell_type)
  result2 <- extract_correlations(data2, species, cell_type)
  
  correlations_file1[[species]] <- result1$correlations
  correlations_file2[[species]] <- result2$correlations
}

#fig 4G
output_file <- sprintf("%s_%s_tss_comparison.jpg", output_prefix, cell_type)
cat(sprintf("Generating plot: %s\n", output_file))

jpeg(file = output_file, width = 8, height = 8, units = "in", res = 600)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(2, 2, 3, 1), 
    cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.0)

plot_order <- c("Parent_Human", "Parent_Chimp", "Hybrid_Human", "Hybrid_Chimp")
color_pairs <- list(
  Parent_Human = c("coral4", "pink"),
  Parent_Chimp = c("steelblue4", "lightblue"),
  Hybrid_Human = c("coral4", "pink"),
  Hybrid_Chimp = c("steelblue4", "lightblue")
)

for (species in plot_order) {
  conditions <- c(suffix1, suffix2)
  
  correlations_list <- list(
    correlations_file1[[species]],
    correlations_file2[[species]]
  )
  
  plot_conditions_bar(
    cell_type, 
    conditions, 
    correlations_list,
    gsub("_", " ", species), 
    color_pairs[[species]]
  )
}

mtext(sprintf("%s: %s vs %s", cell_type, suffix1, suffix2), 
      outer = TRUE, cex = 1.8, font = 2, line = 0.5)

dev.off()