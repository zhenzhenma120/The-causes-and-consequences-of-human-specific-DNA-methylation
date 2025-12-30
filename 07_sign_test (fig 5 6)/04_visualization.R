library(dplyr)
library(plotrix)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript plot_pathways_fdr.R <output_dir>")
}

output_dir <- args[1]
setwd("SignTest_inputs")

meth_results <- readRDS(paste0(output_dir, "/10_meth_results_fdr_optimized.rds"))
discord_results <- readRDS(paste0(output_dir, "/12_discord_results_fdr_optimized.rds"))
fdr_eligible_pathways <- readRDS(paste0(output_dir, "/08_fdr_eligible_pathways.rds"))


get_fdr_asterisks <- function(fdr_human, fdr_chimp) {
  asterisks <- ""
  direction <- "none"
  
  fdr_to_asterisks <- function(fdr) {
    if (is.na(fdr) || fdr > 0.25) return("")
    if (fdr <= 0.05) return("***")
    if (fdr <= 0.1) return("**")
    return("*")  # fdr <= 0.25
  }
  
  if (!is.na(fdr_human) && fdr_human <= 0.25) {
    asterisks <- fdr_to_asterisks(fdr_human)
    direction <- "human"
  }
  
  if (!is.na(fdr_chimp) && fdr_chimp <= 0.25) {
    chimp_asterisks <- fdr_to_asterisks(fdr_chimp)
    
    if (is.na(fdr_human) || fdr_chimp < fdr_human) {
      asterisks <- chimp_asterisks
      direction <- "chimp"
    }
  }
  
  return(list(asterisks = asterisks, direction = direction))
}

get_label_cex <- function(label, base_cex = 0.8, max_chars = 50) {
  n_chars <- nchar(label)
  if (n_chars <= max_chars) {
    return(base_cex)
  } else {
    shrink_factor <- max_chars / n_chars
    return(max(base_cex * shrink_factor, 0.8))
  }
}

cat("Filtering to FDR-eligible pathways only...\n")
cat(sprintf("Total FDR-eligible pathways: %d\n", length(fdr_eligible_pathways)))

meth_results <- meth_results %>% filter(pathway %in% fdr_eligible_pathways)
discord_results <- discord_results %>% filter(pathway %in% fdr_eligible_pathways)

cat(sprintf("Methylation results after filtering: %d pathways\n", nrow(meth_results)))
cat(sprintf("Discord results after filtering: %d pathways\n", nrow(discord_results)))

# Significant methylation pathways (FDR-eligible, p < 0.05, FDR <= 0.25)
sig_meth <- meth_results %>%
  filter(pmin(p_human, p_chimp) < 0.05) %>%
  filter(pmin(fdr_human, fdr_chimp, na.rm = TRUE) <= 0.25) %>%
  pull(pathway)

# Significant discord pathways (FDR-eligible, p < 0.05, FDR <= 0.25)
sig_discord <- discord_results %>%
  filter(pmin(p_human, p_chimp) < 0.05) %>%
  filter(pmin(fdr_human, fdr_chimp, na.rm = TRUE) <= 0.25) %>%
  pull(pathway)

# Pathways significant in BOTH meth and discord
both_sig <- intersect(sig_meth, sig_discord)


create_plot <- function(plot_data, title, filename) {
  pdf(paste0(output_dir, "/", filename, ".pdf"), 
      width = 5, height = 4.5)
  par(mar = c(4, 4, 3, 1))
  y_pos <- seq(nrow(plot_data), 1)
  
  max_abs <- max(abs(c(plot_data$meth_chimp, plot_data$meth_human,
                       plot_data$discord_chimp, plot_data$discord_human)), na.rm = TRUE)
  x_limit <- max_abs * 1.15
  
  plot(NULL, xlim = c(-x_limit, x_limit),
       ylim = c(0.5, nrow(plot_data) + 0.5), ylab = "", xlab = "Number of genes", 
       yaxt = "n",xaxt = "n", bty = "o", main = "", font = 2, font.lab = 2)
  
  abline(v = 0, col = "black", lwd = 1)
  
  grid_step <- 10
  max_grid <- ceiling(x_limit / grid_step) * grid_step
  grid_seq <- seq(-max_grid , max_grid, by = grid_step)
  grid_seq <- grid_seq[grid_seq != 0]  
  abline(v = grid_seq, col = "gray90", lty = 1)
  axis_step <- 20  
  axis_seq <- seq(-max_grid, max_grid, by = axis_step)
  axis(1, at = axis_seq)
  for(i in 1:nrow(plot_data)) {
    # faint bars (meth)
    rect(plot_data$meth_chimp[i], y_pos[i] - 0.125, 0, y_pos[i] + 0.125, 
         col = "#d1eaff", border = NA)
    rect(0, y_pos[i] - 0.125, plot_data$meth_human[i], y_pos[i] + 0.125, 
         col = "#fae1e1", border = NA)
    
    # bright bars (discord)
    rect(plot_data$discord_chimp[i], y_pos[i] - 0.125, 0, y_pos[i] + 0.125, 
         col = "#61a0ff", border = NA)
    rect(0, y_pos[i] - 0.125, plot_data$discord_human[i], y_pos[i] + 0.125, 
         col = "#f59090", border = NA)
    
    if (plot_data$meth_sig_direction[i] == "chimp" && plot_data$meth_asterisks[i] != "") {
      text(plot_data$meth_chimp[i] - 2.25, y_pos[i], 
           plot_data$meth_asterisks[i], 
           col = "#4a7fa8", cex = 0.8, font = 2, srt = 90, adj = c(0.5, 1))
    } else if (plot_data$meth_sig_direction[i] == "human" && plot_data$meth_asterisks[i] != "") {
      text(plot_data$meth_human[i] + 4, y_pos[i], 
           plot_data$meth_asterisks[i], 
           col = "#c96b6b", cex = 0.8, font = 2, srt = 90, adj = c(0.5, 0))
    }
    
    if (plot_data$discord_sig_direction[i] == "chimp" && plot_data$discord_asterisks[i] != "") {
      text(plot_data$discord_chimp[i] - 3.25, y_pos[i], 
           plot_data$discord_asterisks[i], 
           col = "#2c5fa8", cex = 0.8, font = 2, srt = 90, adj = c(0.5, 1))
    } else if (plot_data$discord_sig_direction[i] == "human" && plot_data$discord_asterisks[i] != "") {
      text(plot_data$discord_human[i] + 5, y_pos[i], 
           plot_data$discord_asterisks[i], 
           col = "#c93838", cex = 0.8, font = 2, srt = 90, adj = c(0.5, 0))
    }
  }
  
  for(i in 1:nrow(plot_data)) {
    label <- plot_data$pathway[i]
    label <- gsub(" - Homo sapiens \\(human\\)", "", label)
    
    label_cex <- get_label_cex(label, base_cex = 0.8)
    
    label_width <- strwidth(label, cex = label_cex)
    label_height <- strheight(label, cex = label_cex)
    
    x_center <- 0
    y_center <- y_pos[i] + 0.45
    padding <- 0.1
    
    symbols(x_center, y_center, rectangles = matrix(c(label_width + padding + 6.5, 
                                                      label_height + padding), 1, 2),
            inches = FALSE, add = TRUE, bg = "white", fg = "gray70", lwd = 1)
    
    text(x_center, y_center, label, cex = label_cex, family = "Helvetica", font = 2)
  }
  
  dev.off()
  
  plot_obj <- list(
    plot_data = plot_data,
    title = title,
    y_pos = y_pos,
    max_abs = max_abs
  )
  saveRDS(plot_obj, paste0(output_dir, "/", filename, ".rds"))
}

# ============================================================================
# Create pathway plots
# Priority 1: Pathways where BOTH meth and discord FDR <= 0.25 and nominal p < 0.05 (ranked by min discord p)
# Priority 2: Remaining significant pathways with discord FDR < 0.25 and both meth and discord nominal p < 0.05 (ranked by min discord p)
# ============================================================================

if (length(both_sig) > 0 || length(sig_meth) > 0 || length(sig_discord) > 0) {
  
  if (length(both_sig) > 0) {
    both_sig_data <- data.frame(pathway = both_sig) %>%
      left_join(meth_results %>% select(pathway, 
                                        p_human_meth = p_human, p_chimp_meth = p_chimp,
                                        fdr_human_meth = fdr_human, fdr_chimp_meth = fdr_chimp), 
                by = "pathway") %>%
      left_join(discord_results %>% select(pathway, 
                                           p_human_discord = p_human, p_chimp_discord = p_chimp,
                                           fdr_human_discord = fdr_human, fdr_chimp_discord = fdr_chimp), 
                by = "pathway") %>%
      mutate(min_p_discord = pmin(p_human_discord, p_chimp_discord, na.rm = TRUE)) %>%
      arrange(min_p_discord) %>%
      mutate(category = "both")
  } else {
    both_sig_data <- data.frame()
  }
  
  remaining_sig <- unique(c(sig_meth, sig_discord))
  remaining_sig <- setdiff(remaining_sig, both_sig)
  
  if (length(remaining_sig) > 0) {
    remaining_data <- data.frame(pathway = remaining_sig) %>%
      left_join(meth_results %>% select(pathway, 
                                        p_human_meth = p_human, p_chimp_meth = p_chimp,
                                        fdr_human_meth = fdr_human, fdr_chimp_meth = fdr_chimp), 
                by = "pathway") %>%
      left_join(discord_results %>% select(pathway, 
                                           p_human_discord = p_human, p_chimp_discord = p_chimp,
                                           fdr_human_discord = fdr_human, fdr_chimp_discord = fdr_chimp), 
                by = "pathway") %>%
      mutate(min_p_discord = pmin(p_human_discord, p_chimp_discord, na.rm = TRUE)) %>%
      arrange(min_p_discord) %>%
      mutate(category = "single")
  } else {
    remaining_data <- data.frame()
  }
  
  if (nrow(both_sig_data) > 0 && nrow(remaining_data) > 0) {
    combined_pathways <- rbind(both_sig_data, remaining_data)
  } else if (nrow(both_sig_data) > 0) {
    combined_pathways <- both_sig_data
  } else {
    combined_pathways <- remaining_data
  }
  

  max_pathways <- 6
  if (nrow(combined_pathways) > max_pathways) {
    n_both <- sum(combined_pathways$category == "both")
    if (n_both >= max_pathways) {
      combined_pathways <- combined_pathways %>% filter(category == "both") %>% head(max_pathways)
    } else {
      n_remaining <- max_pathways - n_both
      both_part <- combined_pathways %>% filter(category == "both")
      single_part <- combined_pathways %>% filter(category == "single") %>% head(n_remaining)
      combined_pathways <- rbind(both_part, single_part)
    }
  }
  
  if (nrow(combined_pathways) > 0) {
    plot_data <- data.frame(
      pathway = combined_pathways$pathway,
      category = combined_pathways$category,
      meth_chimp = -meth_results$n_chimp[match(combined_pathways$pathway, meth_results$pathway)],
      meth_human = meth_results$n_human[match(combined_pathways$pathway, meth_results$pathway)],
      discord_chimp = -discord_results$n_chimp[match(combined_pathways$pathway, discord_results$pathway)],
      discord_human = discord_results$n_human[match(combined_pathways$pathway, discord_results$pathway)]
    )
    
    plot_data$meth_asterisks <- ""
    plot_data$meth_sig_direction <- "none"
    for (i in 1:nrow(plot_data)) {
      pw <- plot_data$pathway[i]
      meth_row <- meth_results[meth_results$pathway == pw, ]
      if (nrow(meth_row) > 0) {
        fdr_info <- get_fdr_asterisks(meth_row$fdr_human, meth_row$fdr_chimp)
        plot_data$meth_asterisks[i] <- fdr_info$asterisks
        plot_data$meth_sig_direction[i] <- fdr_info$direction
      }
    }
    
    plot_data$discord_asterisks <- ""
    plot_data$discord_sig_direction <- "none"
    for (i in 1:nrow(plot_data)) {
      pw <- plot_data$pathway[i]
      discord_row <- discord_results[discord_results$pathway == pw, ]
      if (nrow(discord_row) > 0) {
        fdr_info <- get_fdr_asterisks(discord_row$fdr_human, discord_row$fdr_chimp)
        plot_data$discord_asterisks[i] <- fdr_info$asterisks
        plot_data$discord_sig_direction[i] <- fdr_info$direction
      }
    }
    
    plot_data$meth_chimp[is.na(plot_data$meth_chimp)] <- 0
    plot_data$meth_human[is.na(plot_data$meth_human)] <- 0
    plot_data$discord_chimp[is.na(plot_data$discord_chimp)] <- 0
    plot_data$discord_human[is.na(plot_data$discord_human)] <- 0
    
    n_both <- sum(plot_data$category == "both")
    n_single <- sum(plot_data$category == "single")
    title <- sprintf("Pathway Enrichment (p<0.05 & FDR<=0.25)\nBoth sig: %d | Single sig: %d", 
                     n_both, n_single)
    
    create_plot(plot_data, title, "pathway_enrichment_fdr_combined")
    
    cat(sprintf("\nPlot saved to: %s/pathway_enrichment_fdr_combined.pdf\n", output_dir))
        output_table <- combined_pathways %>%
      select(pathway, category, min_p_discord) %>%
      left_join(meth_results %>% select(pathway, n_genes_tested, n_genes_with_data,
                                        meth_n_human = n_human, meth_n_chimp = n_chimp,
                                        meth_p_human = p_human, meth_p_chimp = p_chimp,
                                        meth_fdr_human = fdr_human, meth_fdr_chimp = fdr_chimp),
                by = "pathway") %>%
      left_join(discord_results %>% select(pathway, 
                                           discord_n_genes_with_data = n_genes_with_data,
                                           discord_n_human = n_human, discord_n_chimp = n_chimp,
                                           discord_p_human = p_human, discord_p_chimp = p_chimp,
                                           discord_fdr_human = fdr_human, discord_fdr_chimp = fdr_chimp),
                by = "pathway")
    
    write.table(output_table, 
                paste0(output_dir, "/pathway_enrichment_fdr_combined.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
  } else {
    cat("No pathways remaining after filtering.\n")
  }
  
} else {
  cat("No significant pathways found with p < 0.05 and FDR <= 0.25.\n")
}
