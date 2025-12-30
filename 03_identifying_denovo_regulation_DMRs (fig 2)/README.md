# Identification of denovo DMR

Differentially Methylated Region (DMR) detection using changepoint analysis and downstream genomic enrichment analysis.

## Overview

1. **Identifies DMRs denovo using results from 01_02_betabinomial_regulation_classification** - Genomic regions enriched for specific methylation regulatory categories (pure_cis, pure_trans, cis_plus_trans, cis_x_trans, conserved) using PELT changepoint detection
2. **Enrichment in genomic features** - Calculates and visualize relative enrichment of DMRs in genomic features (UTRs, exons, enhancers, etc.)
3. **Visualization** - Figure 2 heatmaps and barplots 

## File Structure

```
├── 01_dmr_fixed_penalty.R       # Fixed penalty method
├── 02_dmr_binary_crops.R        # Binary scoring + CROPS optimization
├── 03_dmr_gradient_crops.R      # Gradient scoring + CROPS optimization
└── 04_visualization.R    # Genomic enrichment heatmap + barplot
```

## Input Format

CSV file with columns:
- `chr`: Chromosome number (numeric, e.g., 1, 2, ..., 22)
- `pos`: Genomic position (integer)
- `classification`: One of "pure_cis", "pure_trans", "cis_plus_trans", "cis_x_trans", "conserved"
- `parent_HC_diff`: Methylation difference between parents (-1 to 1)

```
{basedir}/cistrans_results_{celltype}_5bp_coverage/01_individual_cpg_analysis/{celltype}_individual_cpg_results.csv
```

## Output Format

### Per-celltype files

**CSV**: `{celltype}_dmrs_{method}_{maxgap}.csv`
| Column | Description |
|--------|-------------|
| chr | Chromosome (e.g., chr1) |
| start | DMR start position |
| end | DMR end position |
| n_cpgs | Number of CpGs in DMR |
| target_class | Regulatory category |
| target_proportion | Proportion of target class sites |
| mean_effect | Mean absolute methylation difference |
| median_effect | Median methylation difference |
| effect_sd | Standard deviation of effects |
| dominant_direction | "hyper", "hypo", or "stable" |
| directional_consistency | Proportion of effects in dominant direction |
| optimal_penalty | Penalty value used|

### Configurable Parameters

```r
CONFIG <- list(
  min_sites = 5,                      # Min CpGs per DMR
  target_proportion_threshold = 0.5,  # Min proportion of target class
  min_effect_threshold = 0.05,        # Min effect size for non-conserved
  conserved_effect_range = 0.05,      # Max effect for conserved sites
  
  fixed_penalty = 0.75,               # For fixed method
  pen_value_range = c(0.01, 2.0),     # CROPS search range
  
  gradient_proportion_threshold = 0.2,  # Relaxed for gradient
  gradient_edge_threshold = 0.75        # Edge trimming threshold
)
```
## Output Directory Structure

```
09_dmr_fixed/
├── DA_dmrs_fixed_500.csv
├── DA_dmrs_fixed_500.rds
├── CNCC_dmrs_fixed_500.csv
├── ...
└── dmr_list_fixed_500.rds

09_dmr_binary_crops/
├── DA_dmrs_binary_500.csv
├── ...

09_dmr_gradient_crops/
├── DA_dmrs_gradient_500.csv
├── ...
```

### Visualization

**DMR inputs**:
```
{input_dir}/{celltype}_changepoint_dmrs_simple_{maxgap}_.rds
{input_dir}/{celltype}_changepoint_dmrs_simple_{maxgap}_.csv
```

**BED inputs** (genomic features):
```
{gene_features_dir}/gene.bed, exon.bed, cds.bed, intron.bed, utr5.bed, utr3.bed
{ccre_dir}/GRCh38-ELS.bed (enhancers)
{ccre_dir}/GRCh38-CTCF.bed
```

The enrichment is calculated as:

```
Expected = Total_DMRs × (Feature_Length / Length of union of all genomic features being tested)
Observed = Number of DMRs overlapping feature
Log2FC = log2(Observed / Expected)
```