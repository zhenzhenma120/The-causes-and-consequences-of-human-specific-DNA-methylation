# S1_S2_betabinomial_regulation_classification

A two-model beta-binomial framework for classifying cis and trans regulatory effects on DNA methylation from allele-specific data in hybrid systems.


## Overview

This pipeline analyzes methylation differences between parental and hybrid alleles to classify CpG sites and genomic regions into regulatory categories:
- **pure_cis**: Only hybrid allele effect significant (cis-acting factors), Significant allelic differences in both systems with no significant interaction, indicating consistent cis-acting effects
- **pure_trans**: Only interaction effect significant (trans-acting factors), Significant allelic difference in parents but not hybrids, indicating trans-acting regulation 
- **cis_plus_trans**: Both significant, effects in same direction, Significant allelic differences in both systems with significant interaction, where effects are in the same direction 
- **cis_x_trans**: Both significant, effects in opposite directions, Significant allelic differences in both systems with significant interaction, where effects occur in opposite directions
- **conserved**: Neither effect significant, No significant allelic differences in either generation (neither parental nor hybrid allelic tests significant)

### Input Data Format

**Methylation files (DSS input format, can be created from bismark.cov.gz files, tab-separated; can also load BSSeq object if already created) for individual CpG site classification:**
```
chr    pos    N    X
1      1000   50   25
1      1020   45   20
```
- `chr`: Chromosome (with or without "chr" prefix)
- `pos`: Genomic position
- `N`: Total read coverage
- `X`: Methylated read count

**BED files in addition to Methylation files (standard format) for analysis of pooled methylation data across predefined regions:**
```
chr1    1000    2000    region_1
chr1    3000    4000    region_2
```
- 0-based coordinates
- Name column optional

### Metadata Sample Naming Convention
Example metadata:
- `species`: H for Human and C for Chimp
- `system`: Par for parental and Hyb for hybrid

```
sample.id,filename,cell.type,species,system,replicate
CpG_6401,18146XD-64-01_S0_L001_R1_001_val_1_bismark_bt2_pe_DSS.txt,CNCC,C,Par,1
CpG_6402,18146XD-64-02_S0_L001_R1_001_val_1_bismark_bt2_pe_DSS.txt,CNCC,C,Par,2
CpG_6403,18146XD-64-03_S0_L001_R1_001_val_1_bismark_bt2_pe_DSS.txt,CNCC,H,Par,1
CpG_6404_Halt,18146XD-64-04_S0_L001_R1_001_Complete_val_1_bismark_bt2_pe_sorted_Halt_DSS.txt,CNCC,C,Hyb,1
CpG_6404_Href,18146XD-64-04_S0_L001_R1_001_Complete_val_1_bismark_bt2_pe_sorted_Href_DSS.txt,CNCC,H,Hyb,1
CpG_6405_Halt,18146XD-64-05_S0_L001_R1_001_val_1_bismark_bt2_pe_sorted_Halt_DSS.txt,CNCC,C,Hyb,2
CpG_6405_Href,18146XD-64-05_S0_L001_R1_001_val_1_bismark_bt2_pe_sorted_Href_DSS.txt,CNCC,H,Hyb,2
CpG_6209,18146XD-62-09_S13_L005_R1_001_val_1_bismark_bt2_pe_DSS.txt,SKM,H,Par,1
CpG_6210,18146XD-62-10_S14_L005_R1_001_val_1_bismark_bt2_pe_DSS.txt,SKM,C,Par,1
CpG_6211_Halt,18146XD-62-11_S0_L001_R1_001_val_1_bismark_bt2_pe_sorted_Halt_DSS.txt,SKM,C,Hyb,1
CpG_6211_Href,18146XD-62-11_S0_L001_R1_001_val_1_bismark_bt2_pe_sorted_Href_DSS.txt,SKM,H,Hyb,1
CpG_6212_Halt,18146XD-62-12_S0_L001_R1_001_val_1_bismark_bt2_pe_sorted_Halt_DSS.txt,SKM,C,Hyb,2
CpG_6212_Href,18146XD-62-12_S0_L001_R1_001_val_1_bismark_bt2_pe_sorted_Href_DSS.txt,SKM,H,Hyb,2
```

Sample names are then generated as shown followed by a second column with filename:
- `parH*` - Parental H allele (e.g., parH_rep1, parH_rep2)
- `parC*` - Parental C allele (e.g., parC_rep1)
- `hybH*` - Hybrid H allele (e.g., hybH_rep1)
- `hybC*` - Hybrid C allele (e.g., hybC_rep1, hybC_rep2)

## Pipeline Files

| File | Description |
|------|-------------|
| `00_setup_script.sh` | Validate inputs |
| `01_data_processing.R` | Data loading and filtering |
| `02_regional_analysis.R` | Regional aggregation and analysis |
| `03_visualization.R` | Plotting functions |
| `04_main_multibed_pipeline.R` | Main pipeline for pooled methylation across predefined regions |
| `05_CIS_TRANS_CORE.R` | Core beta-binomial fitting |
| `06_RUN_ANALYSIS.R` | Cell type-specific analysis script |
| `07_submit_all_celltypes.sh` | analysis for the current paper |

## Quick Start


### Configure Analysis

Edit `06_RUN_ANALYSIS.R`:
- Set paths to metadata 
- Set paths to BED files if doing regional analysis
- Adjust input parameters (below)

### Run Analysis

```bash
Rscript 06_RUN_ANALYSIS.R <celltype>
```

```bash
./07_submit_all_celltypes.sh
```

## Input params

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha` | 0.05 | Significance threshold |
| `min_coverage_per_sample` | 5 | Minimum reads per sample per CpG |
| `min_cpgs_per_region` | 3 | Minimum CpGs required per region |
| `min_total_coverage` | 50 | Minimum total coverage per region |

## Output

```
cistrans_results_<CELLTYPE>/
├── 01_individual_cpg_analysis/
│   ├── individual_cpg_results.csv
│   ├── individual_cpg_scatter_plot.png
│   ├── individual_cpg_log2fc_plot.png
│   ├── individual_cpg_summary_plot.png
│   └── individual_analysis_complete.RData
├── 02_regional_<bed_name>/
│   ├── regional_results.csv
│   ├── heterogeneity_analysis.csv
│   └── plots/
├── 03_cross_bed_comparison/
│   └── (comparison plots if multiple BED files)
├── cross_bed_summary.csv
└── report.txt
```
### Individual CpG Results (individual_cpg_results.csv)
| Column | Description |
|--------|-------------|
| chr, pos | Genomic coordinates |
| classification | Regulatory category |
| p_hybrid | P-value for hybrid H vs C effect |
| p_interaction | P-value for allele × generation interaction |
| parent_HC_diff | Methylation difference in parents (H-C) |
| hybrid_HC_diff | Methylation difference in hybrids (H-C) |

### Regional Results (regional_results.csv)
| Column | Description |
|--------|-------------|
| region_id, region_name | Region identifiers |
| chr, start, end | Genomic coordinates |
| meets_requirements | Whether region passed min_coverage_per_sample, min_cpgs_per_region, and min_total_coverage filtering |
| classification | Regulatory category |
| n_cpgs | Number of CpGs in region |
| support_fraction | Fraction of CpGs supporting classification |
| confidence_level | high/medium/low confidence based on support fraction |

## Citations
Hao Wu, Chi Wang, Zhijin Wu (2013):,
               A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data., 
                Biostatistics, 14(2):232-43. doi:10.1093/biostatistics/kxs033

Hao Feng, Karen Conneely, Hao Wu (2014):,
               A bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data., 
                Nucleic acids research, 42(8):e69. doi:10.1093/nar/gku154

Hao Wu, Tianlei Xu, Hao Feng, Li Chen, Ben Li, Bing Yao, Zhaohui Qin, Peng Jin and Karen N. Conneely (2015):, 
               Detection of differentially methylated regions from whole-genome bisulfite sequencing data without replicates.,
                Nucleic acids research. doi: 10.1093/nar/gkv715

Yongseok Park, Hao Wu(2016):,
               Differential methylation analysis for BS-seq data under general experimental design.,
                Bioinformatics. doi:10.1093/bioinformatics/btw026
