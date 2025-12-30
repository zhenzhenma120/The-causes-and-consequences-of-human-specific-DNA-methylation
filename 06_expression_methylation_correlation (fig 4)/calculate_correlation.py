#!/usr/bin/env python3

"""
This script is adapted from CpGtools from https://cpgtools.readthedocs.io/en/latest/demo/CpG_density_gene_centered.html
Citation:
Wei T, Nie J, Larson NB, Ye Z, Eckel Passow JE, Robertson KD, Kocher JA, Wang L. CpGtools: A Python Package for DNA Methylation Analysis. Bioinformatics. 2019 Dec 6 Epub 2019 Dec 06

Modified to calculate Spearman correlation between methylation levels
and gene expression in regions around the TSS (transcription start site).

"""

import sys
import os
import collections
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED

__author__ = "Modified from Liguo Wang"

def read_expression_data(expression_file):
    return pd.read_csv(expression_file, index_col=0)

def parse_coordinate(coord_str):
    try:
        return int(float(coord_str))
    except ValueError:
        return int(coord_str)

def read_region_bed_with_names(bed_file):
    regions = []
    for line in ireader.reader(bed_file):
        if line.startswith('#') or line.strip() == '':
            continue
        fields = line.strip().split('\t')
        if len(fields) >= 4:
            chrom = fields[0]
            start = parse_coordinate(fields[1])
            end = parse_coordinate(fields[2])
            name = fields[3]
            strand = fields[5] if len(fields) >= 6 else '+'
            regions.append((chrom, start, end, strand, name))
    return regions

def read_CpG_bed_fixed(bed_file):

    import tempfile
    import os
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed')
    
    try:
        with open(bed_file, 'r') as infile:
            for line in infile:
                if line.startswith('#') or line.strip() == '':
                    temp_file.write(line)
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    fields[1] = str(parse_coordinate(fields[1]))  # start
                    fields[2] = str(parse_coordinate(fields[2]))  # end
                
                temp_file.write('\t'.join(fields) + '\n')
        
        temp_file.close()
        
        cpg_ranges = read_CpG_bed(temp_file.name)
        os.unlink(temp_file.name)
        
        return cpg_ranges
        
    except Exception as e:
        temp_file.close()
        if os.path.exists(temp_file.name):
            os.unlink(temp_file.name)
        raise e

def get_methylation_in_region(chrom, start, end, cpg_ranges):
    """Extract methylation values from a genomic region using IntervalTree"""
    methylation_values = []
    
    if chrom in cpg_ranges:
        # Query the interval tree for overlapping CpG sites
        overlapping_intervals = cpg_ranges[chrom].find(start, end)
        
        for interval in overlapping_intervals:
            pos = interval.start  # CpG position
            beta_value = interval.value  # methylation beta value
            
            if start <= pos <= end:
                methylation_values.append((pos, beta_value))
    
    return methylation_values

def bin_methylation_data_tss_fixed(region_list, cpg_ranges, upstream_bins, downstream_bins, upstream_size, downstream_size):
    """
    Bin methylation data around TSS with fixed binning structure:
    - User-specified upstream bins (e.g., 20 bins for 1kb = 50bp per bin)
    - User-specified downstream bins (e.g., 20 bins for 1kb = 50bp per bin)
    
    This ensures TSS is always at the boundary between upstream and downstream
    """
    total_bins = upstream_bins + downstream_bins
    
    print(f"Fixed TSS binning structure:")
    print(f"  Upstream bins: 0-{upstream_bins-1} ({upstream_bins} bins, {upstream_size}bp)")
    print(f"  TSS at bin boundary: {upstream_bins}")
    print(f"  Downstream bins: {upstream_bins}-{total_bins-1} ({downstream_bins} bins, {downstream_size}bp)")
    print(f"  Total bins: {total_bins}")
    
    binned_data = {i: [] for i in range(total_bins)}
    
    for chrom, tx_start, tx_end, strand, gene_name in region_list:
        # Define TSS based on strand
        if strand == '+':
            tss = tx_start
            upstream_region_start = max(tss - upstream_size, 0)
            upstream_region_end = tss
            downstream_region_start = tss
            downstream_region_end = tss + downstream_size
        else:  # negative strand - TSS is at tx_end
            tss = tx_end
            upstream_region_start = tss
            upstream_region_end = tss + upstream_size
            downstream_region_start = max(tss - downstream_size, 0)
            downstream_region_end = tss
        
        # Process upstream region
        upstream_meth = get_methylation_in_region(chrom, upstream_region_start, upstream_region_end, cpg_ranges)
        upstream_length = upstream_region_end - upstream_region_start
        upstream_bin_size = upstream_length / upstream_bins if upstream_bins > 0 else 1
        
        for pos, beta_value in upstream_meth:
            if strand == '+':
                relative_pos = pos - upstream_region_start
            else:  # For negative strand, reverse direction
                relative_pos = upstream_region_end - pos
            
            local_bin_idx = min(int(relative_pos / upstream_bin_size), upstream_bins - 1)
            if 0 <= local_bin_idx < upstream_bins:
                binned_data[local_bin_idx].append((gene_name, beta_value))
        
        # Process downstream region
        downstream_meth = get_methylation_in_region(chrom, downstream_region_start, downstream_region_end, cpg_ranges)
        downstream_length = downstream_region_end - downstream_region_start
        downstream_bin_size = downstream_length / downstream_bins if downstream_bins > 0 else 1
        
        for pos, beta_value in downstream_meth:
            if strand == '+':
                relative_pos = pos - downstream_region_start
            else:  # For negative strand, reverse direction
                relative_pos = downstream_region_end - pos
            
            local_bin_idx = min(int(relative_pos / downstream_bin_size), downstream_bins - 1)
            global_bin_idx = upstream_bins + local_bin_idx
            if upstream_bins <= global_bin_idx < total_bins:
                binned_data[global_bin_idx].append((gene_name, beta_value))
    
    return binned_data, upstream_bins, downstream_bins, total_bins

def get_significance_marker(p_value):
    """
    Convert p-value to significance marker
    * for p < 0.05
    ** for p < 0.01
    *** for p < 0.001
    """
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return ""

def calculate_correlation_profile(binned_methylation, expression_df, expression_column):
    """Calculate Spearman correlation for each bin"""
    correlations = {}
    p_values = {}
    bin_counts = {}
    
    for bin_idx, methylation_data in binned_methylation.items():
        bin_counts[bin_idx] = len(methylation_data)
        
        if len(methylation_data) < 3:  # Need at least 3 points for correlation
            correlations[bin_idx] = 0.0
            p_values[bin_idx] = 1.0
            continue
        
        # Extract gene names and methylation values
        gene_names = [item[0] for item in methylation_data]
        meth_values = [item[1] for item in methylation_data]
        
        # Get corresponding expression values
        expr_values = []
        valid_meth = []
        
        for gene_name, meth_val in zip(gene_names, meth_values):
            if gene_name in expression_df.index:
                expr_val = expression_df.loc[gene_name, expression_column]
                if not pd.isna(expr_val):
                    expr_values.append(expr_val)
                    valid_meth.append(meth_val)
        
        if len(valid_meth) >= 3:
            try:
                corr, p_value = spearmanr(valid_meth, expr_values)
                correlations[bin_idx] = corr if not np.isnan(corr) else 0.0
                p_values[bin_idx] = p_value if not np.isnan(p_value) else 1.0
            except:
                correlations[bin_idx] = 0.0
                p_values[bin_idx] = 1.0
        else:
            correlations[bin_idx] = 0.0
            p_values[bin_idx] = 1.0
    
    return correlations, p_values, bin_counts

def add_tss_markers_to_plot(ROUT, upstream_bins, downstream_bins, total_bins):
    """Add TSS marker and region labels to the correlation plot"""
    tss_bin = upstream_bins  # TSS is at the boundary
    
    # Add vertical line for TSS
    print(f'abline(v={tss_bin}, col="darkgreen", lty="solid", lwd=2.5)', file=ROUT)
    
    # Add TSS label at top
    print(f'text(x={tss_bin}, y=0.95, labels="TSS", col="darkgreen", cex=1.0, pos=3, font=2)', file=ROUT)
    
    # Add region labels at bottom
    upstream_mid = upstream_bins / 2
    downstream_mid = upstream_bins + (downstream_bins / 2)
    
    print(f'text(x={upstream_mid}, y=-0.95, labels="Upstream", cex=0.9, col="blue", font=2)', file=ROUT)
    print(f'text(x={downstream_mid}, y=-0.95, labels="Downstream", cex=0.9, col="purple", font=2)', file=ROUT)

def main():
    usage = "%prog [options]" + "\n"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-m", "--methylation_files", action="store", type="string", dest="methylation_files",
                     help="Text file listing methylation BED6+ files (one per line)")
    parser.add_option("-r", "--region", action="store", type="string", dest="region_file",
                     help="BED3+ file of genomic regions with gene names in 4th column")
    parser.add_option("-e", "--expression", action="store", type="string", dest="expression_file",
                     help="CSV file with gene expression data")
    parser.add_option("-d", "--design", action="store", type="string", dest="design_file",
                     help="CSV design matrix mapping files to colors and line styles")
    parser.add_option("--upstream_bins", action="store", type="int", dest="upstream_bins", default=20,
                     help="Number of bins for upstream region. default=%default")
    parser.add_option("--downstream_bins", action="store", type="int", dest="downstream_bins", default=20,
                     help="Number of bins for downstream region. default=%default")
    parser.add_option("--downstream", action="store", type="int", dest="downstream_size", default=1000,
                     help="Size of downstream region from TSS. default=%default (bp)")
    parser.add_option("--upstream", action="store", type="int", dest="upstream_size", default=1000,
                     help="Size of upstream region from TSS. default=%default (bp)")
    parser.add_option("-o", "--output", action="store", type='string', dest="out_file",
                     help="The prefix of the output file.")
    
    (options, args) = parser.parse_args()
    
    if not all([options.methylation_files, options.region_file, options.expression_file, 
                options.design_file, options.out_file]):
        print(__doc__)
        parser.print_help()
        sys.exit(1)
    
    design_df = pd.read_csv(options.design_file)
    expression_df = read_expression_data(options.expression_file)
    print(f"Available expression columns: {list(expression_df.columns)}")
    printlog("Reading BED file: \"%s\"" % options.region_file)
    region_list = read_region_bed_with_names(options.region_file)
    print(f"Read {len(region_list)} gene regions")
    
    with open(options.methylation_files, 'r') as f:
        methylation_files = [line.strip() for line in f if line.strip()]
    
    FOUT = open(options.out_file + '.txt', 'w')
    ROUT = open(options.out_file + '.r', 'w')
    
    print("Species_System\tCell_Type\tBin\tRegion_Type\tSpearman_Correlation\tP_Value\tSignificance\tData_Points", file=FOUT)
    
    print('.libPaths("/home/users/zzma/R/x86_64-pc-linux-gnu-library/4.3")', file=ROUT)
    print('library(RColorBrewer)', file=ROUT)
    print('pdf(file=\"%s\", width=10, height=6)' % (options.out_file + '.pdf'), file=ROUT)
    print('par(mar=c(5,5,4,2))', file=ROUT)
    
    cell_type_colors = {}
    species_line_types = {'Parent_Human': 1, 'Parent_Chimp': 2, 'Hybrid_Human': 3, 'Hybrid_Chimp': 4}
    
    unique_cell_types = design_df['cell_type'].unique()
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'grey']
    for i, cell_type in enumerate(unique_cell_types):
        cell_type_colors[cell_type] = colors[i % len(colors)]
    
    all_correlations = {}
    first_plot = True
    upstream_bins = downstream_bins = total_bins = None
    
    valid_files = 0
    for meth_file in methylation_files:
        printlog("Processing methylation file: \"%s\"" % meth_file)
        
        meth_basename = os.path.basename(meth_file)
        design_row = design_df[design_df['methylation_file'] == meth_basename]
        if design_row.empty:
            print(f"Warning: No design info found for {meth_file}", file=sys.stderr)
            continue
        
        design_info = design_row.iloc[0]
        expression_column = design_info['expression_column']
        
        if expression_column not in expression_df.columns:
            print(f"Warning: Expression column {expression_column} not found in expression data. Skipping {meth_file}", file=sys.stderr)
            continue
        
        try:
            cpg_ranges = read_CpG_bed_fixed(meth_file)
        except Exception as e:
            print(f"Warning: Could not read {meth_file}: {e}", file=sys.stderr)
            continue
        
        # Bin methylation data around TSS with fixed range
        binned_methylation, up_bins, down_bins, t_bins = bin_methylation_data_tss_fixed(
            region_list, cpg_ranges, options.upstream_bins, options.downstream_bins,
            options.upstream_size, options.downstream_size)
        
        if upstream_bins is None:
            upstream_bins, downstream_bins, total_bins = up_bins, down_bins, t_bins
        
        # Calc correlation and p
        correlations, p_vals, bin_counts = calculate_correlation_profile(binned_methylation, expression_df, expression_column)
        
        data_key = f"{design_info['species_system']}_{design_info['cell_type']}"
        data_key_safe = data_key.replace("-", "_").replace(" ", "_")
        
        all_correlations[data_key_safe] = {
            'correlations': correlations,
            'color': cell_type_colors[design_info['cell_type']],
            'line_type': species_line_types.get(design_info['species_system'], 1),
            'label': f"{design_info['species_system']} {design_info['cell_type']}"
        }
        
        for bin_idx in sorted(correlations.keys()):
            corr_val = correlations[bin_idx]
            p_val = p_vals[bin_idx]
            sig_marker = get_significance_marker(p_val)
            data_points = bin_counts[bin_idx]
            
            if bin_idx < upstream_bins:
                region_type = "Upstream"
            else:
                region_type = "Downstream"
            
            print(f"{design_info['species_system']}\t{design_info['cell_type']}\t{bin_idx}\t{region_type}\t{corr_val}\t{p_val}\t{sig_marker}\t{data_points}", file=FOUT)
        
        valid_files += 1
    
    print(f"Successfully processed {valid_files} methylation files")
if __name__ == '__main__':
    main()
