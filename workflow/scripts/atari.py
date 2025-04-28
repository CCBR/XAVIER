#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Tuple, Optional
from tqdm import tqdm
from rich.console import Console
from matplotlib.ticker import FuncFormatter

__version__="0.2.0"

console = Console()
ascii_banner = f"""
[cyan]
░█████╗░████████╗░█████╗░██████╗░██╗
██╔══██╗╚══██╔══╝██╔══██╗██╔══██╗██║
███████║░░░██║░░░███████║██████╔╝██║
██╔══██║░░░██║░░░██╔══██║██╔══██╗██║
██║░░██║░░░██║░░░██║░░██║██║░░██║██║
╚═╝░░╚═╝░░░╚═╝░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═╝
[/cyan]

[bold]A[/bold]lterna[bold]T[/bold]e [bold]A[/bold]llele [bold]R[/bold]ead v[bold]I[/bold]sualizer 

[magenta] [bold] Version:  {__version__} [/bold] [/magenta] 
"""

description=f"""
A tool to count reference and alternate alleles at genomic positions using samtools via pysam.
This script analyzes BAM files and reports the number of reads supporting reference and alternate
alleles at each position within specified regions, then visualizes the results.
"""


console.print(ascii_banner)
console.print(description)

usage_epilog=f"""
Examples:
  # Basic usage with BED file
  python allele_counter_visualizer.py \\
    --bam sample1.bam sample2.bam \\
    --reference genome.fa \\
    --bed regions.bed \\
    --window 100 \\
    --output results.tsv

  # Full pipeline with visualization and gene annotations
  python allele_counter_visualizer.py \\
    --bam sample1.bam sample2.bam \\
    --reference genome.fa \\
    --bed regions.bed \\
    --window 100 \\
    --output results.tsv \\
    --plot \\
    --plot-output coverage_plot.png \\
    --gtf annotations.gtf \\
    --verbose
    
  # Generate detailed base-level counts
  python allele_counter_visualizer.py \\
    --bam sample1.bam \\
    --reference genome.fa \\
    --bed regions.bed \\
    --window 100 \\
    --output detailed_results.tsv \\
    --detailed \\
    --verbose

  # Analyze a specific region with advanced visualization options
  python allele_counter_visualizer.py \\
    --bam sample.bam \\
    --reference genome.fa \\
    --chromosome chr20 \\
    --start 1000000 \\
    --end 1100000 \\
    --window 100 \\
    --output chr20_results.tsv \\
    --plot \\
    --plot-output chr20_coverage.png \\
    --gtf annotations.gtf \\
    --figsize 14,10 \\
    --dpi 600 \\
    --max-genes 15

  # Multi-sample comparison with quality filters
  python allele_counter_visualizer.py \\
    --bam tumor.bam normal.bam \\
    --reference genome.fa \\
    --bed hotspot_regions.bed \\
    --window 50 \\
    --min-mapping-quality 30 \\
    --min-base-quality 30 \\
    --min-depth 20 \\
    --output comparison_results.tsv \\
    --plot \\
    --plot-output sample_comparison.png
    """

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Count reference and alternate alleles from BAM files and visualize results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=usage_epilog)
    
    
    # Required arguments
    parser.add_argument('-b', '--bam', required=True, nargs='+',
                        help='Input BAM file(s). Must be indexed.')
    parser.add_argument('-r', '--reference', required=True,
                        help='Reference genome in FASTA format. Must be indexed.')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file name (TSV format).')
    parser.add_argument('-w', '--window', required=True, type=int, 
                    help='Window size (bp) for calculating mean ref/alt depths in non-overlapping windows.')
    
    # Region specification (mutually exclusive)
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument('-B', '--bed',
                            help='BED file with regions to analyze.')
    region_group.add_argument('-c', '--chromosome',
                            help='Chromosome to analyze (e.g., chr1).')
    
    # Optional parameters for manual region specification
    parser.add_argument('-s', '--start', type=int,
                        help='Start position (1-based, inclusive). Required if --chromosome is used.')
    parser.add_argument('-e', '--end', type=int,
                        help='End position (1-based, inclusive). Required if --chromosome is used.')
    
    # Optional filters
    parser.add_argument('-q', '--min-mapping-quality', type=int, default=20,
                        help='Minimum mapping quality score (default: 20).')
    parser.add_argument('-Q', '--min-base-quality', type=int, default=20,
                        help='Minimum base quality score (default: 20).')
    parser.add_argument('-d', '--min-depth', type=int, default=0,
                        help='Minimum read depth to report a position (default: 0).')
    
    # Output options
    parser.add_argument('-D', '--detailed', action='store_true',
                        help='Include detailed counts for each alternate base (A,C,G,T,DEL,INS).')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print progress information.')
    parser.add_argument('--no-progress-bar', action='store_true',
                        help='Disable progress bar.')
    
    # Visualization options
    visualization_group = parser.add_argument_group('Visualization Options')
    visualization_group.add_argument('--plot', action='store_true',
                        help='Generate a visualization plot of the results.')
    visualization_group.add_argument('--plot-output',
                        help='Output file for visualization (e.g., plot.png). Required if --plot is specified.')
    visualization_group.add_argument('--dpi', type=int, default=300,
                        help='DPI (resolution) for output image (default: 300).')
    visualization_group.add_argument('--figsize', type=str, default="12,8",
                        help='Figure size in inches, format "width,height" (default: "12,8").')
    visualization_group.add_argument('--no-show-plot', action='store_true',
                        help='Do not display plot (only save if --plot-output is specified).')
    visualization_group.add_argument('--gtf',
                        help='GTF file with gene annotations for visualization.')
    visualization_group.add_argument('--max-genes', type=int, default=10,
                        help='Maximum number of genes to display in each region (default: 10).')
    
    args = parser.parse_args()

    # Validate arguments
    if args.chromosome and (args.start is None or args.end is None):
        parser.error("--start and --end are required when --chromosome is specified.")
    
    if args.start is not None and args.end is not None and args.start > args.end:
        parser.error("Start position must be less than or equal to end position.")
    
    for bam_file in args.bam:
        if not os.path.exists(bam_file):
            parser.error(f"BAM file not found: {bam_file}")
        if not os.path.exists(f"{bam_file}.bai") and not os.path.exists(f"{os.path.splitext(bam_file)[0]}.bai"):
            parser.error(f"Index not found for BAM file: {bam_file}")
    
    if not os.path.exists(args.reference):
        parser.error(f"Reference file not found: {args.reference}")
    if not os.path.exists(f"{args.reference}.fai"):
        parser.error(f"Index not found for reference file: {args.reference}")
    
    if args.bed and not os.path.exists(args.bed):
        parser.error(f"BED file not found: {args.bed}")
    
    if args.plot and not args.plot_output and args.no_show_plot:
        parser.error("--plot-output is required when --plot and --no-show-plot are both specified.")
    
    if args.gtf and not os.path.exists(args.gtf):
        parser.error(f"GTF file not found: {args.gtf}")
    
    return args


def parse_bed_file(bed_file: str) -> List[Dict[str, any]]:
    """Parse regions from a BED file."""
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
                
            regions.append({
                'chrom': fields[0],
                'start': int(fields[1]) + 1,  # Convert from 0-based to 1-based
                'end': int(fields[2])
            })
    
    return regions


def get_total_regions_size(regions: List[Dict[str, any]]) -> int:
    """Calculate total size of all regions."""
    total_size = 0
    for region in regions:
        total_size += region['end'] - region['start'] + 1
    return total_size


def count_alleles(bam_file: str, reference_fasta: str, chromosome: str, start_pos: int, end_pos: int, 
                  min_mapping_quality: int = 20, min_base_quality: int = 20, 
                  progress_callback=None) -> List[Dict[str, any]]:
    """
    Count reference and alternate alleles at each position in a region.
    
    Args:
        bam_file: Path to the BAM file
        reference_fasta: Path to the reference genome
        chromosome: Chromosome name
        start_pos: Start position (1-based)
        end_pos: End position (1-based)
        min_mapping_quality: Minimum mapping quality
        min_base_quality: Minimum base quality
        progress_callback: Function to call for progress updates with position as argument
    """
    # Open files
    samfile = pysam.AlignmentFile(bam_file, "rb")
    fastafile = pysam.FastaFile(reference_fasta)
    results = []
    
    # Check if chromosome exists in both BAM and reference
    try:
        test_fetch = samfile.fetch(chromosome, start_pos-1, start_pos)
        next(test_fetch, None)
    except (ValueError, StopIteration):
        return results
    
    try:
        fastafile.fetch(chromosome, start_pos-1, start_pos)
    except (ValueError, KeyError):
        return results
    
    # Process each position
    for pos in range(start_pos, end_pos + 1):
        # Call progress callback for each position
        if progress_callback:
            progress_callback()
            
        # Get reference base
        try:
            ref_base = fastafile.fetch(chromosome, pos-1, pos).upper()
        except (ValueError, KeyError):
            continue
        
        # Initialize counters
        base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'DEL': 0, 'INS': 0}
        coverage = 0
        
        # Pileup at the current position
        for pileupcolumn in samfile.pileup(chromosome, pos-1, pos, truncate=True, min_base_quality=min_base_quality):
            if pileupcolumn.pos != pos-1:
                continue
                
            coverage = pileupcolumn.n
            
            # Process each read
            for pileupread in pileupcolumn.pileups:
                # Skip low quality reads
                if pileupread.alignment.mapping_quality < min_mapping_quality:
                    continue
                
                # Count deletions
                if pileupread.is_del:
                    base_counts['DEL'] += 1
                    continue
                
                # Count insertions
                if pileupread.indel > 0:
                    base_counts['INS'] += 1
                
                # Get base
                if pileupread.query_position is None:
                    continue
                    
                base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                
                # Count base
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1
        
        # Calculate reference and alternate allele counts
        ref_count = base_counts.get(ref_base, 0)
        total_alt_count = sum(count for base, count in base_counts.items() 
                             if base != ref_base and base != 'N' and count > 0)
        
        # Add result
        result = {
            'Bam_file': os.path.basename(bam_file),
            'CHROM': chromosome,
            'POS': pos,
            'REF': ref_base,
            'DP': coverage,
            'Ref_reads': ref_count,
            'Alt_reads': total_alt_count,
            'A_count': base_counts['A'],
            'C_count': base_counts['C'],
            'G_count': base_counts['G'],
            'T_count': base_counts['T'],
            'DEL_count': base_counts['DEL'],
            'INS_count': base_counts['INS'],
            'N_count': base_counts['N'],
        }
        
        results.append(result)
    
    # Close files
    samfile.close()
    fastafile.close()
    
    return results


def process_regions(bam_files: List[str], reference_fasta: str, regions: List[Dict[str, any]], 
                    min_mapping_quality: int, min_base_quality: int, min_depth: int, 
                    verbose: bool, show_progress: bool) -> pd.DataFrame:
    """Process all specified regions across all BAM files."""
    all_results = []
    
    # Calculate the total work - number of positions to process
    total_positions = get_total_regions_size(regions) * len(bam_files)
    positions_processed = 0
    
    # Create progress bar based on total positions
    progress_bar = None
    if show_progress:
        progress_bar = tqdm(
            total=total_positions,
            desc="Processing genomic positions",
            unit="pos",
            bar_format="{desc}: {percentage:3.1f}% |{bar}| {n_fmt}/{total_fmt} positions"
        )
    
    # Process each BAM file and region
    for i, bam_file in enumerate(bam_files, 1):
        if verbose:
            print(f"Processing BAM file {i}/{len(bam_files)}: {bam_file}", file=sys.stderr)
        
        for j, region in enumerate(regions, 1):
            if verbose and not show_progress:
                print(f"  Region {j}/{len(regions)}: {region['chrom']}:{region['start']}-{region['end']}", file=sys.stderr)
            
            # Update progress function
            def update_progress():
                nonlocal positions_processed
                positions_processed += 1
                if progress_bar:
                    progress_bar.update(1)
            
            # Process the region
            results = count_alleles(
                bam_file, 
                reference_fasta, 
                region['chrom'], 
                region['start'], 
                region['end'], 
                min_mapping_quality, 
                min_base_quality,
                update_progress
            )
            
            all_results.extend(results)
    
    # Close progress bar
    if progress_bar:
        progress_bar.close()
    
    # Convert results to DataFrame
    if not all_results:
        return pd.DataFrame()
        
    df = pd.DataFrame(all_results)
    
    # Filter by minimum depth
    if min_depth > 0:
        df = df[df['DP'] >= min_depth]
    
    return df

def process_windows(df: pd.DataFrame, window_size: int) -> pd.DataFrame:
    """
    Group results into non-overlapping windows and calculate mean values.
    
    Args:
        df: DataFrame with position-level results
        window_size: Size of the windows in base pairs
    
    Returns:
        DataFrame with window-level summaries
    """
    if df.empty:
        return df
    
    # Add a window ID column based on position
    df['Window'] = (df['POS'] - 1) // window_size
    
    # Group by BAM file, chromosome, and window
    window_groups = df.groupby(['Bam_file', 'CHROM', 'Window'])
    
    # Calculate window summaries
    window_results = []
    
    for (bam, chrom, window), group in window_groups:
        # Calculate window bounds
        start_pos = window * window_size + 1
        end_pos = (window + 1) * window_size
        
        # Calculate means
        mean_ref_reads = group['Ref_reads'].mean()
        mean_alt_reads = group['Alt_reads'].mean()
        mean_depth = group['DP'].mean()
        
        # Add detailed base counts if they exist
        result = {
            'Bam_file': bam,
            'CHROM': chrom,
            'Window_start': start_pos,
            'Window_end': end_pos,
            'Mean_depth': mean_depth,
            'Mean_ref_reads': mean_ref_reads,
            'Mean_alt_reads': mean_alt_reads,
            'Positions_count': len(group)
        }
        
        # Include mean base counts if detailed output is requested
        if 'A_count' in group.columns:
            result.update({
                'Mean_A_count': group['A_count'].mean(),
                'Mean_C_count': group['C_count'].mean(),
                'Mean_G_count': group['G_count'].mean(),
                'Mean_T_count': group['T_count'].mean(),
                'Mean_DEL_count': group['DEL_count'].mean(),
                'Mean_INS_count': group['INS_count'].mean(),
                'Mean_N_count': group['N_count'].mean(),
            })
        
        window_results.append(result)
    
    return pd.DataFrame(window_results)

def parse_gtf_file(gtf_file: str) -> Dict[str, List[Dict]]:
    """Parse gene and exon information from a GTF file.
    
    Args:
        gtf_file: Path to the GTF file
        
    Returns:
        Dictionary mapping chromosome names to lists of features
    """
    features_by_chrom = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Parse GTF line
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])  # GTF is 1-based
            end = int(fields[4])
            strand = fields[6]
            
            # Only process genes and exons
            if feature_type not in ['gene', 'exon', 'transcript', 'CDS']:
                continue
            
            # Parse attributes
            attr_str = fields[8]
            attrs = {}
            
            # Extract key attributes
            for attr_pair in attr_str.split(';'):
                attr_pair = attr_pair.strip()
                if not attr_pair:
                    continue
                
                # Handle different GTF formats (Ensembl vs UCSC)
                if '=' in attr_pair:
                    key, value = attr_pair.split('=', 1)
                    value = value.strip('"')
                elif ' ' in attr_pair:
                    key, value = attr_pair.split(' ', 1)
                    value = value.strip('"')
                else:
                    continue
                
                attrs[key] = value
            
            # Get gene name and ID
            gene_name = attrs.get('gene_name', attrs.get('gene_id', 'Unknown'))
            gene_id = attrs.get('gene_id', 'Unknown')
            
            # Store feature information
            feature = {
                'type': feature_type,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_name': gene_name,
                'gene_id': gene_id
            }
            
            # Add transcript and exon IDs if available
            if feature_type in ['transcript', 'exon', 'CDS']:
                feature['transcript_id'] = attrs.get('transcript_id', 'Unknown')
                
            # Add to chromosome dictionary
            if chrom not in features_by_chrom:
                features_by_chrom[chrom] = []
            
            features_by_chrom[chrom].append(feature)
    
    return features_by_chrom

def plot_gene_annotations(ax, features, start_pos, end_pos, max_genes=10):
    """Plot gene annotations on the given axes.
    
    Args:
        ax: Matplotlib axes to plot on
        features: List of gene features for the current chromosome
        start_pos: Start position for the current view
        end_pos: End position for the current view
        max_genes: Maximum number of genes to display
        
    Returns:
        Boolean indicating if a legend was added
    """
    # Filter features for the visible region
    visible_genes = {}
    visible_transcripts = {}
    visible_exons = {}
    
    for feature in features:
        # Skip features outside the visible region
        if feature['end'] < start_pos or feature['start'] > end_pos:
            continue
            
        if feature['type'] == 'gene':
            gene_id = feature['gene_id']
            visible_genes[gene_id] = feature
            
        elif feature['type'] == 'transcript':
            transcript_id = feature['transcript_id']
            gene_id = feature['gene_id']
            
            if gene_id not in visible_transcripts:
                visible_transcripts[gene_id] = []
                
            visible_transcripts[gene_id].append(feature)
            
        elif feature['type'] in ['exon', 'CDS']:
            transcript_id = feature['transcript_id']
            gene_id = feature['gene_id']
            
            if gene_id not in visible_exons:
                visible_exons[gene_id] = []
                
            visible_exons[gene_id].append(feature)
    
    # Limit the number of genes if there are too many
    gene_ids = list(visible_genes.keys())
    if len(gene_ids) > max_genes:
        # Sort genes by length (show longest genes)
        gene_ids.sort(key=lambda g: visible_genes[g]['end'] - visible_genes[g]['start'], reverse=True)
        gene_ids = gene_ids[:max_genes]
        
    # No genes to display
    if not gene_ids:
        ax.text(0.5, 0.5, 'No gene annotations in view', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_yticks([])
        return False
        
    # Plot genes
    y_pos = 0
    y_positions = {}
    colors = plt.cm.tab10.colors
    
    # For legend - collect items to add to the legend
    legend_elements = []
    
    for i, gene_id in enumerate(gene_ids):
        gene = visible_genes[gene_id]
        color = colors[i % len(colors)]
        y_pos = i
        y_positions[gene_id] = y_pos
        
        # Plot gene as a thick line
        gene_start = max(gene['start'], start_pos)
        gene_end = min(gene['end'], end_pos)
        
        # Gene body
        gene_line = ax.plot([gene_start, gene_end], [y_pos, y_pos], 
                color=color, linewidth=2, solid_capstyle='butt')[0]
        
        # Direction arrow
        if gene['strand'] == '+':
            ax.arrow(gene_end - (gene_end - gene_start) * 0.05, y_pos, 
                    (gene_end - gene_start) * 0.03, 0, 
                    head_width=0.2, head_length=(gene_end - gene_start) * 0.02, 
                    fc=color, ec=color)
        else:
            ax.arrow(gene_start + (gene_end - gene_start) * 0.05, y_pos, 
                    -(gene_end - gene_start) * 0.03, 0, 
                    head_width=0.2, head_length=(gene_end - gene_start) * 0.02, 
                    fc=color, ec=color)
        
        # Add to legend items for this gene
        if i == 0:  # Only add one gene example to the legend
            legend_elements.append(plt.Line2D([0], [0], color=color, lw=2, 
                                            label='Gene body'))
            # Add arrow style to legend 
            legend_elements.append(plt.Line2D([0], [0], marker='>',
                                         color='w', markerfacecolor=color, 
                                         markersize=8, label='Direction'))
        
        exon_legend_added = False
        cds_legend_added = False
        
        # Plot exons if available
        if gene_id in visible_exons:
            exons = visible_exons[gene_id]
            for exon in exons:
                exon_start = max(exon['start'], start_pos)
                exon_end = min(exon['end'], end_pos)
                
                if exon['type'] == 'CDS':
                    # CDS as taller box
                    cds_patch = plt.Rectangle(
                        (exon_start, y_pos - 0.3), 
                        exon_end - exon_start, 
                        0.6, 
                        facecolor=color, 
                        edgecolor='none', 
                        alpha=0.8)
                    ax.add_patch(cds_patch)
                    
                    # Add to legend once
                    if i == 0 and not cds_legend_added:
                        legend_elements.append(plt.Rectangle((0, 0), 1, 1, 
                                              facecolor=color, alpha=0.8, 
                                              label='CDS (coding)'))
                        cds_legend_added = True
                else:
                    # Regular exon as box
                    exon_patch = plt.Rectangle(
                        (exon_start, y_pos - 0.2), 
                        exon_end - exon_start, 
                        0.4, 
                        facecolor=color, 
                        edgecolor='none', 
                        alpha=0.6)
                    ax.add_patch(exon_patch)
                    
                    # Add to legend once
                    if i == 0 and not exon_legend_added:
                        legend_elements.append(plt.Rectangle((0, 0), 1, 1, 
                                              facecolor=color, alpha=0.6, 
                                              label='Exon'))
                        exon_legend_added = True
        
        # Add gene name
        ax.text(gene_start, y_pos + 0.4, gene['gene_name'], 
               fontsize=8, ha='left', va='bottom', 
               bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0))
    
    # Set y-axis limits and labels
    ax.set_ylim(-1, len(gene_ids))
    ax.set_yticks([])
    
    # Format x-axis in Mbp
    def format_mbp(x, pos):
        return f'{x/1_000_000:.2f}'
    
    ax.xaxis.set_major_formatter(FuncFormatter(format_mbp))
    ax.set_xlabel('Position (Mbp)')
    
    # Add legend for gene features
    if legend_elements:
        ax.legend(handles=legend_elements, 
                 loc='upper right', 
                 fontsize='small',
                 framealpha=0.7,
                 title="Gene Features")
        return True
    
    return False

def plot_allele_counts_with_annotations(df, gtf_features=None, output_file=None, dpi=300, 
                                      figsize=(12, 8), show_plot=True, max_genes=10,
                                      min_coverage_threshold=None, highlight_regions=None):
    """
    Generate an improved visualization of allele counts with gene annotations.
    
    Args:
        df: DataFrame with allele count data
        gtf_features: Dictionary of gene features by chromosome from parse_gtf_file
        output_file: Path to save the plot (if None, will display plot)
        dpi: Resolution for the output image
        figsize: Figure size as a tuple (width, height) in inches
        show_plot: Whether to display the plot (in addition to saving it)
        max_genes: Maximum number of genes to display per region
        min_coverage_threshold: Optional minimum coverage threshold to display as horizontal line
        highlight_regions: Optional list of dict with 'start', 'end', and 'label' keys for regions to highlight
    """
    # Check if the required columns exist
    required_cols = ['CHROM', 'Window_start', 'Window_end', 'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads']
    if not all(col in df.columns for col in required_cols):
        print(f"DataFrame must contain columns: {', '.join(required_cols)}", file=sys.stderr)
        return
    
    # Get unique BAM files and chromosomes
    bam_files = df['Bam_file'].unique()
    chromosomes = df['CHROM'].unique()
    
    # Set up the plot
    num_bams = len(bam_files)
    num_chroms = len(chromosomes)
    
    # Create a custom color palette with more distinct colors
    palette = {
        'depth': '#1f77b4',  # Blue
        'ref': '#2ca02c',    # Green
        'alt': '#d62728',    # Red
        'ratio': '#9467bd'   # Purple
    }
    
    # Adjust layout to include ONE gene annotation panel per chromosome
    # Modified to have a single gene annotation panel per chromosome
    if num_bams == 1 and num_chroms == 1:
        fig = plt.figure(figsize=figsize)
        gs = plt.GridSpec(2 + num_bams, 1, height_ratios=[3] * num_bams + [1, 1], figure=fig)
        
        # Coverage plot
        ax_cov = fig.add_subplot(gs[0])  
        # Ratio plot
        ax_ratio = fig.add_subplot(gs[1], sharex=ax_cov)
        # Single gene annotation plot
        ax_genes = fig.add_subplot(gs[2], sharex=ax_cov)
        
        axes_cov = [[ax_cov]]
        axes_ratio = [[ax_ratio]]
        axes_genes = [[ax_genes]]  # One genes panel
    elif num_bams >= 1 and num_chroms == 1:
        # One chromosome, multiple BAMs
        fig = plt.figure(figsize=(figsize[0], figsize[1] * (num_bams + 0.5)))
        # Each BAM gets a coverage and ratio plot, plus one shared gene plot at the bottom
        gs = plt.GridSpec(2 * num_bams + 1, 1, 
                         height_ratios=[3, 1] * num_bams + [1], 
                         figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            # Coverage plot for this BAM
            ax_cov = fig.add_subplot(gs[2*i])
            # Ratio plot for this BAM
            ax_ratio = fig.add_subplot(gs[2*i+1], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
        
        # Single gene annotation at the bottom for this chromosome
        ax_genes = fig.add_subplot(gs[2*num_bams], sharex=axes_cov[0][0])
        # Store in a 2D array for consistency
        axes_genes = [[ax_genes]]
    elif num_bams == 1 and num_chroms > 1:
        # Multiple chromosomes, one BAM
        fig = plt.figure(figsize=(figsize[0], figsize[1] * num_chroms))
        gs = plt.GridSpec(3 * num_chroms, 1, height_ratios=[3, 1, 1] * num_chroms, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        axes_genes = []
        
        for j in range(num_chroms):
            # Coverage plot
            ax_cov = fig.add_subplot(gs[3*j])
            # Ratio plot
            ax_ratio = fig.add_subplot(gs[3*j+1], sharex=ax_cov)
            # Gene annotation for this chromosome
            ax_genes = fig.add_subplot(gs[3*j+2], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
            axes_genes.append([ax_genes])
    else:
        # Multiple BAMs and chromosomes
        fig = plt.figure(figsize=(figsize[0] * num_chroms, figsize[1] * (num_bams + 0.5)))
        
        # Calculate heights: each BAM gets a coverage and ratio plot, 
        # plus one gene plot per chromosome
        gs = plt.GridSpec(2 * num_bams + 1, num_chroms, 
                         height_ratios=[3, 1] * num_bams + [1], 
                         figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            axes_cov_row = []
            axes_ratio_row = []
            
            for j in range(num_chroms):
                # Coverage plot
                ax_cov = fig.add_subplot(gs[2*i, j])
                # Ratio plot
                ax_ratio = fig.add_subplot(gs[2*i+1, j], sharex=ax_cov)
                
                axes_cov_row.append(ax_cov)
                axes_ratio_row.append(ax_ratio)
            
            axes_cov.append(axes_cov_row)
            axes_ratio.append(axes_ratio_row)
        
        # Create one gene annotation row per chromosome at the bottom
        axes_genes = [[fig.add_subplot(gs[2*num_bams, j], sharex=axes_cov[0][j])] 
                      for j in range(num_chroms)]
    
    # Format function for x-axis (millions of base pairs)
    def format_mbp(x, pos):
        return f'{x/1_000_000:.2f}'
    
    # Process each chromosome to determine plot ranges
    chrom_data = {}
    for j, chrom in enumerate(chromosomes):
        # Get all data for this chromosome across all BAMs
        chrom_subset = df[df['CHROM'] == chrom]
        if chrom_subset.empty:
            continue
            
        x_min = chrom_subset['Window_start'].min()
        x_max = chrom_subset['Window_end'].max()
        chrom_data[chrom] = {
            'x_min': x_min,
            'x_max': x_max
        }
    
    # Create each subplot
    for i, bam in enumerate(bam_files):
        for j, chrom in enumerate(chromosomes):
            subset = df[(df['Bam_file'] == bam) & (df['CHROM'] == chrom)]
            
            if subset.empty:
                continue
                
            # Calculate window midpoints for x-axis
            subset['window_mid'] = (subset['Window_start'] + subset['Window_end']) / 2
            
            # Calculate alt/total ratio
            subset['alt_ratio'] = (subset['Mean_alt_reads'] / (subset['Mean_depth'].replace(0, np.nan))) * 100  # Convert to percentage
            subset['alt_ratio'] = subset['alt_ratio'].fillna(0)  # Replace NaNs with 0
            
            # Get current axes
            ax_cov = axes_cov[i][j] if num_bams > 1 or num_chroms > 1 else axes_cov[0][0]
            ax_ratio = axes_ratio[i][j] if num_bams > 1 or num_chroms > 1 else axes_ratio[0][0]
            # For genes, use the chromosome index only since we share gene plots across BAMs
            ax_genes = axes_genes[j][0] if num_chroms > 1 else axes_genes[0][0]
            
            # Get axis boundaries for highlighting
            x_min = subset['Window_start'].min()
            x_max = subset['Window_end'].max()
            
            # Plot coverage data with enhanced styling
            # Add shaded background for alt reads area
            ax_cov.set_facecolor('#f9f9f9')  # Light background
            ax_cov.axvspan(x_min, x_max, color='#fff3f3', alpha=0.5)  # Light red shading for alt reads area
            
            # Add grid for better readability
            ax_cov.grid(True, linestyle='--', alpha=0.3)
            
            # Draw coverage threshold if specified
            if min_coverage_threshold is not None:
                ax_cov.axhline(y=min_coverage_threshold, color='#aaaaaa', linestyle='--', alpha=0.7, 
                              label=f'Min. Coverage ({min_coverage_threshold}x)')
            
            # Plot main coverage data with both lines and markers
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], '-', 
                       color=palette['depth'], label='Mean Depth', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], 'o', 
                       color=palette['depth'], markersize=4, alpha=0.5)
            
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], '-', 
                       color=palette['ref'], label='Ref Reads', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], 's', 
                       color=palette['ref'], markersize=4, alpha=0.5)
            
            # Calculate reasonable y-max (max depth + 10%)
            y_max = max(subset['Mean_depth'].max(), subset['Mean_ref_reads'].max()) * 1.1
            
            # Add alternate allele subplot with different scale and clearer styling
            alt_color = palette['alt']
            ax2 = ax_cov.twinx()
            # Add light red background to make it clear this is a different axis
            ax2.patch.set_alpha(0.0)  # Make background transparent
            
            # Add lines and markers for alt reads
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '-', 
                    color=alt_color, linewidth=1.5, label='Alt Reads')
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '^', 
                    color=alt_color, markersize=4, alpha=0.5)
            
            # Set alt reads scale
            alt_max = max(subset['Mean_alt_reads'].max() * 1.5, 5)  # Ensure non-zero max
            ax2.set_ylim(0, alt_max)
            
            # Styling for alt reads axis - Fix for Issue 2: No box around the label
            alt_label = 'Alternate Reads (Right Axis)'
            ax2.set_ylabel(alt_label, color=alt_color, fontweight='bold')
            ax2.tick_params(axis='y', labelcolor=alt_color)
            
            # Add annotation to clarify alt reads is on right axis
            ax_cov.annotate('Alt Reads →', xy=(0.98, 0.5), xycoords='axes fraction',
                          ha='right', va='center', color=alt_color, fontsize=9)
            
            # Plot ratio data
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], '-', 
                        color=palette['ratio'], label='Alt/Total Ratio', linewidth=1.5)
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], 'o', 
                        color=palette['ratio'], markersize=4, alpha=0.5)
            ax_ratio.set_ylim(0, min(100.0, subset['alt_ratio'].max() * 1.5 or 10.0))
            ax_ratio.set_ylabel('Alt/Total Ratio (%)')
            ax_ratio.grid(True, linestyle='--', alpha=0.3)
            
            # Highlight regions of interest if specified
            if highlight_regions:
                for region in highlight_regions:
                    if region['start'] <= x_max and region['end'] >= x_min:  # Only if overlaps visible region
                        # Add transparent highlighting to all subplots
                        for ax in [ax_cov, ax_ratio]:
                            ax.axvspan(region['start'], region['end'], 
                                     color='yellow', alpha=0.2, label=region.get('label', 'ROI'))
                        
                        # Add label above the highlighted region
                        if 'label' in region:
                            ax_cov.annotate(region['label'], 
                                          xy=((region['start'] + region['end'])/2, 0.95),
                                          xycoords=('data', 'axes fraction'),
                                          ha='center', va='center',
                                          bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3))
            
            # Set titles and labels for coverage plot
            ax_cov.set_title(f"{os.path.basename(bam)} - {chrom}", fontweight='bold')
            ax_cov.set_ylabel('Mean Read Count (Left Axis)', color='black', fontweight='bold')
            ax_cov.set_ylim(0, y_max)
            
            # Remove x-label from coverage and ratio plots
            ax_cov.set_xlabel('')
            ax_cov.tick_params(axis='x', labelbottom=False)
            ax_ratio.set_xlabel('')
            ax_ratio.tick_params(axis='x', labelbottom=False)
            
            # Add legend to first plot only with better positioning
            if i == 0 and j == 0:
                lines1, labels1 = ax_cov.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines3, labels3 = ax_ratio.get_legend_handles_labels()
                # Use a more visible location for the legend
                leg = ax_cov.legend(lines1 + lines2 + lines3, 
                                  labels1 + labels2 + labels3, 
                                  loc='upper right', 
                                  framealpha=0.9, 
                                  fancybox=True,
                                  shadow=True)
                leg.get_frame().set_edgecolor('#888888')
    
    # Plot gene annotations at the bottom, once per chromosome, after all data plots
    for j, chrom in enumerate(chromosomes):
        if chrom not in chrom_data:
            continue
            
        # Get gene annotation axis - one per chromosome
        ax_genes = axes_genes[j][0] if num_chroms > 1 else axes_genes[0][0]
        
        # Use chromosome bounds for annotations
        min_pos = chrom_data[chrom]['x_min']
        max_pos = chrom_data[chrom]['x_max']
        
        # Add gene annotations if available
        if gtf_features and chrom in gtf_features:
            gene_legend_added = plot_gene_annotations(ax_genes, gtf_features[chrom], min_pos, max_pos, max_genes)
            # Format x-axis in Mbp
            ax_genes.xaxis.set_major_formatter(FuncFormatter(format_mbp))
            ax_genes.set_xlabel('Position (Mbp)')
        else:
            ax_genes.text(0.5, 0.5, 'No gene annotations available', 
                         ha='center', va='center', transform=ax_genes.transAxes)
            ax_genes.set_yticks([])
            ax_genes.set_xlabel('Position (bp)')
    
    plt.tight_layout()
    
    # Save the plot if output file is specified
    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Plot saved to {output_file}", file=sys.stderr)
    
    # Display the plot if requested
    if show_plot:
        plt.show()
    else:
        plt.close()

def plot_allele_counts(df, output_file=None, dpi=300, figsize=(12, 8), show_plot=True,
                              min_coverage_threshold=None, highlight_regions=None):
    """
    Generate an improved visualization of allele counts with better visual clarity.
    
    Args:
        df: DataFrame with allele count data
        output_file: Path to save the plot (if None, will display plot)
        dpi: Resolution for the output image
        figsize: Figure size as a tuple (width, height) in inches
        show_plot: Whether to display the plot (in addition to saving it) 
        min_coverage_threshold: Optional minimum coverage threshold to display as horizontal line
        highlight_regions: Optional list of dict with 'start', 'end', and 'label' keys for regions to highlight
    """
    # Check if the required columns exist
    required_cols = ['CHROM', 'Window_start', 'Window_end', 'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads']
    if not all(col in df.columns for col in required_cols):
        print(f"DataFrame must contain columns: {', '.join(required_cols)}", file=sys.stderr)
        return
    
    # Get unique BAM files and chromosomes
    bam_files = df['Bam_file'].unique()
    chromosomes = df['CHROM'].unique()
    
    # Set up the plot
    num_bams = len(bam_files)
    num_chroms = len(chromosomes)
    
    # Create a custom color palette with more distinct colors
    palette = {
        'depth': '#1f77b4',  # Blue
        'ref': '#2ca02c',    # Green
        'alt': '#d62728',    # Red
        'ratio': '#9467bd'   # Purple
    }
    
    # Adjust layout to include ratio subplot
    if num_bams == 1 and num_chroms == 1:
        fig = plt.figure(figsize=figsize)
        gs = plt.GridSpec(2, 1, height_ratios=[3, 1], figure=fig)
        
        # Coverage plot
        ax_cov = fig.add_subplot(gs[0])
        # Ratio plot
        ax_ratio = fig.add_subplot(gs[1], sharex=ax_cov)
        
        axes_cov = [[ax_cov]]
        axes_ratio = [[ax_ratio]]
    elif num_bams == 1:
        fig = plt.figure(figsize=(figsize[0], figsize[1]*num_chroms))
        gs = plt.GridSpec(2*num_chroms, 1, height_ratios=[3, 1] * num_chroms, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for j in range(num_chroms):
            # Coverage plot
            ax_cov = fig.add_subplot(gs[2*j])
            # Ratio plot
            ax_ratio = fig.add_subplot(gs[2*j+1], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
    elif num_chroms == 1:
        fig = plt.figure(figsize=(figsize[0], figsize[1]*num_bams))
        gs = plt.GridSpec(2*num_bams, 1, height_ratios=[3, 1] * num_bams, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            # Coverage plot
            ax_cov = fig.add_subplot(gs[2*i])
            # Ratio plot
            ax_ratio = fig.add_subplot(gs[2*i+1], sharex=ax_cov)
            
            axes_cov.append([ax_cov])
            axes_ratio.append([ax_ratio])
    else:
        fig = plt.figure(figsize=(figsize[0]*num_chroms/2, figsize[1]*num_bams))
        gs = plt.GridSpec(2*num_bams, num_chroms, height_ratios=[3, 1] * num_bams, figure=fig)
        
        axes_cov = []
        axes_ratio = []
        
        for i in range(num_bams):
            axes_cov_row = []
            axes_ratio_row = []
            
            for j in range(num_chroms):
                # Coverage plot
                ax_cov = fig.add_subplot(gs[2*i, j])
                # Ratio plot
                ax_ratio = fig.add_subplot(gs[2*i+1, j], sharex=ax_cov)
                
                axes_cov_row.append(ax_cov)
                axes_ratio_row.append(ax_ratio)
            
            axes_cov.append(axes_cov_row)
            axes_ratio.append(axes_ratio_row)
    
    # Format function for x-axis (millions of base pairs)
    def format_mbp(x, pos):
        return f'{x/1_000_000:.2f}'
    
    # Create each subplot
    for i, bam in enumerate(bam_files):
        for j, chrom in enumerate(chromosomes):
            subset = df[(df['Bam_file'] == bam) & (df['CHROM'] == chrom)]
            
            if subset.empty:
                continue
                
            # Calculate window midpoints for x-axis
            subset['window_mid'] = (subset['Window_start'] + subset['Window_end']) / 2
            
            # Calculate alt/total ratio
            subset['alt_ratio'] = (subset['Mean_alt_reads'] / (subset['Mean_depth'].replace(0, np.nan))) * 100
            subset['alt_ratio'] = subset['alt_ratio'].fillna(0)  # Replace NaNs with 0
            
            # Get current axes
            ax_cov = axes_cov[i][j] if num_bams > 1 or num_chroms > 1 else axes_cov[0][0]
            ax_ratio = axes_ratio[i][j] if num_bams > 1 or num_chroms > 1 else axes_ratio[0][0]
            
            # Get axis boundaries for highlighting
            x_min = subset['Window_start'].min()
            x_max = subset['Window_end'].max()
            
            # Plot coverage data with enhanced styling
            # Add shaded background for alt reads area
            ax_cov.set_facecolor('#f9f9f9')  # Light background
            ax_cov.axvspan(x_min, x_max, color='#fff3f3', alpha=0.5)  # Light red shading for alt reads area
            
            # Add grid for better readability
            ax_cov.grid(True, linestyle='--', alpha=0.3)
            
            # Draw coverage threshold if specified
            if min_coverage_threshold is not None:
                ax_cov.axhline(y=min_coverage_threshold, color='#aaaaaa', linestyle='--', alpha=0.7, 
                              label=f'Min. Coverage ({min_coverage_threshold}x)')
            
            # Plot main coverage data with both lines and markers
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], '-', 
                       color=palette['depth'], label='Mean Depth', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_depth'], 'o', 
                       color=palette['depth'], markersize=4, alpha=0.5)
            
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], '-', 
                       color=palette['ref'], label='Ref Reads', linewidth=1.5)
            ax_cov.plot(subset['window_mid'], subset['Mean_ref_reads'], 's', 
                       color=palette['ref'], markersize=4, alpha=0.5)
            
            # Calculate reasonable y-max (max depth + 10%)
            y_max = max(subset['Mean_depth'].max(), subset['Mean_ref_reads'].max()) * 1.1
            
            # Add alternate allele subplot with different scale and clearer styling
            alt_color = palette['alt']
            ax2 = ax_cov.twinx()
            # Add light red background to make it clear this is a different axis
            ax2.patch.set_alpha(0.0)  # Make background transparent
            
            # Add lines and markers for alt reads
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '-', 
                    color=alt_color, linewidth=1.5, label='Alt Reads')
            ax2.plot(subset['window_mid'], subset['Mean_alt_reads'], '^', 
                    color=alt_color, markersize=4, alpha=0.5)
            
            # Set alt reads scale
            alt_max = max(subset['Mean_alt_reads'].max() * 1.5, 5)  # Ensure non-zero max
            ax2.set_ylim(0, alt_max)
            
            # Styling for alt reads axis
            alt_label = 'Alternate Reads (Right Axis)'
            ax2.set_ylabel(alt_label, color=alt_color, fontweight='bold')
            ax2.tick_params(axis='y', labelcolor=alt_color)
            
            # Draw a border around the right axis label to make it stand out
            ax2.yaxis.label.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor=alt_color, pad=5, boxstyle='round'))
            
            # Add annotation to clarify alt reads is on right axis
            ax_cov.annotate('Alt Reads →', xy=(0.98, 0.5), xycoords='axes fraction',
                          ha='right', va='center', color=alt_color, fontsize=9,
                          bbox=dict(facecolor='white', alpha=0.8, edgecolor=alt_color))
            
            # Plot ratio data
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], '-', 
                        color=palette['ratio'], label='Alt/Total Ratio', linewidth=1.5)
            ax_ratio.plot(subset['window_mid'], subset['alt_ratio'], 'o', 
                        color=palette['ratio'], markersize=4, alpha=0.5)
            ax_ratio.set_ylim(0, min(100.0, subset['alt_ratio'].max() * 1.5 or 10.0))
            ax_ratio.set_ylabel('Alt/Total Ratio')
            ax_ratio.grid(True, linestyle='--', alpha=0.3)
            
            # Highlight regions of interest if specified
            if highlight_regions:
                for region in highlight_regions:
                    if region['start'] <= x_max and region['end'] >= x_min:  # Only if overlaps visible region
                        # Add transparent highlighting to all subplots
                        for ax in [ax_cov, ax_ratio]:
                            ax.axvspan(region['start'], region['end'], 
                                     color='yellow', alpha=0.2, label=region.get('label', 'ROI'))
                        
                        # Add label above the highlighted region
                        if 'label' in region:
                            ax_cov.annotate(region['label'], 
                                          xy=((region['start'] + region['end'])/2, 0.95),
                                          xycoords=('data', 'axes fraction'),
                                          ha='center', va='center',
                                          bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3))
            
            # Set titles and labels for coverage plot
            ax_cov.set_title(f"{os.path.basename(bam)} - {chrom}", fontweight='bold')
            ax_cov.set_ylabel('Mean Read Count (Left Axis)', color='black', fontweight='bold')
            ax_cov.set_ylim(0, y_max)
            
            # Remove x-label from coverage plot
            ax_cov.set_xlabel('')
            ax_cov.tick_params(axis='x', labelbottom=False)
            
            # Format x-axis in Mbp for ratio plot
            ax_ratio.xaxis.set_major_formatter(FuncFormatter(format_mbp))
            ax_ratio.set_xlabel('Position (Mbp)')
            
            # Add legend to first plot only with better positioning
            if i == 0 and j == 0:
                lines1, labels1 = ax_cov.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                lines3, labels3 = ax_ratio.get_legend_handles_labels()
                # Use a more visible location for the legend
                leg = ax_cov.legend(lines1 + lines2 + lines3, 
                                  labels1 + labels2 + labels3, 
                                  loc='upper right', 
                                  framealpha=0.9, 
                                  fancybox=True,
                                  shadow=True)
                leg.get_frame().set_edgecolor('#888888')
    
    plt.tight_layout()
    
    # Save the plot if output file is specified
    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        print(f"Plot saved to {output_file}", file=sys.stderr)
    
    # Display the plot if requested
    if show_plot:
        plt.show()
    else:
        plt.close()


def main():
    """Main function."""
    args = parse_arguments()
    
    # Parse regions
    if args.bed:
        regions = parse_bed_file(args.bed)
    else:
        regions = [{
            'chrom': args.chromosome,
            'start': args.start,
            'end': args.end
        }]
    
    if args.verbose:
        console.print(f"Reference genome: [bold]{args.reference}[/bold]", style="green")
        console.print(f"Number of BAM files: [bold]{len(args.bam)}[/bold]", style="green")
        console.print(f"Number of regions: [bold]{len(regions)}[/bold]", style="green")
        console.print(f"Total positions to analyze: [bold]{get_total_regions_size(regions) * len(args.bam):,}[/bold]", style="green")
    
    # Process regions
    results_df = process_regions(
        args.bam,
        args.reference,
        regions,
        args.min_mapping_quality,
        args.min_base_quality,
        args.min_depth,
        args.verbose,
        not args.no_progress_bar
    )
    
    if results_df.empty:
        console.print("No positions found with the specified criteria.", style="bold red")
        sys.exit(1)

    # Process windows if requested
    if args.window:
        if args.verbose:
            console.print(f"Grouping results into windows of size [bold]{args.window}[/bold] bp", style="blue")
        
        # Create windowed summary data (this replaces the position-level data)
        window_df = process_windows(results_df, args.window)
        
        if window_df.empty:
            console.print("No windows found with the specified criteria.", style="bold red")
            sys.exit(1)
        
        # Define window-specific columns
        if args.detailed:
            columns = [
                'Bam_file', 'CHROM', 'Window_start', 'Window_end', 'Positions_count',
                'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads',
                'Mean_A_count', 'Mean_C_count', 'Mean_G_count', 'Mean_T_count', 
                'Mean_DEL_count', 'Mean_INS_count', 'Mean_N_count'
            ]
        else:
            columns = [
                'Bam_file', 'CHROM', 'Window_start', 'Window_end', 'Positions_count',
                'Mean_depth', 'Mean_ref_reads', 'Mean_alt_reads'
            ]
            
        results_df = window_df
    else:
        # Position-level output (original behavior)
        if args.detailed:
            columns = [
                'Bam_file', 'CHROM', 'POS', 'REF', 'DP', 'Ref_reads', 'Alt_reads',
                'A_count', 'C_count', 'G_count', 'T_count', 'DEL_count', 'INS_count', 'N_count'
            ]
        else:
            columns = ['Bam_file', 'CHROM', 'POS', 'REF', 'DP', 'Ref_reads', 'Alt_reads']

    # Write output TSV
    results_df[columns].to_csv(args.output, sep='\t', index=False)

    if args.verbose:
        console.print(f"[bold green]Results written to {args.output}[/bold green]")
        if args.window:
            console.print(f"Total windows processed: [bold]{len(results_df)}[/bold]", style="blue")
        else:
            console.print(f"Total positions processed: [bold]{len(results_df)}[/bold]", style="blue")
    
    # Generate visualization if requested
    if args.plot:
        if args.verbose:
            console.print("Generating visualization...", style="magenta")
        
        # Parse GTF file if provided
        gtf_features = None
        if args.gtf:
            if args.verbose:
                console.print(f"Parsing gene annotations from [bold]{args.gtf}[/bold]", style="magenta")
            gtf_features = parse_gtf_file(args.gtf)
        
        # Parse figure size
        try:
            figsize = tuple(float(x) for x in args.figsize.split(','))
            if len(figsize) != 2:
                figsize = (12, 8)
        except:
            figsize = (12, 8)
            
        # Plot with window data and gene annotations
        if gtf_features:
            plot_allele_counts_with_annotations(
                results_df,
                gtf_features=gtf_features,
                output_file=args.plot_output,
                dpi=args.dpi,
                figsize=figsize,
                show_plot=not args.no_show_plot,
                max_genes=args.max_genes
            )
        else:
            # Fall back to original plotting function if no GTF
            plot_allele_counts(
                results_df,
                output_file=args.plot_output,
                dpi=args.dpi,
                figsize=figsize,
                show_plot=not args.no_show_plot
            )
        
        if args.verbose and args.plot_output:
            console.print(f"[bold green]Visualization saved to {args.plot_output}[/bold green]")
    
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        console.print("\n[bold red]Process interrupted by user[/bold red]")
        sys.exit(1)
    except Exception as e:
        console.print(f"[bold red]Error: {str(e)}[/bold red]")
        if os.environ.get('ATARI_DEBUG'):
            import traceback
            traceback.print_exc()
        sys.exit(1)