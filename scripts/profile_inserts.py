#!/usr/bin/env python3

import os
import sys
import gzip
import glob
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from fuzzywuzzy import fuzz
from tqdm import tqdm
import re
from pathlib import Path

# Set up plotting style
plt.style.use('default')
sns.set_palette("husl")

def GC(sequence):
    """Calculate GC content percentage."""
    return gc_fraction(sequence) * 100

def setup_logging():
    """Set up logging to both file and console."""
    log_dir = "results/profile_inserts"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "profile_inserts.log")
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger("profile_inserts")

def load_probes(csv_path):
    """Load upstream and downstream probes from CSV."""
    logger = logging.getLogger("profile_inserts")
    try:
        df = pd.read_csv(csv_path)
        if 'upstream' not in df.columns or 'downstream' not in df.columns:
            raise ValueError("CSV must contain 'upstream' and 'downstream' columns")
        
        probes = []
        for idx, row in df.iterrows():
            probes.append({
                'upstream': row['upstream'].upper(),
                'downstream': row['downstream'].upper(),
                'name': f"probe_{idx+1}"
            })
        
        logger.info(f"Loaded {len(probes)} probe pairs from {csv_path}")
        return probes
    except Exception as e:
        logger.error(f"Error loading probes from {csv_path}: {e}")
        raise

def find_probe_positions(sequence, probe, min_ratio=80):
    """Find probe positions in sequence using fuzzy matching."""
    positions = []
    probe_len = len(probe)
    
    # Check forward strand
    for i in range(len(sequence) - probe_len + 1):
        window = sequence[i:i + probe_len]
        ratio = fuzz.ratio(window, probe)
        if ratio >= min_ratio:
            positions.append((i, i + probe_len, ratio, 'forward'))
    
    # Check reverse complement
    rc_probe = str(Seq(probe).reverse_complement())
    for i in range(len(sequence) - probe_len + 1):
        window = sequence[i:i + probe_len]
        ratio = fuzz.ratio(window, rc_probe)
        if ratio >= min_ratio:
            positions.append((i, i + probe_len, ratio, 'reverse'))
    
    return positions

def extract_inserts(fastq_path, probes, min_ratio=80):
    """Extract inserts from FASTQ file using probe pairs."""
    logger = logging.getLogger("profile_inserts")
    inserts = []
    total_reads = 0
    matched_reads = 0
    
    # Count total reads first
    with gzip.open(fastq_path, 'rt') as f:
        total_reads = sum(1 for _ in f) // 4
    
    logger.info(f"Processing {total_reads} reads from {fastq_path}")
    
    with gzip.open(fastq_path, 'rt') as f:
        with tqdm(total=total_reads, desc="Processing reads") as pbar:
            while True:
                header = f.readline()
                if not header:
                    break
                
                seq = f.readline().strip()
                f.readline()  # '+'
                f.readline()  # quality
                
                total_reads += 1
                pbar.update(1)
                
                # Find inserts for each probe pair
                for probe in probes:
                    upstream_positions = find_probe_positions(seq, probe['upstream'], min_ratio)
                    downstream_positions = find_probe_positions(seq, probe['downstream'], min_ratio)
                    
                    total_combinations = len(upstream_positions) * len(downstream_positions)
                    proper_orientation_count = 0
                    
                    for up_start, up_end, up_ratio, up_strand in upstream_positions:
                        for down_start, down_end, down_ratio, down_strand in downstream_positions:
                            # Only extract if both probes are in forward strand (5'→3' direction)
                            # and downstream comes after upstream
                            if (up_strand == 'forward' and down_strand == 'forward' and 
                                down_start > up_end):
                                insert_seq = seq[up_end:down_start]
                                if len(insert_seq) > 0:
                                    inserts.append({
                                        'sequence': insert_seq,
                                        'length': len(insert_seq),
                                        'probe_name': probe['name'],
                                        'up_ratio': up_ratio,
                                        'down_ratio': down_ratio,
                                        'up_strand': up_strand,
                                        'down_strand': down_strand,
                                        'read_id': header.strip()[1:]
                                    })
                                    matched_reads += 1
                                    proper_orientation_count += 1
                    
                    if total_combinations > 0:
                        logger.debug(f"Probe pair {probe['name']}: {proper_orientation_count}/{total_combinations} combinations had proper orientation")
    
    logger.info(f"Extracted {len(inserts)} inserts from {matched_reads} matched reads")
    return inserts

def calculate_qc_metrics(inserts):
    """Calculate comprehensive QC metrics."""
    logger = logging.getLogger("profile_inserts")
    
    if not inserts:
        logger.warning("No inserts found for QC analysis")
        return {}
    
    # Basic statistics
    lengths = [insert['length'] for insert in inserts]
    sequences = [insert['sequence'] for insert in inserts]
    
    metrics = {
        'total_inserts': len(inserts),
        'length_stats': {
            'mean': np.mean(lengths),
            'median': np.median(lengths),
            'std': np.std(lengths),
            'min': min(lengths),
            'max': max(lengths),
            'q25': np.percentile(lengths, 25),
            'q75': np.percentile(lengths, 75)
        },
        'gc_content': np.mean([GC(seq) for seq in sequences]),
        'length_distribution': Counter(lengths),
        'probe_matches': Counter([insert['probe_name'] for insert in inserts]),
        'strand_combinations': Counter([
            f"{insert['up_strand']}-{insert['down_strand']}" 
            for insert in inserts
        ]),
        'match_quality': {
            'up_ratio_mean': np.mean([insert['up_ratio'] for insert in inserts]),
            'down_ratio_mean': np.mean([insert['down_ratio'] for insert in inserts])
        }
    }
    
    # Sequence composition
    all_bases = ''.join(sequences)
    base_counts = Counter(all_bases)
    total_bases = len(all_bases)
    metrics['base_composition'] = {
        base: count/total_bases for base, count in base_counts.items()
    }
    
    # Duplicate analysis
    seq_counts = Counter(sequences)
    metrics['unique_sequences'] = len(seq_counts)
    metrics['duplicate_rate'] = 1 - (len(seq_counts) / len(sequences))
    
    logger.info(f"QC metrics calculated for {len(inserts)} inserts")
    return metrics

def create_qc_plots(inserts, metrics, output_dir):
    """Create comprehensive QC plots."""
    logger = logging.getLogger("profile_inserts")
    
    # Set up the plotting grid
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Length distribution histogram
    plt.subplot(3, 4, 1)
    lengths = [insert['length'] for insert in inserts]
    plt.hist(lengths, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Insert Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Insert Length Distribution')
    plt.axvline(metrics['length_stats']['mean'], color='red', linestyle='--', 
                label=f"Mean: {metrics['length_stats']['mean']:.1f}")
    plt.axvline(metrics['length_stats']['median'], color='green', linestyle='--', 
                label=f"Median: {metrics['length_stats']['median']:.1f}")
    plt.legend()
    
    # 2. GC content distribution
    plt.subplot(3, 4, 2)
    gc_contents = [GC(insert['sequence']) for insert in inserts]
    plt.hist(gc_contents, bins=30, alpha=0.7, edgecolor='black')
    plt.xlabel('GC Content (%)')
    plt.ylabel('Frequency')
    plt.title('GC Content Distribution')
    plt.axvline(metrics['gc_content'], color='red', linestyle='--', 
                label=f"Mean: {metrics['gc_content']:.1f}%")
    plt.legend()
    
    # 3. Probe match quality
    plt.subplot(3, 4, 3)
    up_ratios = [insert['up_ratio'] for insert in inserts]
    down_ratios = [insert['down_ratio'] for insert in inserts]
    plt.scatter(up_ratios, down_ratios, alpha=0.6)
    plt.xlabel('Upstream Probe Match Ratio (%)')
    plt.ylabel('Downstream Probe Match Ratio (%)')
    plt.title('Probe Match Quality')
    
    # 4. Strand combination analysis
    plt.subplot(3, 4, 4)
    strand_combos = [f"{insert['up_strand']}-{insert['down_strand']}" for insert in inserts]
    combo_counts = Counter(strand_combos)
    plt.bar(combo_counts.keys(), combo_counts.values())
    plt.xlabel('Strand Combination')
    plt.ylabel('Count')
    plt.title('Strand Combination Analysis')
    plt.xticks(rotation=45)
    
    # 5. Base composition
    plt.subplot(3, 4, 5)
    bases = list(metrics['base_composition'].keys())
    compositions = list(metrics['base_composition'].values())
    plt.bar(bases, compositions)
    plt.xlabel('Base')
    plt.ylabel('Frequency')
    plt.title('Base Composition')
    
    # 6. Probe performance
    plt.subplot(3, 4, 6)
    probe_counts = metrics['probe_matches']
    plt.bar(probe_counts.keys(), probe_counts.values())
    plt.xlabel('Probe')
    plt.ylabel('Insert Count')
    plt.title('Probe Performance')
    plt.xticks(rotation=45)
    
    # 7. Length vs GC content
    plt.subplot(3, 4, 7)
    plt.scatter(lengths, gc_contents, alpha=0.6)
    plt.xlabel('Insert Length (bp)')
    plt.ylabel('GC Content (%)')
    plt.title('Length vs GC Content')
    
    # 8. Cumulative length distribution
    plt.subplot(3, 4, 8)
    sorted_lengths = sorted(lengths)
    cumulative = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
    plt.plot(sorted_lengths, cumulative)
    plt.xlabel('Insert Length (bp)')
    plt.ylabel('Cumulative Fraction')
    plt.title('Cumulative Length Distribution')
    
    # 9. Quality score distribution
    plt.subplot(3, 4, 9)
    quality_scores = [(insert['up_ratio'] + insert['down_ratio'])/2 for insert in inserts]
    plt.hist(quality_scores, bins=30, alpha=0.7, edgecolor='black')
    plt.xlabel('Average Match Quality (%)')
    plt.ylabel('Frequency')
    plt.title('Match Quality Distribution')
    
    # 10. Length statistics box plot
    plt.subplot(3, 4, 10)
    plt.boxplot(lengths)
    plt.ylabel('Insert Length (bp)')
    plt.title('Length Statistics')
    
    # 11. Duplicate analysis
    plt.subplot(3, 4, 11)
    seq_counts = Counter([insert['sequence'] for insert in inserts])
    duplicate_counts = [count for count in seq_counts.values() if count > 1]
    if duplicate_counts:
        plt.hist(duplicate_counts, bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Duplicate Count')
        plt.ylabel('Frequency')
        plt.title('Duplicate Distribution')
    else:
        plt.text(0.5, 0.5, 'No duplicates found', ha='center', va='center', 
                transform=plt.gca().transAxes)
        plt.title('Duplicate Distribution')
    
    # 12. Summary statistics
    plt.subplot(3, 4, 12)
    plt.axis('off')
    summary_text = f"""
    Summary Statistics:
    
    Total Inserts: {metrics['total_inserts']}
    Mean Length: {metrics['length_stats']['mean']:.1f} bp
    Median Length: {metrics['length_stats']['median']:.1f} bp
    GC Content: {metrics['gc_content']:.1f}%
    Unique Sequences: {metrics['unique_sequences']}
    Duplicate Rate: {metrics['duplicate_rate']:.2%}
    Mean Match Quality: {metrics['match_quality']['up_ratio_mean']:.1f}%
    """
    plt.text(0.1, 0.9, summary_text, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    plot_path = os.path.join(output_dir, "qc_plots.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"QC plots saved to {plot_path}")

def save_inserts_fasta(inserts, output_dir):
    """Save inserts to FASTA file."""
    fasta_path = os.path.join(output_dir, "extracted_inserts.fasta")
    
    with open(fasta_path, 'w') as f:
        for i, insert in enumerate(inserts):
            header = f">{insert['probe_name']}_insert_{i+1}_len{insert['length']}_up{insert['up_ratio']:.0f}_down{insert['down_ratio']:.0f}"
            f.write(f"{header}\n{insert['sequence']}\n")
    
    logging.getLogger("profile_inserts").info(f"Saved {len(inserts)} inserts to {fasta_path}")

def save_qc_report(metrics, output_dir):
    """Save detailed QC report."""
    report_path = os.path.join(output_dir, "qc_report.txt")
    
    with open(report_path, 'w') as f:
        f.write("PROFILE INSERTS QC REPORT\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("SUMMARY STATISTICS:\n")
        f.write(f"Total inserts: {metrics['total_inserts']}\n")
        f.write(f"Mean length: {metrics['length_stats']['mean']:.2f} bp\n")
        f.write(f"Median length: {metrics['length_stats']['median']:.2f} bp\n")
        f.write(f"Standard deviation: {metrics['length_stats']['std']:.2f} bp\n")
        f.write(f"Min length: {metrics['length_stats']['min']} bp\n")
        f.write(f"Max length: {metrics['length_stats']['max']} bp\n")
        f.write(f"GC content: {metrics['gc_content']:.2f}%\n")
        f.write(f"Unique sequences: {metrics['unique_sequences']}\n")
        f.write(f"Duplicate rate: {metrics['duplicate_rate']:.2%}\n\n")
        
        f.write("LENGTH STATISTICS:\n")
        f.write(f"Q25: {metrics['length_stats']['q25']:.2f} bp\n")
        f.write(f"Q75: {metrics['length_stats']['q75']:.2f} bp\n\n")
        
        f.write("PROBE PERFORMANCE:\n")
        for probe, count in metrics['probe_matches'].items():
            f.write(f"{probe}: {count} inserts\n")
        f.write("\n")
        
        f.write("STRAND COMBINATIONS:\n")
        for combo, count in metrics['strand_combinations'].items():
            f.write(f"{combo}: {count} inserts\n")
        f.write("\n")
        
        f.write("BASE COMPOSITION:\n")
        for base, freq in metrics['base_composition'].items():
            f.write(f"{base}: {freq:.3f}\n")
        f.write("\n")
        
        f.write("MATCH QUALITY:\n")
        f.write(f"Mean upstream ratio: {metrics['match_quality']['up_ratio_mean']:.2f}%\n")
        f.write(f"Mean downstream ratio: {metrics['match_quality']['down_ratio_mean']:.2f}%\n")
    
    logging.getLogger("profile_inserts").info(f"QC report saved to {report_path}")

def main():
    """Main function to run the profile inserts analysis."""
    logger = setup_logging()
    logger.info("Starting profile inserts analysis")
    
    try:
        # Check input files
        data_dir = "data/profileinserts"
        if not os.path.exists(data_dir):
            raise FileNotFoundError(f"Data directory {data_dir} not found")
        
        csv_files = glob.glob(os.path.join(data_dir, "*.csv"))
        fastq_files = glob.glob(os.path.join(data_dir, "*.fastq.gz"))
        
        if not csv_files:
            raise FileNotFoundError(f"No CSV files found in {data_dir}")
        if not fastq_files:
            raise FileNotFoundError(f"No FASTQ files found in {data_dir}")
        
        csv_path = csv_files[0]
        fastq_path = fastq_files[0]
        
        logger.info(f"Using CSV: {csv_path}")
        logger.info(f"Using FASTQ: {fastq_path}")
        
        # Load probes
        probes = load_probes(csv_path)
        
        # Extract inserts
        inserts = extract_inserts(fastq_path, probes)
        
        if not inserts:
            logger.warning("No inserts found!")
            return
        
        # Calculate QC metrics
        metrics = calculate_qc_metrics(inserts)
        
        # Create output directory
        output_dir = "results/profile_inserts"
        os.makedirs(output_dir, exist_ok=True)
        
        # Save results
        save_inserts_fasta(inserts, output_dir)
        save_qc_report(metrics, output_dir)
        create_qc_plots(inserts, metrics, output_dir)
        
        logger.info("Profile inserts analysis completed successfully")
        print(f"Analysis complete! Results saved to {output_dir}")
        
    except Exception as e:
        logger.exception(f"Error in profile inserts analysis: {e}")
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 