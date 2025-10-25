#!/usr/bin/env python3

import subprocess
import sys
import glob
import os
import shutil
import logging
import pysam
import random
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import math
import tempfile
from pathlib import Path

# Use a built-in Matplotlib style ("ggplot") for consistency
plt.style.use("ggplot")

# Try to import scipy for Z‐test p‐value and KDE
try:
    from scipy.stats import norm, gaussian_kde, beta, binom
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
    print("Warning: scipy is not installed. Z-test p-values and KDE may be limited.")

# Base directories
DATA_DIR         = os.path.join("data", "ep-library-profile")
BASE_WORK_DIR    = os.path.join(DATA_DIR, "tmp")
BASE_RESULTS_DIR = os.path.join("results", "ep-library-profile")

def setup_logging(log_dir):
    """
    Configure logging so that INFO (and above) go into run.log inside log_dir.
    Before configuring, remove any existing handlers to avoid carryover.
    """
    # Remove any existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "run.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            # INFO logs go into file only
        ]
    )

def run_minimap2(reads_input, ref_fasta, out_prefix, work_dir):
    """
    Runs minimap2 on reads_input (either one FASTQ or a glob pattern).
    Saves SAM to work_dir/out_prefix.sam.
    Returns the full path to the SAM file.
    """
    # Expand glob if pattern contains wildcard; otherwise assume it's a single path
    fastq_files = glob.glob(reads_input) if ("*" in reads_input or "?" in reads_input) else [reads_input]
    if not fastq_files:
        logging.error(f"No FASTQ files found matching {reads_input}")
        raise FileNotFoundError(f"No FASTQ files found matching {reads_input}")

    sam_output = os.path.join(work_dir, f"{out_prefix}.sam")
    cmd = ["minimap2", "-ax", "map-ont", ref_fasta, *fastq_files]
    logging.info(f"Running minimap2: {' '.join(cmd)} → {sam_output}")
    with open(sam_output, "w") as out_sam:
        subprocess.run(cmd, stdout=out_sam, check=True)
    return sam_output

def load_single_sequence(fasta_path):
    """
    Loads exactly one sequence from a FASTA. Raises if not exactly one record.
    Returns (sequence_string, sequence_id).
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) != 1:
        logging.error(f"Expected exactly 1 record in {fasta_path}, found {len(records)}.")
        raise ValueError(f"Expected exactly 1 record in {fasta_path}, found {len(records)}.")
    seq_id = records[0].id
    seq_str = str(records[0].seq)
    logging.info(f"Loaded sequence {seq_id} from {fasta_path}")
    return seq_str, seq_id

def create_multi_fasta(chunks, work_dir, out_fasta_prefix="plasmid_chunks"):
    """
    Writes each chunk as a separate FASTA entry into work_dir/out_fasta_prefix.fasta:
      >chunk_1
      ACTG...
      >chunk_2
      TTAG...
    Returns the path to the created FASTA.
    """
    out_fasta = os.path.join(work_dir, f"{out_fasta_prefix}.fasta")
    logging.info(f"Writing {len(chunks)} chunks to {out_fasta}")
    with open(out_fasta, "w") as f:
        for i, seq in enumerate(chunks):
            f.write(f">chunk_{i+1}\n")
            f.write(str(seq) + "\n")
    return out_fasta

def calculate_background_from_plasmid(sam_plasmid, plasmid_seq, target_start, target_length):
    """
    Calculate background mismatch statistics from full plasmid alignment, excluding target region.
    
    Args:
        sam_plasmid: Path to SAM file with full plasmid alignment
        plasmid_seq: Full plasmid sequence
        target_start: Start position of target region (0-based)
        target_length: Length of target region
    
    Returns:
        tuple: (total_mismatches, total_covered_bases, mapped_reads)
    """
    target_end = target_start + target_length
    
    # Initialize counters
    total_mismatches = 0
    total_covered_bases = 0
    mapped_reads = 0
    
    samfile = pysam.AlignmentFile(sam_plasmid, "r")
    for read in samfile.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        
        mapped_reads += 1
        
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if read_pos is not None and ref_pos is not None and 0 <= ref_pos < len(plasmid_seq):
                # Skip positions within the target region
                if target_start <= ref_pos < target_end:
                    continue
                
                total_covered_bases += 1
                if read.query_sequence[read_pos].upper() != plasmid_seq[ref_pos].upper():
                    total_mismatches += 1
    
    samfile.close()
    
    return total_mismatches, total_covered_bases, mapped_reads

def calculate_rolling_mutation_rate_plasmid(sam_plasmid, plasmid_seq, window_size=20):
    """
    Calculate rolling mutation rate across the entire plasmid with a specified window size.
    
    Args:
        sam_plasmid: Path to SAM file with full plasmid alignment
        plasmid_seq: Full plasmid sequence
        window_size: Size of the rolling window (default: 20 bp)
    
    Returns:
        tuple: (positions, rolling_rates) for the entire plasmid
    """
    # Calculate per-position mutation rates for the entire plasmid
    plasmid_mismatches = [0] * len(plasmid_seq)
    plasmid_coverage = [0] * len(plasmid_seq)
    
    samfile = pysam.AlignmentFile(sam_plasmid, "r")
    for read in samfile.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if read_pos is not None and ref_pos is not None and 0 <= ref_pos < len(plasmid_seq):
                plasmid_coverage[ref_pos] += 1
                if read.query_sequence[read_pos].upper() != plasmid_seq[ref_pos].upper():
                    plasmid_mismatches[ref_pos] += 1
    
    samfile.close()
    
    # Calculate per-position mutation rates
    plasmid_rates = []
    for cov, mis in zip(plasmid_coverage, plasmid_mismatches):
        if cov > 0:
            plasmid_rates.append(mis / cov)
        else:
            plasmid_rates.append(0.0)
    
    # Calculate rolling average
    rolling_rates = []
    positions = []
    
    for i in range(len(plasmid_rates) - window_size + 1):
        window_rates = plasmid_rates[i:i + window_size]
        rolling_rates.append(np.mean(window_rates))
        positions.append(i + window_size // 2 + 1)  # Center position of window
    
    return positions, rolling_rates

def compute_mismatch_stats_sam(sam_file, refs_dict):
    """
    For each reference in refs_dict:
      - Count how many reads mapped to it
      - Count per‐position coverage
      - Count per‐position mismatches
      - Compute total mismatches, total covered bases, average mismatch rate

    Returns a dict keyed by reference name, each containing:
        pos_rates             -> list of mismatch rates per position
        cov                   -> list of coverage per position
        mismatch              -> list of mismatch counts per position
        avg_mismatch_rate     -> float (raw per‐base fraction)
        total_mismatches      -> int
        total_covered_bases   -> int
        mapped_reads          -> int
    """
    mismatch_data = {
        name: {
            "pos_rates": [0.0] * len(seq),
            "cov": [0] * len(seq),
            "mismatch": [0] * len(seq),
            "avg_mismatch_rate": 0.0,
            "total_mismatches": 0,
            "total_covered_bases": 0,
            "mapped_reads": 0,
        }
        for name, seq in refs_dict.items()
    }

    logging.info(f"Computing mismatch stats for {sam_file}")
    samfile = pysam.AlignmentFile(sam_file, "r")
    for read in samfile.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        ref_name = samfile.get_reference_name(read.reference_id)
        if ref_name not in mismatch_data:
            continue

        mismatch_data[ref_name]["mapped_reads"] += 1
        info = mismatch_data[ref_name]
        ref_seq = refs_dict[ref_name]

        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if read_pos is not None and ref_pos is not None and 0 <= ref_pos < len(ref_seq):
                info["cov"][ref_pos] += 1
                if read.query_sequence[read_pos].upper() != ref_seq[ref_pos].upper():
                    info["mismatch"][ref_pos] += 1

    samfile.close()

    for name, info in mismatch_data.items():
        total_mis = 0
        total_cov = 0
        pos_rates = []
        for cov, mis in zip(info["cov"], info["mismatch"]):
            if cov > 0:
                rate = mis / cov
                pos_rates.append(rate)
                total_mis += mis
                total_cov += cov
            else:
                pos_rates.append(0.0)

        info["pos_rates"] = pos_rates
        info["total_mismatches"] = total_mis
        info["total_covered_bases"] = total_cov
        info["avg_mismatch_rate"] = (total_mis / total_cov) if total_cov > 0 else 0.0

        logging.info(
            f"{name}: mapped_reads={info['mapped_reads']}, "
            f"mismatches={total_mis}, covered_bases={total_cov}, "
            f"avg_rate={info['avg_mismatch_rate']:.6f}"
        )

    return mismatch_data

def z_test_two_proportions(mis1, cov1, mis2, cov2):
    """
    Performs a two‐proportion Z‐test:
      H0: p1 == p2
      H1: p1 != p2

    Returns (z_statistic, p_value). If scipy is unavailable, p_value=None.
    """
    if cov1 == 0 or cov2 == 0:
        return 0.0, 1.0

    p1 = mis1 / cov1
    p2 = mis2 / cov2
    p  = (mis1 + mis2) / (cov1 + cov2)
    denom = math.sqrt(p * (1 - p) * (1 / cov1 + 1 / cov2))
    if denom == 0:
        return 0.0, 1.0

    z_stat = (p1 - p2) / denom
    if HAVE_SCIPY:
        p_val = 2 * (1 - norm.cdf(abs(z_stat)))
        logging.info(f"Z‐test: z={z_stat:.4f}, p‐value={p_val:.4e}")
    else:
        p_val = None
        logging.info(f"Z‐test: z={z_stat:.4f}, p‐value=(scipy unavailable)")

    return z_stat, p_val

def run_nanofilt_filtering(input_fastq, quality_threshold, output_fastq):
    """
    Run NanoFilt to filter FASTQ file by quality score threshold and minimum length.
    
    Args:
        input_fastq: Path to input FASTQ.gz file
        quality_threshold: Quality score threshold (integer)
        output_fastq: Path to output filtered FASTQ.gz file
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Use gunzip to decompress, pipe to NanoFilt with length filter, then compress output
        cmd = f"gunzip -c {input_fastq} | NanoFilt -q {quality_threshold} -l 30 | gzip > {output_fastq}"
        logging.info(f"Running NanoFilt with quality threshold {quality_threshold} and min length 30bp: {cmd}")
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"NanoFilt failed with return code {result.returncode}: {result.stderr}")
            return False
        
        # Check if output file was created and has content
        if os.path.exists(output_fastq) and os.path.getsize(output_fastq) > 0:
            logging.info(f"Successfully created filtered FASTQ: {output_fastq}")
            return True
        else:
            logging.error(f"Output file {output_fastq} was not created or is empty")
            return False
            
    except Exception as e:
        logging.error(f"Error running NanoFilt: {e}")
        return False

def calculate_mutation_rate_for_quality(fastq_path, quality_threshold, work_dir, ref_hit_fasta, plasmid_fasta):
    """
    Calculate comprehensive AA mutation analysis for a specific quality threshold.
    
    Args:
        fastq_path: Path to input FASTQ.gz file
        quality_threshold: Quality score threshold
        work_dir: Working directory for temporary files
        ref_hit_fasta: Path to reference hit FASTA
        plasmid_fasta: Path to plasmid FASTA
    
    Returns:
        dict: Comprehensive results with error estimates or None if failed
    """
    try:
        # Create filtered FASTQ file
        filtered_fastq = os.path.join(work_dir, f"filtered_q{quality_threshold}.fastq.gz")
        
        if not run_nanofilt_filtering(fastq_path, quality_threshold, filtered_fastq):
            return None
        
        # Load sequences
        hit_seq, hit_id = load_single_sequence(ref_hit_fasta)
        plasmid_seq, plasmid_id = load_single_sequence(plasmid_fasta)
        
        # Find hit region in plasmid
        idx = plasmid_seq.upper().find(hit_seq.upper())
        if idx == -1:
            logging.error("Gene region not found in plasmid")
            return None
        
        # Align filtered reads to hit region
        sam_hit = run_minimap2(filtered_fastq, ref_hit_fasta, f"hit_q{quality_threshold}", work_dir)
        
        # Align filtered reads to full plasmid for background calculation
        sam_plasmid = run_minimap2(filtered_fastq, plasmid_fasta, f"plasmid_q{quality_threshold}", work_dir)
        
        # Calculate background rate from full plasmid alignment, excluding target region
        bg_mis, bg_cov, bg_reads = calculate_background_from_plasmid(sam_plasmid, plasmid_seq, idx, len(hit_seq))
        
        # Calculate hit region mutation rate
        mismatch_hit = compute_mismatch_stats_sam(sam_hit, {hit_id: hit_seq})
        hit_info = mismatch_hit[hit_id]
        hit_mis = hit_info["total_mismatches"]
        hit_cov = hit_info["total_covered_bases"]
        
        # Extract Q-score statistics for both hit and background regions
        hit_qscore_stats = extract_qscores_from_sam(sam_hit)
        bg_qscore_stats = extract_qscores_from_sam(sam_plasmid)
        
        # Check if it's a protein-coding sequence
        is_protein = True
        seq_upper = hit_seq.upper()
        if len(seq_upper) % 3 != 0:
            is_protein = False
        elif "*" in str(Seq(seq_upper).translate(to_stop=False))[:-1]:
            is_protein = False
        
        # Run comprehensive analysis if protein-coding
        if is_protein:
            results = comprehensive_aa_mutation_analysis(
                hit_mis, hit_cov, bg_mis, bg_cov, hit_seq, 
                quality_threshold=quality_threshold, n_trials=10000,
                hit_qscore_stats=hit_qscore_stats, bg_qscore_stats=bg_qscore_stats,
                sam_hit=sam_hit, sam_plasmid=sam_plasmid, hit_seq=hit_seq, plasmid_seq=plasmid_seq
            )
            return results
        else:
            # For non-protein sequences, return basic info
            return {
                'mean_aa_mutations': 0.0,
                'std_aa_mutations': 0.0,
                'ci_lower': 0.0,
                'ci_upper': 0.0,
                'hit_rate': hit_mis / hit_cov if hit_cov > 0 else 0,
                'bg_rate': bg_mis / bg_cov if bg_cov > 0 else 0,
                'net_rate': max((hit_mis / hit_cov) - (bg_mis / bg_cov), 0) if hit_cov > 0 and bg_cov > 0 else 0,
                'mappable_bases': hit_cov,
                'quality_threshold': quality_threshold,
                'is_protein': False
            }
        
    except Exception as e:
        logging.error(f"Error calculating mutation rate for quality {quality_threshold}: {e}")
        return None

def find_optimal_qscore(qc_results):
    """
    Find the Q-score threshold with the lowest net mutation rate error.
    
    Args:
        qc_results: List of comprehensive analysis results
    
    Returns:
        tuple: (optimal_qscore, optimal_result, error_comparison)
    """
    logging.info("=== FINDING OPTIMAL Q-SCORE THRESHOLD ===")
    
    if not qc_results:
        return None, None, None
    
    # Find Q-score with minimum net mutation rate error
    min_error = float('inf')
    optimal_result = None
    optimal_qscore = None
    
    error_comparison = []
    
    for result in qc_results:
        qscore = result['quality_threshold']
        # Use weighted error for optimal Q-score selection
        net_rate_error = result.get('net_weighted_error', result['net_rate_error'])
        mappable_bases = result['mappable_bases']
        
        error_comparison.append({
            'qscore': qscore,
            'net_rate_error': net_rate_error,
            'net_weighted_error': result.get('net_weighted_error', 0.0),
            'mappable_bases': mappable_bases,
            'aa_mutations': result['mean_aa_mutations'],
            'aa_error': result['std_aa_mutations']
        })
        
        logging.info(f"Q{qscore}: net_weighted_error={net_rate_error:.6f}, mappable_bases={mappable_bases}")
        
        if net_rate_error < min_error:
            min_error = net_rate_error
            optimal_result = result
            optimal_qscore = qscore
    
    logging.info(f"OPTIMAL Q-SCORE: Q{optimal_qscore} (lowest net mutation rate error: {min_error:.6f})")
    logging.info(f"Optimal result: AA mutations = {optimal_result['mean_aa_mutations']:.4f} ± {optimal_result['std_aa_mutations']:.4f}")
    
    return optimal_qscore, optimal_result, error_comparison

def run_qc_analysis(fastq_path, results_dir):
    """
    Run simple QC analysis using segmentation-based error estimation.
    
    Args:
        fastq_path: Path to input FASTQ.gz file
        results_dir: Directory to save QC results
    """
    logging.info("Starting simple QC analysis with segmentation-based error estimation")
    
    # Define quality thresholds to test
    quality_thresholds = [10, 12, 14, 16, 18, 20, 22, 24, 26]
    
    # Paths to reference files
    ref_hit_fasta = os.path.join(DATA_DIR, "region_of_interest.fasta")
    plasmid_fasta = os.path.join(DATA_DIR, "plasmid.fasta")
    
    # Segment the input FASTQ file
    logging.info("Segmenting FASTQ file into 10 parts for error estimation...")
    segment_files = segment_fastq_file(fastq_path, n_segments=10)
    
    if not segment_files:
        logging.error("Failed to segment FASTQ file")
        return
    
    # Create temporary work directory for QC analysis
    with tempfile.TemporaryDirectory() as qc_work_dir:
        logging.info(f"Using temporary work directory: {qc_work_dir}")
        
        # Calculate results for each quality threshold
        qc_results = []
        successful_thresholds = []
        
        for q_threshold in quality_thresholds:
            logging.info(f"Processing quality threshold: {q_threshold}")
            result = run_segmented_analysis(
                segment_files, q_threshold, qc_work_dir, ref_hit_fasta, plasmid_fasta
            )
            
            if result is not None:
                qc_results.append(result)
                successful_thresholds.append(q_threshold)
                logging.info(f"Quality {q_threshold}: AA mutations = {result['mean_aa_mutations']:.4f} ± {result['std_aa_mutations']:.4f}, "
                           f"mappable bases = {result['total_mappable_bases']}")
            else:
                logging.warning(f"Failed to calculate mutation rate for quality threshold {q_threshold}")
        
        # Find optimal Q-score threshold (lowest empirical error)
        optimal_qscore, optimal_result = find_optimal_qscore_simple(qc_results)
        
        # Create QC plots
        if len(qc_results) >= 2:
            create_simple_qc_plots(successful_thresholds, qc_results, results_dir, optimal_qscore, optimal_result)
        else:
            logging.warning("Insufficient data points for QC plots (need at least 2)")
        
        # Save optimal Q-score information
        if optimal_qscore is not None:
            optimal_qscore_path = os.path.join(results_dir, "optimal_qscore_analysis.txt")
            with open(optimal_qscore_path, 'w') as f:
                f.write("=== OPTIMAL Q-SCORE ANALYSIS (PRECISION-WEIGHTED) ===\n")
                f.write(f"Optimal Q-score threshold: {optimal_qscore}\n")
                f.write(f"Precision-weighted score: {(1.0 / optimal_result['std_aa_mutations']) * optimal_qscore:.6f}\n" if optimal_result['std_aa_mutations'] > 0 else "Precision-weighted score: inf (perfect precision)\n")
                f.write(f"Empirical error (std): {optimal_result['std_aa_mutations']:.6f}\n")
                f.write(f"AA mutations per gene: {optimal_result['mean_aa_mutations']:.4f} ± {optimal_result['std_aa_mutations']:.4f}\n")
                f.write(f"95% Confidence Interval: [{optimal_result['ci_lower']:.4f}, {optimal_result['ci_upper']:.4f}]\n")
                f.write(f"Total mappable bases: {optimal_result['total_mappable_bases']}\n")
                f.write(f"Number of segments: {optimal_result['n_segments']}\n")
                f.write("\n=== ALL Q-SCORE COMPARISON ===\n")
                f.write("Q-score\tEmpirical_Error\tPrecision_Score\tMappable_Bases\tAA_Mutations\tCI_Lower\tCI_Upper\n")
                for result in qc_results:
                    precision_score = (1.0 / result['std_aa_mutations']) * result['quality_threshold'] if result['std_aa_mutations'] > 0 else float('inf')
                    f.write(f"{result['quality_threshold']}\t{result['std_aa_mutations']:.6f}\t{precision_score:.6f}\t{result['total_mappable_bases']}\t{result['mean_aa_mutations']:.4f}\t{result['ci_lower']:.4f}\t{result['ci_upper']:.4f}\n")
            
            logging.info(f"Optimal Q-score analysis saved to: {optimal_qscore_path}")
        
        # Clean up segment files
        import shutil
        segment_dir = os.path.dirname(segment_files[0])
        if os.path.exists(segment_dir):
            shutil.rmtree(segment_dir)
            logging.info(f"Cleaned up segment directory: {segment_dir}")
        
        # Return both QC results and optimal Q-score for use in main analysis
        return qc_results, optimal_qscore

def find_optimal_qscore_simple(qc_results):
    """
    Find the Q-score threshold with the highest precision-weighted score.
    Precision-weighted score = (1 / standard_deviation) * q_score
    
    Args:
        qc_results: List of segmentation analysis results
    
    Returns:
        tuple: (optimal_qscore, optimal_result)
    """
    logging.info("=== FINDING OPTIMAL Q-SCORE THRESHOLD (PRECISION-WEIGHTED) ===")
    
    if not qc_results:
        return None, None
    
    # Find Q-score with highest precision-weighted score
    max_score = -1
    optimal_result = None
    optimal_qscore = None
    
    logging.info("Q-score\tEmpirical_Error\tPrecision_Score\tMappable_Bases")
    logging.info("-" * 60)
    
    for result in qc_results:
        qscore = result['quality_threshold']
        empirical_error = result['std_aa_mutations']
        mappable_bases = result['total_mappable_bases']
        
        # Calculate precision-weighted score: (1/sd) * q_score
        if empirical_error > 0:
            precision_score = (1.0 / empirical_error) * qscore
        else:
            precision_score = float('inf')  # Perfect precision
        
        logging.info(f"Q{qscore}\t{empirical_error:.6f}\t{precision_score:.6f}\t{mappable_bases}")
        
        if precision_score > max_score:
            max_score = precision_score
            optimal_result = result
            optimal_qscore = qscore
    
    logging.info("-" * 60)
    logging.info(f"OPTIMAL Q-SCORE: Q{optimal_qscore} (highest precision-weighted score: {max_score:.6f})")
    logging.info(f"Optimal result: AA mutations = {optimal_result['mean_aa_mutations']:.4f} ± {optimal_result['std_aa_mutations']:.4f}")
    
    return optimal_qscore, optimal_result

def create_simple_qc_plots(quality_thresholds, qc_results, results_dir, optimal_qscore=None, optimal_result=None):
    """
    Create simple QC plots with empirical error bars.
    
    Args:
        quality_thresholds: List of quality score thresholds
        qc_results: List of segmentation analysis results
        results_dir: Directory to save the plots
        optimal_qscore: Optimal Q-score threshold (optional)
        optimal_result: Optimal result data (optional)
    """
    try:
        # Extract data for plotting
        aa_mutations = [r['mean_aa_mutations'] for r in qc_results]
        aa_errors = [r['std_aa_mutations'] for r in qc_results]
        aa_ci_lower = [r['ci_lower'] for r in qc_results]
        aa_ci_upper = [r['ci_upper'] for r in qc_results]
        mappable_bases = [r['total_mappable_bases'] for r in qc_results]
        
        # Create main QC plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Top plot: AA mutations per gene with empirical error bars
        color1 = '#2E8B57'
        ax1.errorbar(quality_thresholds, aa_mutations, yerr=aa_errors, 
                    fmt='o', capsize=5, capthick=2, markersize=8, 
                    color=color1, ecolor=color1, alpha=0.8, label='Mean ± Empirical Std')
        
        # Add confidence intervals as shaded area
        ax1.fill_between(quality_thresholds, aa_ci_lower, aa_ci_upper, 
                        alpha=0.3, color=color1, label='95% Confidence Interval')
        
        # Highlight optimal Q-score
        if optimal_qscore is not None:
            ax1.axvline(x=optimal_qscore, color='red', linestyle='--', alpha=0.7, 
                       label=f'Optimal Q{optimal_qscore}')
        
        ax1.set_xlabel('Quality Score Threshold', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Estimated AA Mutations per Gene', fontsize=12, fontweight='bold', color=color1)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1.set_title('AA Mutations per Gene vs Quality Score Filter (Segmentation-Based Error)', 
                     fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(frameon=False, fontsize=10)
        
        # Add data point labels
        for i, (q, aa_mut, aa_err) in enumerate(zip(quality_thresholds, aa_mutations, aa_errors)):
            ax1.annotate(f'Q{q}\n{aa_mut:.3f}±{aa_err:.3f}', 
                        (q, aa_mut), xytext=(5, 5), 
                        textcoords='offset points', fontsize=8, alpha=0.8, color=color1)
        
        # Bottom plot: Mappable bases and AA mutations per gene
        color2 = '#FF6B6B'
        color3 = '#4169E1'
        
        # Mappable bases (left y-axis)
        ax2_twin = ax2.twinx()
        ax2_twin.scatter(quality_thresholds, mappable_bases, 
                        s=100, alpha=0.7, color=color2, edgecolors='black', 
                        linewidth=1, marker='s', label='Mappable Bases')
        ax2_twin.set_ylabel('Number of Mappable Bases', fontsize=12, fontweight='bold', color=color2)
        ax2_twin.tick_params(axis='y', labelcolor=color2)
        
        # AA mutations per gene with error bars (right y-axis)
        ax2.errorbar(quality_thresholds, aa_mutations, yerr=aa_errors,
                    fmt='^', capsize=5, capthick=2, markersize=8,
                    color=color3, ecolor=color3, alpha=0.8, label='AA Mutations ± Empirical Error')
        ax2.set_ylabel('Estimated AA Mutations per Gene', fontsize=12, fontweight='bold', color=color3)
        ax2.tick_params(axis='y', labelcolor=color3)
        ax2.set_xlabel('Quality Score Threshold', fontsize=12, fontweight='bold')
        ax2.set_title('Mappable Bases and AA Mutations per Gene vs Quality Score Filter', 
                     fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        # Add legends
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2_twin.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', frameon=False, fontsize=10)
        
        # Add data point labels for mappable bases
        for i, (q, bases) in enumerate(zip(quality_thresholds, mappable_bases)):
            ax2_twin.annotate(f'{bases}', (q, bases), xytext=(5, -15), 
                             textcoords='offset points', fontsize=8, alpha=0.8, color=color2)
        
        plt.tight_layout()
        
        # Save the plot
        project_name = os.path.basename(results_dir)
        qc_plot_path = os.path.join(results_dir, f"qc_plot_{project_name}.png")
        fig.savefig(qc_plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        logging.info(f"QC plot saved to: {qc_plot_path}")
        
        # Save data as CSV
        qc_data_path = os.path.join(results_dir, "simple_qc_data.csv")
        with open(qc_data_path, 'w') as f:
            f.write("quality_threshold,mean_aa_mutations,std_aa_mutations,ci_lower,ci_upper,")
            f.write("total_mappable_bases,n_segments\n")
            
            for q, r in zip(quality_thresholds, qc_results):
                f.write(f"{q},{r['mean_aa_mutations']:.6f},{r['std_aa_mutations']:.6f},")
                f.write(f"{r['ci_lower']:.6f},{r['ci_upper']:.6f},")
                f.write(f"{r['total_mappable_bases']},{r['n_segments']}\n")
        
        logging.info(f"Simple QC data saved to: {qc_data_path}")
        
    except Exception as e:
        logging.error(f"Error creating simple QC plots: {e}")
    """
    Create comprehensive QC plots with error bars and uncertainty quantification.
    
    Args:
        quality_thresholds: List of quality score thresholds
        qc_results: List of comprehensive analysis results
        results_dir: Directory to save the plots
        optimal_qscore: Optimal Q-score threshold (optional)
        optimal_result: Optimal result data (optional)
    """
    try:
        # Extract data for plotting
        aa_mutations = [r['mean_aa_mutations'] for r in qc_results]
        aa_errors = [r['std_aa_mutations'] for r in qc_results]
        aa_ci_lower = [r['ci_lower'] for r in qc_results]
        aa_ci_upper = [r['ci_upper'] for r in qc_results]
        mappable_bases = [r['mappable_bases'] for r in qc_results]
        net_rates = [r['net_rate'] for r in qc_results]
        net_rate_errors = [r['net_rate_error'] for r in qc_results]
        
        # Create main QC plot with error bars
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
        
        # Top plot: AA mutations per gene with error bars
        color1 = '#2E8B57'
        ax1.errorbar(quality_thresholds, aa_mutations, yerr=aa_errors, 
                    fmt='o', capsize=5, capthick=2, markersize=8, 
                    color=color1, ecolor=color1, alpha=0.8, label='Mean ± Std')
        
        # Add confidence intervals as shaded area
        ax1.fill_between(quality_thresholds, aa_ci_lower, aa_ci_upper, 
                        alpha=0.3, color=color1, label='95% Confidence Interval')
        
        ax1.set_xlabel('Quality Score Threshold', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Estimated AA Mutations per Gene', fontsize=12, fontweight='bold', color=color1)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1.set_title('AA Mutations per Gene vs Quality Score Filter (with Error Propagation)', 
                     fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(frameon=False, fontsize=10)
        
        # Add data point labels
        for i, (q, aa_mut, aa_err) in enumerate(zip(quality_thresholds, aa_mutations, aa_errors)):
            ax1.annotate(f'Q{q}\n{aa_mut:.3f}±{aa_err:.3f}', 
                        (q, aa_mut), xytext=(5, 5), 
                        textcoords='offset points', fontsize=8, alpha=0.8, color=color1)
        
        # Bottom plot: Mappable bases and AA mutations per gene
        color2 = '#FF6B6B'
        color3 = '#4169E1'
        
        # Mappable bases (left y-axis)
        ax2_twin = ax2.twinx()
        ax2_twin.scatter(quality_thresholds, mappable_bases, 
                        s=100, alpha=0.7, color=color2, edgecolors='black', 
                        linewidth=1, marker='s', label='Mappable Bases')
        ax2_twin.set_ylabel('Number of Mappable Bases', fontsize=12, fontweight='bold', color=color2)
        ax2_twin.tick_params(axis='y', labelcolor=color2)
        
        # AA mutations per gene with error bars (right y-axis)
        ax2.errorbar(quality_thresholds, aa_mutations, yerr=aa_errors,
                    fmt='^', capsize=5, capthick=2, markersize=8,
                    color=color3, ecolor=color3, alpha=0.8, label='AA Mutations ± Error')
        ax2.set_ylabel('Estimated AA Mutations per Gene', fontsize=12, fontweight='bold', color=color3)
        ax2.tick_params(axis='y', labelcolor=color3)
        ax2.set_xlabel('Quality Score Threshold', fontsize=12, fontweight='bold')
        ax2.set_title('Mappable Bases and AA Mutations per Gene vs Quality Score Filter', 
                     fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        # Add legends
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2_twin.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper right', frameon=False, fontsize=10)
        
        # Add data point labels for mappable bases
        for i, (q, bases) in enumerate(zip(quality_thresholds, mappable_bases)):
            ax2_twin.annotate(f'{bases}', (q, bases), xytext=(5, -15), 
                             textcoords='offset points', fontsize=8, alpha=0.8, color=color2)
        
        plt.tight_layout()
        
        # Save the comprehensive plot
        qc_plot_path = os.path.join(results_dir, "comprehensive_qc_analysis.png")
        fig.savefig(qc_plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        logging.info(f"Comprehensive QC plot saved to: {qc_plot_path}")
        
        # Create error analysis plot
        create_error_analysis_plot(quality_thresholds, qc_results, results_dir)
        
        # Save comprehensive data as CSV
        qc_data_path = os.path.join(results_dir, "comprehensive_qc_data.csv")
        with open(qc_data_path, 'w') as f:
            f.write("quality_threshold,mean_aa_mutations,std_aa_mutations,ci_lower,ci_upper,")
            f.write("mappable_bases,hit_rate,hit_rate_ci_lower,hit_rate_ci_upper,")
            f.write("bg_rate,bg_rate_ci_lower,bg_rate_ci_upper,net_rate,net_rate_error,")
            f.write("lambda_bp,lambda_error,alignment_error,")
            f.write("hit_qscore_mean,hit_qscore_std,hit_qscore_uncertainty,")
            f.write("bg_qscore_mean,bg_qscore_std,bg_qscore_uncertainty,")
            f.write("hit_weighted_rate,hit_weighted_error,bg_weighted_rate,bg_weighted_error,")
            f.write("net_weighted_rate,net_weighted_error,lambda_bp_weighted,lambda_error_weighted\n")
            
            for q, r in zip(quality_thresholds, qc_results):
                f.write(f"{q},{r['mean_aa_mutations']:.6f},{r['std_aa_mutations']:.6f},")
                f.write(f"{r['ci_lower']:.6f},{r['ci_upper']:.6f},")
                f.write(f"{r['mappable_bases']},{r['hit_rate']:.6f},")
                f.write(f"{r['hit_rate_ci'][0]:.6f},{r['hit_rate_ci'][1]:.6f},")
                f.write(f"{r['bg_rate']:.6f},{r['bg_rate_ci'][0]:.6f},{r['bg_rate_ci'][1]:.6f},")
                f.write(f"{r['net_rate']:.6f},{r['net_rate_error']:.6f},")
                f.write(f"{r['lambda_bp']:.6f},{r['lambda_error']:.6f},{r['alignment_error']:.6f},")
                
                # Q-score information
                hit_qscore_mean = r['hit_qscore_stats']['mean_qscore'] if r['hit_qscore_stats'] else 0.0
                hit_qscore_std = r['hit_qscore_stats']['std_qscore'] if r['hit_qscore_stats'] else 0.0
                bg_qscore_mean = r['bg_qscore_stats']['mean_qscore'] if r['bg_qscore_stats'] else 0.0
                bg_qscore_std = r['bg_qscore_stats']['std_qscore'] if r['bg_qscore_stats'] else 0.0
                
                f.write(f"{hit_qscore_mean:.2f},{hit_qscore_std:.2f},{r['hit_qscore_uncertainty']:.6f},")
                f.write(f"{bg_qscore_mean:.2f},{bg_qscore_std:.2f},{r['bg_qscore_uncertainty']:.6f},")
                f.write(f"{r.get('hit_weighted_rate', 0.0):.6f},{r.get('hit_weighted_error', 0.0):.6f},")
                f.write(f"{r.get('bg_weighted_rate', 0.0):.6f},{r.get('bg_weighted_error', 0.0):.6f},")
                f.write(f"{r.get('net_weighted_rate', 0.0):.6f},{r.get('net_weighted_error', 0.0):.6f},")
                f.write(f"{r.get('lambda_bp_weighted', 0.0):.6f},{r.get('lambda_error_weighted', 0.0):.6f}\n")
        
        logging.info(f"Comprehensive QC data saved to: {qc_data_path}")
        
    except Exception as e:
        logging.error(f"Error creating comprehensive QC plots: {e}")

def create_error_analysis_plot(quality_thresholds, qc_results, results_dir):
    """
    Create a detailed error analysis plot showing different sources of uncertainty.
    
    Args:
        quality_thresholds: List of quality score thresholds
        qc_results: List of comprehensive analysis results
        results_dir: Directory to save the plot
    """
    try:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Extract error components
        aa_std = [r['std_aa_mutations'] for r in qc_results]
        net_rate_errors = [r['net_rate_error'] for r in qc_results]
        lambda_errors = [r['lambda_error'] for r in qc_results]
        alignment_errors = [r['alignment_error'] for r in qc_results]
        mappable_bases = [r['mappable_bases'] for r in qc_results]
        
        # Plot 1: AA mutation uncertainty vs quality threshold
        ax1.plot(quality_thresholds, aa_std, 'o-', color='#2E8B57', linewidth=2, markersize=6)
        ax1.set_xlabel('Quality Score Threshold')
        ax1.set_ylabel('AA Mutation Standard Deviation')
        ax1.set_title('AA Mutation Uncertainty vs Quality Filter')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Net rate error vs quality threshold
        ax2.plot(quality_thresholds, net_rate_errors, 's-', color='#FF6B6B', linewidth=2, markersize=6)
        ax2.set_xlabel('Quality Score Threshold')
        ax2.set_ylabel('Net Mutation Rate Error')
        ax2.set_title('Net Rate Error vs Quality Filter')
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Lambda error vs quality threshold
        ax3.plot(quality_thresholds, lambda_errors, '^-', color='#4169E1', linewidth=2, markersize=6)
        ax3.set_xlabel('Quality Score Threshold')
        ax3.set_ylabel('Lambda Error (mutations per copy)')
        ax3.set_title('Lambda Error vs Quality Filter')
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Alignment error vs mappable bases
        ax4.scatter(mappable_bases, alignment_errors, s=100, alpha=0.7, color='#FF8C00')
        ax4.set_xlabel('Mappable Bases')
        ax4.set_ylabel('Alignment Error (1/√reads)')
        ax4.set_title('Alignment Error vs Read Count')
        ax4.grid(True, alpha=0.3)
        
        # Add quality threshold labels to scatter plot
        for i, q in enumerate(quality_thresholds):
            ax4.annotate(f'Q{q}', (mappable_bases[i], alignment_errors[i]), 
                        xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        plt.tight_layout()
        
        # Save error analysis plot
        error_plot_path = os.path.join(results_dir, "error_analysis.png")
        fig.savefig(error_plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        logging.info(f"Error analysis plot saved to: {error_plot_path}")
        
    except Exception as e:
        logging.error(f"Error creating error analysis plot: {e}")

def create_qc_plot(quality_thresholds, aa_mutations, mappable_bases, results_dir):
    """
    Create a dual-axis plot showing quality score threshold vs AA mutations per gene and mappable bases.
    
    Args:
        quality_thresholds: List of quality score thresholds
        aa_mutations: List of corresponding AA mutations per gene
        mappable_bases: List of corresponding mappable bases
        results_dir: Directory to save the plot
    """
    try:
        # Create the plot with dual y-axes
        fig, ax1 = plt.subplots(figsize=(12, 8))
        
        # Left y-axis: AA mutations per gene
        color1 = '#2E8B57'
        ax1.set_xlabel('Quality Score Threshold', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Estimated AA Mutations per Gene', fontsize=12, fontweight='bold', color=color1)
        ax1.scatter(quality_thresholds, aa_mutations, 
                   s=100, alpha=0.7, color=color1, edgecolors='black', linewidth=1, label='AA Mutations per Gene')
        ax1.tick_params(axis='y', labelcolor=color1)
        
        # Right y-axis: Mappable bases
        ax2 = ax1.twinx()
        color2 = '#FF6B6B'
        ax2.set_ylabel('Number of Mappable Bases', fontsize=12, fontweight='bold', color=color2)
        ax2.scatter(quality_thresholds, mappable_bases, 
                   s=100, alpha=0.7, color=color2, edgecolors='black', linewidth=1, marker='s', label='Mappable Bases')
        ax2.tick_params(axis='y', labelcolor=color2)
        
        # Customize the plot
        ax1.set_title('AA Mutations per Gene and Mappable Bases vs Quality Score Filter', fontsize=14, fontweight='bold')
        
        # Add grid for better readability
        ax1.grid(True, alpha=0.3)
        
        # Customize ticks and spines
        ax1.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
        ax1.tick_params(axis='both', which='minor', direction='in', length=3)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        
        # Add data point labels for AA mutations
        for i, (q, aa_mut) in enumerate(zip(quality_thresholds, aa_mutations)):
            ax1.annotate(f'Q{q}', (q, aa_mut), xytext=(5, 5), 
                        textcoords='offset points', fontsize=9, alpha=0.8, color=color1)
        
        # Add data point labels for mappable bases
        for i, (q, bases) in enumerate(zip(quality_thresholds, mappable_bases)):
            ax2.annotate(f'{reads}', (q, reads), xytext=(5, -15), 
                        textcoords='offset points', fontsize=8, alpha=0.8, color=color2)
        
        # Add legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', frameon=False, fontsize=10)
        
        # Save the plot
        qc_plot_path = os.path.join(results_dir, "qc_mutation_rate_vs_quality.png")
        fig.savefig(qc_plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        logging.info(f"QC plot saved to: {qc_plot_path}")
        
        # Also save data as CSV for reference
        qc_data_path = os.path.join(results_dir, "qc_mutation_rate_vs_quality.csv")
        with open(qc_data_path, 'w') as f:
            f.write("quality_threshold,aa_mutations_per_gene,mappable_bases\n")
            for q, aa_mut, bases in zip(quality_thresholds, aa_mutations, mappable_bases):
                f.write(f"{q},{aa_mut:.6f},{bases}\n")
        
        logging.info(f"QC data saved to: {qc_data_path}")
        
    except Exception as e:
        logging.error(f"Error creating QC plot: {e}")

def extract_qscores_from_sam(sam_file):
    """
    Extract Q-scores from SAM file and calculate statistics.
    
    Args:
        sam_file: Path to SAM file
    
    Returns:
        dict: Q-score statistics including mean, std, and per-position averages
    """
    try:
        import pysam
        
        qscores = []
        position_qscores = {}  # position -> list of qscores
        
        with pysam.AlignmentFile(sam_file, "r") as samfile:
            for read in samfile:
                if read.is_unmapped:
                    continue
                
                # Get Q-scores for this read
                read_qscores = read.query_qualities
                if read_qscores is None:
                    continue
                
                # Convert to Q-score values (Phred+33 encoding)
                q_values = [q + 33 for q in read_qscores]
                qscores.extend(q_values)
                
                # Store per-position Q-scores
                for i, q_val in enumerate(q_values):
                    pos = read.reference_start + i
                    if pos not in position_qscores:
                        position_qscores[pos] = []
                    position_qscores[pos].append(q_val)
        
        if not qscores:
            return {
                'mean_qscore': 0.0,
                'std_qscore': 0.0,
                'min_qscore': 0.0,
                'max_qscore': 0.0,
                'position_avg_qscores': {},
                'total_bases': 0
            }
        
        # Calculate statistics
        mean_qscore = np.mean(qscores)
        std_qscore = np.std(qscores)
        min_qscore = np.min(qscores)
        max_qscore = np.max(qscores)
        
        # Calculate per-position average Q-scores
        position_avg_qscores = {}
        for pos, pos_qscores in position_qscores.items():
            position_avg_qscores[pos] = np.mean(pos_qscores)
        
        return {
            'mean_qscore': mean_qscore,
            'std_qscore': std_qscore,
            'min_qscore': min_qscore,
            'max_qscore': max_qscore,
            'position_avg_qscores': position_avg_qscores,
            'total_bases': len(qscores)
        }
        
    except Exception as e:
        logging.error(f"Error extracting Q-scores from {sam_file}: {e}")
        return {
            'mean_qscore': 0.0,
            'std_qscore': 0.0,
            'min_qscore': 0.0,
            'max_qscore': 0.0,
            'position_avg_qscores': {},
            'total_bases': 0
        }

def qscore_uncertainty_factor(qscore):
    """
    Convert Q-score to uncertainty factor.
    
    Args:
        qscore: Q-score value (typically 0-40)
    
    Returns:
        float: Uncertainty factor (0-1, where 1 = maximum uncertainty)
    """
    if qscore <= 0:
        return 1.0  # Maximum uncertainty
    
    # Q-score = -10 * log10(P_error)
    # P_error = 10^(-Q/10)
    # Uncertainty factor = sqrt(P_error) for error propagation
    error_probability = 10**(-qscore/10)
    uncertainty_factor = np.sqrt(error_probability)
    
    return uncertainty_factor

def segment_fastq_file(input_fastq, n_segments=10):
    """
    Segment a FASTQ file into N parts for error estimation.
    
    Args:
        input_fastq: Path to input FASTQ.gz file
        n_segments: Number of segments to create
    
    Returns:
        list: Paths to segmented FASTQ files
    """
    try:
        import gzip
        from itertools import cycle
        
        # Create output directory
        base_name = os.path.splitext(os.path.basename(input_fastq))[0].replace('.fastq', '')
        segment_dir = os.path.join(os.path.dirname(input_fastq), f"{base_name}_segments")
        os.makedirs(segment_dir, exist_ok=True)
        
        # Open output files
        segment_files = []
        file_handles = []
        
        for i in range(n_segments):
            segment_path = os.path.join(segment_dir, f"{base_name}_segment_{i+1}.fastq.gz")
            segment_files.append(segment_path)
            file_handles.append(gzip.open(segment_path, 'wt'))
        
        # Read and distribute reads
        read_count = 0
        with gzip.open(input_fastq, 'rt') as infile:
            current_read = []
            
            for line in infile:
                current_read.append(line)
                
                # Complete read (4 lines)
                if len(current_read) == 4:
                    # Write to current segment
                    segment_idx = read_count % n_segments
                    for line in current_read:
                        file_handles[segment_idx].write(line)
                    
                    read_count += 1
                    current_read = []
        
        # Close all files
        for fh in file_handles:
            fh.close()
        
        logging.info(f"Segmented {read_count} reads into {n_segments} files")
        return segment_files
        
    except Exception as e:
        logging.error(f"Error segmenting FASTQ file: {e}")
        return []

def run_segmented_analysis(segment_files, quality_threshold, work_dir, ref_hit_fasta, plasmid_fasta):
    """
    Run mutation rate analysis on each segment and calculate empirical error.
    
    Args:
        segment_files: List of segmented FASTQ files
        quality_threshold: Quality score threshold
        work_dir: Working directory
        ref_hit_fasta: Path to reference hit FASTA
        plasmid_fasta: Path to plasmid FASTA
    
    Returns:
        dict: Results with empirical error estimates
    """
    try:
        segment_results = []
        
        for i, segment_file in enumerate(segment_files):
            logging.info(f"Processing segment {i+1}/{len(segment_files)}")
            
            # Filter segment with NanoFilt
            filtered_segment = os.path.join(work_dir, f"segment_{i+1}_q{quality_threshold}.fastq.gz")
            if not run_nanofilt_filtering(segment_file, quality_threshold, filtered_segment):
                logging.warning(f"Failed to filter segment {i+1}")
                continue
            
            # Load sequences
            hit_seq, hit_id = load_single_sequence(ref_hit_fasta)
            plasmid_seq, plasmid_id = load_single_sequence(plasmid_fasta)
            
            # Find hit region in plasmid
            idx = plasmid_seq.upper().find(hit_seq.upper())
            if idx == -1:
                logging.error(f"Gene region not found in plasmid for segment {i+1}")
                continue
            
            # Align filtered reads to hit region
            sam_hit = run_minimap2(filtered_segment, ref_hit_fasta, f"hit_segment_{i+1}_q{quality_threshold}", work_dir)
            
            # Align filtered reads to full plasmid for background calculation
            sam_plasmid = run_minimap2(filtered_segment, plasmid_fasta, f"plasmid_segment_{i+1}_q{quality_threshold}", work_dir)
            
            # Calculate background rate from full plasmid alignment, excluding target region
            bg_mis, bg_cov, bg_reads = calculate_background_from_plasmid(sam_plasmid, plasmid_seq, idx, len(hit_seq))
            
            # Calculate hit region mutation rate
            mismatch_hit = compute_mismatch_stats_sam(sam_hit, {hit_id: hit_seq})
            hit_info = mismatch_hit[hit_id]
            hit_mis = hit_info["total_mismatches"]
            hit_cov = hit_info["total_covered_bases"]
            
            # Calculate rates
            hit_rate = hit_mis / hit_cov if hit_cov > 0 else 0
            bg_rate = bg_mis / bg_cov if bg_cov > 0 else 0
            net_rate = max(hit_rate - bg_rate, 0.0)
            
            # Calculate AA mutations per gene (simplified)
            lambda_bp = net_rate * len(hit_seq)
            aa_mutations = lambda_bp / 3.0  # Approximate: 3 bp per AA
            
            segment_results.append({
                'segment': i+1,
                'hit_rate': hit_rate,
                'bg_rate': bg_rate,
                'net_rate': net_rate,
                'aa_mutations': aa_mutations,
                'mappable_bases': hit_cov,
                'hit_mismatches': hit_mis,
                'hit_coverage': hit_cov
            })
        
        if not segment_results:
            return None
        
        # Calculate empirical statistics
        aa_mutations_list = [r['aa_mutations'] for r in segment_results]
        net_rates_list = [r['net_rate'] for r in segment_results]
        mappable_bases_list = [r['mappable_bases'] for r in segment_results]
        
        mean_aa = np.mean(aa_mutations_list)
        std_aa = np.std(aa_mutations_list, ddof=1)  # Sample standard deviation
        mean_net_rate = np.mean(net_rates_list)
        std_net_rate = np.std(net_rates_list, ddof=1)
        total_mappable_bases = sum(mappable_bases_list)
        
        # Calculate confidence interval using t-distribution
        n_segments = len(segment_results)
        if n_segments > 1:
            # 95% confidence interval
            from scipy.stats import t
            t_val = t.ppf(0.975, n_segments - 1)
            se_aa = std_aa / np.sqrt(n_segments)
            ci_lower = mean_aa - t_val * se_aa
            ci_upper = mean_aa + t_val * se_aa
        else:
            ci_lower = mean_aa
            ci_upper = mean_aa
        
        return {
            'mean_aa_mutations': mean_aa,
            'std_aa_mutations': std_aa,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'mean_net_rate': mean_net_rate,
            'std_net_rate': std_net_rate,
            'total_mappable_bases': total_mappable_bases,
            'n_segments': n_segments,
            'segment_results': segment_results,
            'quality_threshold': quality_threshold
        }
        
    except Exception as e:
        logging.error(f"Error in segmented analysis: {e}")
        return None
    """
    Calculate mismatches weighted by Q-score uncertainty with proper sampling error.
    
    Args:
        sam_file: Path to SAM file
        ref_seq: Reference sequence
        qscore_stats: Q-score statistics from extract_qscores_from_sam
    
    Returns:
        tuple: (weighted_mismatches, total_weighted_coverage, raw_mismatches, raw_coverage, position_weights, position_outcomes)
    """
    try:
        import pysam
        
        weighted_mismatches = 0.0
        total_weighted_coverage = 0.0
        raw_mismatches = 0
        raw_coverage = 0
        
        # Store position-level data for proper sampling error calculation
        position_weights = []
        position_outcomes = []
        
        position_qscores = qscore_stats['position_avg_qscores']
        
        with pysam.AlignmentFile(sam_file, "r") as samfile:
            for read in samfile:
                if read.is_unmapped:
                    continue
                
                # Get aligned pairs (read_pos, ref_pos)
                for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                    if ref_pos is None or read_pos is None:
                        continue
                    
                    if ref_pos >= len(ref_seq):
                        continue
                    
                    # Get base calls
                    read_base = read.query_sequence[read_pos].upper()
                    ref_base = ref_seq[ref_pos].upper()
                    
                    # Skip if either base is N
                    if read_base == 'N' or ref_base == 'N':
                        continue
                    
                    # Get Q-score for this position
                    qscore = position_qscores.get(ref_pos, qscore_stats['mean_qscore'])
                    uncertainty_factor = qscore_uncertainty_factor(qscore)
                    
                    # Weight by uncertainty (lower Q-score = higher uncertainty = lower weight)
                    weight = 1.0 - uncertainty_factor
                    
                    # Store position-level data
                    position_weights.append(weight)
                    position_outcomes.append(1 if read_base != ref_base else 0)
                    
                    # Count coverage
                    total_weighted_coverage += weight
                    raw_coverage += 1
                    
                    # Count mismatches
                    if read_base != ref_base:
                        weighted_mismatches += weight
                        raw_mismatches += 1
        
        return weighted_mismatches, total_weighted_coverage, raw_mismatches, raw_coverage, position_weights, position_outcomes
        
    except Exception as e:
        logging.error(f"Error calculating Q-score weighted mismatches: {e}")
        return 0.0, 0.0, 0, 0, [], []

def calculate_weighted_sampling_error(position_weights, position_outcomes):
    """
    Calculate proper weighted sampling error from position-level data.
    
    Args:
        position_weights: List of weights for each position
        position_outcomes: List of outcomes (0=match, 1=mismatch) for each position
    
    Returns:
        tuple: (weighted_rate, weighted_error)
    """
    if not position_weights or len(position_weights) == 0:
        return 0.0, 0.0
    
    position_weights = np.array(position_weights)
    position_outcomes = np.array(position_outcomes)
    
    # Calculate weighted rate
    weighted_mismatches = np.sum(position_weights * position_outcomes)
    weighted_coverage = np.sum(position_weights)
    
    if weighted_coverage == 0:
        return 0.0, 0.0
    
    weighted_rate = weighted_mismatches / weighted_coverage
    
    # Proper weighted sampling error calculation
    # Var(p̂) = (1/W²) * Σ[w_i² * (y_i - p̂)²]
    # where W = Σw_i, y_i = outcome, p̂ = weighted rate
    
    residuals = position_outcomes - weighted_rate
    weighted_residuals_squared = position_weights**2 * residuals**2
    weighted_error = np.sqrt(np.sum(weighted_residuals_squared) / (weighted_coverage**2))
    
    return weighted_rate, weighted_error

def calculate_qscore_weighted_error_propagation(weighted_mismatches, weighted_coverage, qscore_stats):
    """
    Calculate error propagation for Q-score weighted mutation rates using proper weighted sampling theory.
    
    Args:
        weighted_mismatches: Q-score weighted mismatch count
        weighted_coverage: Q-score weighted coverage
        qscore_stats: Q-score statistics
    
    Returns:
        tuple: (weighted_rate, weighted_error)
    """
    if weighted_coverage == 0:
        return 0.0, 0.0
    
    weighted_rate = weighted_mismatches / weighted_coverage
    
    # Proper weighted sampling error calculation
    # For weighted binomial: Var(p̂) ≈ (1/n²) * Σ[w_i² * p * (1-p)]
    # where n = Σw_i (weighted_coverage)
    
    # Estimate weight variance from Q-score statistics
    mean_qscore = qscore_stats['mean_qscore']
    std_qscore = qscore_stats['std_qscore']
    
    # Convert Q-score statistics to weight statistics
    mean_weight = 1.0 - qscore_uncertainty_factor(mean_qscore)
    
    # Approximate weight variance using delta method
    # If w = 1 - sqrt(10^(-Q/10)), then Var(w) ≈ Var(Q) * (dw/dQ)²
    # dw/dQ = (ln(10)/20) * 10^(-Q/10) * (1/sqrt(10^(-Q/10)))
    #       = (ln(10)/20) * sqrt(10^(-Q/10))
    
    if mean_qscore > 0:
        # Delta method approximation for weight variance
        error_prob = 10**(-mean_qscore/10)
        weight_derivative = (np.log(10)/20) * np.sqrt(error_prob)
        weight_variance = (std_qscore**2) * (weight_derivative**2)
    else:
        weight_variance = 0.0
    
    # Effective sample size for weighted sampling
    # n_eff = (Σw_i)² / Σ(w_i²) ≈ (Σw_i)² / [n * (E[w]² + Var[w])]
    n_positions = qscore_stats.get('total_bases', weighted_coverage)
    if n_positions > 0:
        expected_w_squared = mean_weight**2 + weight_variance
        effective_n = (weighted_coverage**2) / (n_positions * expected_w_squared)
    else:
        effective_n = weighted_coverage / mean_weight
    
    # Weighted sampling error
    # Var(p̂) = p(1-p) / n_eff * [1 + (Var[w]/E[w]²)]
    weight_cv_squared = weight_variance / (mean_weight**2) if mean_weight > 0 else 0
    weighted_error = np.sqrt(weighted_rate * (1 - weighted_rate) / effective_n * (1 + weight_cv_squared))
    
    return weighted_rate, weighted_error

def binomial_confidence_interval(successes, trials, confidence=0.95):
    """
    Calculate confidence interval for binomial proportion using beta distribution.
    
    Args:
        successes: Number of successes
        trials: Number of trials
        confidence: Confidence level (default 0.95 for 95% CI)
    
    Returns:
        tuple: (lower_bound, upper_bound)
    """
    if not HAVE_SCIPY:
        # Simple normal approximation if scipy not available
        p = successes / trials if trials > 0 else 0
        se = np.sqrt(p * (1 - p) / trials) if trials > 0 else 0
        z = norm.ppf(1 - (1 - confidence) / 2)
        return max(0, p - z * se), min(1, p + z * se)
    
    alpha = 1 - confidence
    lower = beta.ppf(alpha/2, successes, trials - successes + 1) if trials > successes else 0
    upper = beta.ppf(1 - alpha/2, successes + 1, trials - successes) if successes > 0 else 0
    return lower, upper

def propagate_mutation_rate_error(hit_mis, hit_cov, bg_mis, bg_cov, hit_qscore_stats=None, bg_qscore_stats=None):
    """
    Calculate error propagation for net mutation rate = hit_rate - bg_rate, including Q-score uncertainty.
    
    Args:
        hit_mis, hit_cov: Hit region mismatches and coverage
        bg_mis, bg_cov: Background mismatches and coverage
        hit_qscore_stats: Q-score statistics for hit region (optional)
        bg_qscore_stats: Q-score statistics for background region (optional)
    
    Returns:
        tuple: (net_rate, net_rate_error)
    """
    if hit_cov == 0 or bg_cov == 0:
        return 0.0, 0.0
    
    hit_rate = hit_mis / hit_cov
    bg_rate = bg_mis / bg_cov
    
    # Binomial standard errors
    hit_se = np.sqrt(hit_rate * (1 - hit_rate) / hit_cov)
    bg_se = np.sqrt(bg_rate * (1 - bg_rate) / bg_cov)
    
    # Add Q-score uncertainty if available
    if hit_qscore_stats:
        hit_qscore_uncertainty = qscore_uncertainty_factor(hit_qscore_stats['mean_qscore'])
        hit_se = np.sqrt(hit_se**2 + hit_qscore_uncertainty**2)
    
    if bg_qscore_stats:
        bg_qscore_uncertainty = qscore_uncertainty_factor(bg_qscore_stats['mean_qscore'])
        bg_se = np.sqrt(bg_se**2 + bg_qscore_uncertainty**2)
    
    # Net rate and error propagation
    net_rate = max(hit_rate - bg_rate, 0.0)
    net_se = np.sqrt(hit_se**2 + bg_se**2)
    
    return net_rate, net_se

def simulate_aa_distribution_with_error(lambda_bp, lambda_error, cds_seq, n_trials=10000):
    """
    Enhanced Monte Carlo simulation that includes uncertainty in lambda_bp.
    
    Args:
        lambda_bp: Mean mutations per copy (basepairs)
        lambda_error: Standard error of lambda_bp
        cds_seq: Coding sequence
        n_trials: Number of Monte Carlo trials
    
    Returns:
        tuple: (mean_aa_mutations, std_aa_mutations, aa_distribution)
    """
    prot_orig = str(Seq(cds_seq).translate(to_stop=False))
    aa_diffs = []
    
    for _ in range(n_trials):
        # Sample lambda from normal distribution with error
        lambda_sample = np.random.normal(lambda_bp, lambda_error)
        lambda_sample = max(lambda_sample, 0)  # Ensure non-negative
        
        # Number of base changes in this trial ~ Poisson(lambda_sample)
        n_bp_mut = np.random.poisson(lambda_sample)
        
        # Make a mutable copy of the CDS
        seq_list = list(cds_seq.upper())
        
        # Introduce exactly n_bp_mut random single‐base substitutions
        for _ in range(n_bp_mut):
            pos = random.randrange(len(seq_list))
            orig_base = seq_list[pos]
            bases = ["A", "T", "C", "G"]
            bases.remove(orig_base)
            seq_list[pos] = random.choice(bases)
        
        # Translate mutated sequence (no early stop)
        mutated_prot = str(Seq("".join(seq_list)).translate(to_stop=False))
        
        # Count how many amino acids differ
        aa_diff = sum(1 for a, b in zip(prot_orig, mutated_prot) if a != b)
        aa_diffs.append(aa_diff)
    
    mean_aa = np.mean(aa_diffs)
    std_aa = np.std(aa_diffs)
    
    return mean_aa, std_aa, aa_diffs

def bootstrap_aa_mutations(hit_mis, hit_cov, bg_mis, bg_cov, cds_seq, n_bootstrap=1000):
    """
    Bootstrap resampling to estimate confidence intervals for AA mutations.
    
    Args:
        hit_mis, hit_cov: Hit region mismatches and coverage
        bg_mis, bg_cov: Background mismatches and coverage
        cds_seq: Coding sequence
        n_bootstrap: Number of bootstrap samples
    
    Returns:
        tuple: (mean_aa_mutations, ci_lower, ci_upper, bootstrap_distribution)
    """
    bootstrap_results = []
    
    for _ in range(n_bootstrap):
        # Resample reads with replacement (binomial resampling)
        hit_mis_boot = np.random.binomial(hit_cov, hit_mis/hit_cov) if hit_cov > 0 else 0
        bg_mis_boot = np.random.binomial(bg_cov, bg_mis/bg_cov) if bg_cov > 0 else 0
        
        # Calculate net rate
        hit_rate_boot = hit_mis_boot / hit_cov if hit_cov > 0 else 0
        bg_rate_boot = bg_mis_boot / bg_cov if bg_cov > 0 else 0
        net_rate_boot = max(hit_rate_boot - bg_rate_boot, 0)
        
        # Calculate AA mutations
        lambda_bp_boot = net_rate_boot * len(cds_seq)
        
        # Quick simulation for bootstrap (fewer trials for speed)
        aa_mut_boot = simulate_aa_distribution(lambda_bp_boot, cds_seq, n_trials=1000)
        bootstrap_results.append(np.mean(aa_mut_boot))
    
    mean_aa = np.mean(bootstrap_results)
    # Use proper percentile calculation for 95% CI
    ci_lower = np.percentile(bootstrap_results, 2.5)
    ci_upper = np.percentile(bootstrap_results, 97.5)
    
    # Additional validation: ensure CI makes sense
    if ci_lower > mean_aa or ci_upper < mean_aa:
        logging.warning(f"Bootstrap CI validation failed: mean={mean_aa:.4f}, CI=[{ci_lower:.4f}, {ci_upper:.4f}]")
        # Use empirical CI if percentile method fails
        sorted_results = np.sort(bootstrap_results)
        n = len(sorted_results)
        ci_lower = sorted_results[int(0.025 * n)]
        ci_upper = sorted_results[int(0.975 * n)]
    
    return mean_aa, ci_lower, ci_upper, bootstrap_results

def comprehensive_aa_mutation_analysis(hit_mis, hit_cov, bg_mis, bg_cov, cds_seq, 
                                       quality_threshold=None, n_trials=10000,
                                       hit_qscore_stats=None, bg_qscore_stats=None,
                                       sam_hit=None, sam_plasmid=None, hit_seq=None, plasmid_seq=None):
    """
    Comprehensive AA mutation analysis with full error propagation including Q-score uncertainty.
    
    Args:
        hit_mis, hit_cov: Hit region mismatches and coverage
        bg_mis, bg_cov: Background mismatches and coverage
        cds_seq: Coding sequence
        quality_threshold: Quality threshold for logging
        n_trials: Number of Monte Carlo trials
        hit_qscore_stats: Q-score statistics for hit region (optional)
        bg_qscore_stats: Q-score statistics for background region (optional)
    
    Returns:
        dict: Comprehensive results with all error estimates including Q-score effects
    """
    logging.info(f"=== COMPREHENSIVE ERROR MODEL ANALYSIS (Q{quality_threshold}) ===")
    
    # 1. Binomial confidence intervals for mutation rates
    logging.info("1. Calculating binomial confidence intervals for mutation rates...")
    hit_rate_ci = binomial_confidence_interval(hit_mis, hit_cov)
    bg_rate_ci = binomial_confidence_interval(bg_mis, bg_cov)
    logging.info(f"   Hit rate CI: [{hit_rate_ci[0]:.6f}, {hit_rate_ci[1]:.6f}]")
    logging.info(f"   Background rate CI: [{bg_rate_ci[0]:.6f}, {bg_rate_ci[1]:.6f}]")
    
    # 2. Error propagation for net mutation rate (including Q-score uncertainty)
    logging.info("2. Propagating errors for net mutation rate (including Q-score uncertainty)...")
    net_rate, net_rate_error = propagate_mutation_rate_error(
        hit_mis, hit_cov, bg_mis, bg_cov, hit_qscore_stats, bg_qscore_stats
    )
    logging.info(f"   Net mutation rate: {net_rate:.6f} ± {net_rate_error:.6f}")
    
    # 3. Calculate lambda_bp with error
    logging.info("3. Calculating lambda_bp (mutations per copy) with error propagation...")
    lambda_bp = net_rate * len(cds_seq)
    lambda_error = net_rate_error * len(cds_seq)
    logging.info(f"   Lambda_bp: {lambda_bp:.6f} ± {lambda_error:.6f} mutations per copy")
    
    # 4. Q-score weighted analysis
    logging.info("4. Calculating Q-score weighted mutation rates...")
    hit_weighted_mis, hit_weighted_cov, hit_raw_mis, hit_raw_cov, hit_weights, hit_outcomes = calculate_qscore_weighted_mismatches(
        sam_hit, hit_seq, hit_qscore_stats
    )
    bg_weighted_mis, bg_weighted_cov, bg_raw_mis, bg_raw_cov, bg_weights, bg_outcomes = calculate_qscore_weighted_mismatches(
        sam_plasmid, plasmid_seq, bg_qscore_stats
    )
    
    # Calculate proper weighted sampling errors
    hit_weighted_rate, hit_weighted_error = calculate_weighted_sampling_error(hit_weights, hit_outcomes)
    bg_weighted_rate, bg_weighted_error = calculate_weighted_sampling_error(bg_weights, bg_outcomes)
    
    # Net weighted rate
    net_weighted_rate = max(hit_weighted_rate - bg_weighted_rate, 0.0)
    net_weighted_error = np.sqrt(hit_weighted_error**2 + bg_weighted_error**2)
    
    logging.info(f"   Hit weighted rate: {hit_weighted_rate:.6f} ± {hit_weighted_error:.6f}")
    logging.info(f"   Background weighted rate: {bg_weighted_rate:.6f} ± {bg_weighted_error:.6f}")
    logging.info(f"   Net weighted rate: {net_weighted_rate:.6f} ± {net_weighted_error:.6f}")
    
    # 5. Calculate AA mutations per gene (simplified - no Monte Carlo)
    logging.info("5. Calculating AA mutations per gene from weighted rates...")
    lambda_bp_weighted = net_weighted_rate * len(cds_seq)
    lambda_error_weighted = net_weighted_error * len(cds_seq)
    
    # Simple AA mutation estimate (mean of Poisson distribution)
    mean_aa = lambda_bp_weighted / 3.0  # Approximate: 3 bp per AA
    std_aa = np.sqrt(lambda_bp_weighted) / 3.0  # Standard deviation of Poisson
    
    logging.info(f"   Lambda_bp (weighted): {lambda_bp_weighted:.6f} ± {lambda_error_weighted:.6f}")
    logging.info(f"   AA mutations per gene: {mean_aa:.4f} ± {std_aa:.4f}")
    
    # 6. Bootstrap confidence intervals
    logging.info("6. Calculating bootstrap confidence intervals (1,000 resamples)...")
    bootstrap_mean, ci_lower, ci_upper, bootstrap_dist = bootstrap_aa_mutations(
        hit_mis, hit_cov, bg_mis, bg_cov, cds_seq
    )
    logging.info(f"   Bootstrap 95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
    
    # 7. Alignment error estimation
    logging.info("7. Calculating alignment error estimation...")
    alignment_error = 1.0 / np.sqrt(hit_cov) if hit_cov > 0 else 1.0
    logging.info(f"   Alignment error: {alignment_error:.6f}")
    
    # 8. Q-score uncertainty factors
    logging.info("8. Calculating Q-score uncertainty factors...")
    hit_qscore_uncertainty = qscore_uncertainty_factor(hit_qscore_stats['mean_qscore']) if hit_qscore_stats else 0.0
    bg_qscore_uncertainty = qscore_uncertainty_factor(bg_qscore_stats['mean_qscore']) if bg_qscore_stats else 0.0
    logging.info(f"   Hit Q-score uncertainty: {hit_qscore_uncertainty:.6f}")
    logging.info(f"   Background Q-score uncertainty: {bg_qscore_uncertainty:.6f}")
    
    results = {
        'mean_aa_mutations': mean_aa,
        'std_aa_mutations': std_aa,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'hit_rate': hit_mis / hit_cov if hit_cov > 0 else 0,
        'hit_rate_ci': hit_rate_ci,
        'bg_rate': bg_mis / bg_cov if bg_cov > 0 else 0,
        'bg_rate_ci': bg_rate_ci,
        'net_rate': net_rate,
        'net_rate_error': net_rate_error,
        'lambda_bp': lambda_bp,
        'lambda_error': lambda_error,
        'alignment_error': alignment_error,
        'hit_qscore_uncertainty': hit_qscore_uncertainty,
        'bg_qscore_uncertainty': bg_qscore_uncertainty,
        'hit_qscore_stats': hit_qscore_stats,
        'bg_qscore_stats': bg_qscore_stats,
        'bootstrap_distribution': bootstrap_dist,
        'quality_threshold': quality_threshold,
        'mappable_bases': hit_cov,
        # Q-score weighted results
        'hit_weighted_rate': hit_weighted_rate,
        'hit_weighted_error': hit_weighted_error,
        'bg_weighted_rate': bg_weighted_rate,
        'bg_weighted_error': bg_weighted_error,
        'net_weighted_rate': net_weighted_rate,
        'net_weighted_error': net_weighted_error,
        'lambda_bp_weighted': lambda_bp_weighted,
        'lambda_error_weighted': lambda_error_weighted,
        'hit_weighted_mismatches': hit_weighted_mis,
        'hit_weighted_coverage': hit_weighted_cov,
        'bg_weighted_mismatches': bg_weighted_mis,
        'bg_weighted_coverage': bg_weighted_cov
    }
    
    if quality_threshold is not None:
        qscore_info = ""
        if hit_qscore_stats:
            qscore_info = f", hit_qscore={hit_qscore_stats['mean_qscore']:.1f}±{hit_qscore_stats['std_qscore']:.1f}"
        if bg_qscore_stats:
            qscore_info += f", bg_qscore={bg_qscore_stats['mean_qscore']:.1f}±{bg_qscore_stats['std_qscore']:.1f}"
        
        logging.info(f"Quality {quality_threshold}: AA mutations = {mean_aa:.4f} ± {std_aa:.4f} "
                    f"(95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]), "
                    f"mappable_bases={hit_cov}, net_rate={net_rate:.6f}±{net_rate_error:.6f}{qscore_info}")
    
    return results

def simulate_aa_distribution(lambda_bp, cds_seq, n_trials=1000):
    """
    Monte Carlo: each trial draws n_bp_mut ~ Poisson(lambda_bp),
    introduces those random single‐base substitutions, translates,
    and returns a list of amino acid differences per trial.
    """
    prot_orig = str(Seq(cds_seq).translate(to_stop=False))
    aa_diffs  = []

    for _ in range(n_trials):
        # Number of base changes in this trial ~ Poisson(lambda_bp)
        n_bp_mut = np.random.poisson(lambda_bp)

        # Make a mutable copy of the CDS
        seq_list = list(cds_seq.upper())

        # Introduce exactly n_bp_mut random single‐base substitutions
        for _ in range(n_bp_mut):
            pos = random.randrange(len(seq_list))
            orig_base = seq_list[pos]
            bases = ["A", "T", "C", "G"]
            bases.remove(orig_base)
            seq_list[pos] = random.choice(bases)

        # Translate mutated sequence (no early stop)
        mutated_prot = str(Seq("".join(seq_list)).translate(to_stop=False))

        # Count how many amino acids differ
        aa_diff = sum(1 for a, b in zip(prot_orig, mutated_prot) if a != b)
        aa_diffs.append(aa_diff)

    return aa_diffs

def run_main_analysis_for_qscore(fastq_path, qscore, qscore_desc, sample_name, work_dir, results_dir, 
                                 chunks, ref_hit_fasta, plasmid_fasta, hit_seq, hit_id, plasmid_seq, idx):
    """
    Run the main mutation rate analysis for a specific Q-score.
    
    Args:
        fastq_path: Path to the FASTQ file to analyze
        qscore: Q-score threshold (None for unfiltered)
        qscore_desc: Description of the Q-score (e.g., "Q18", "unfiltered")
        sample_name: Name of the sample
        work_dir: Working directory for temporary files
        results_dir: Results directory for output files
        chunks: List of plasmid chunks
        ref_hit_fasta: Path to reference hit FASTA
        plasmid_fasta: Path to plasmid FASTA
        hit_seq: Hit sequence
        hit_id: Hit ID
        plasmid_seq: Plasmid sequence
        idx: Index of hit in plasmid
    
    Returns:
        dict: Analysis results
    """
    logging.info(f"Running main analysis for {qscore_desc}...")
    
    # Ensure work directory exists
    os.makedirs(work_dir, exist_ok=True)
    
    # Create subdirectory for this Q-score analysis
    qscore_results_dir = results_dir
    if qscore is not None:
        qscore_results_dir = os.path.join(results_dir, f"q{qscore}_analysis")
        os.makedirs(qscore_results_dir, exist_ok=True)
    
    # Write chunks FASTA & align to background‐chunks
    chunks_fasta = create_multi_fasta(chunks, work_dir)
    sam_chunks   = run_minimap2(fastq_path, chunks_fasta, "plasmid_chunks_alignment", work_dir)

    # Align to hit (target) alone
    sam_hit = run_minimap2(fastq_path, ref_hit_fasta, "hit_alignment", work_dir)

    # Compute mismatch stats for background chunks (for reference, but not used for background rate)
    chunk_refs      = { f"chunk_{i+1}": seq for i, seq in enumerate(chunks) }
    mismatch_chunks = compute_mismatch_stats_sam(sam_chunks, chunk_refs)

    # ----------------------------
    # COMPUTE BASE DISTRIBUTION AT EACH POSITION OF HIT
    # ----------------------------
    base_counts = [
        {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
        for _ in range(len(hit_seq))
    ]
    samfile_hit = pysam.AlignmentFile(sam_hit, "r")
    for read in samfile_hit.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if read_pos is not None and ref_pos is not None and 0 <= ref_pos < len(hit_seq):
                base = read.query_sequence[read_pos].upper()
                if base not in {"A", "C", "G", "T"}:
                    base = "N"
                base_counts[ref_pos][base] += 1
    samfile_hit.close()

    # ----------------------------
    # ALIGN TO FULL PLASMID TO GET COVERAGE
    # ----------------------------
    sam_plasmid = run_minimap2(fastq_path, plasmid_fasta, "plasmid_full_alignment", work_dir)

    # Calculate plasmid coverage
    plasmid_cov = [0] * len(plasmid_seq)
    samfile_full = pysam.AlignmentFile(sam_plasmid, "r")
    for read in samfile_full.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if ref_pos is not None and 0 <= ref_pos < len(plasmid_seq):
                plasmid_cov[ref_pos] += 1
    samfile_full.close()

    # Calculate background rate from full plasmid alignment, excluding target region
    # This avoids artificial junction mismatches from concatenated chunks
    bg_mis, bg_cov, bg_reads = calculate_background_from_plasmid(sam_plasmid, plasmid_seq, idx, len(hit_seq))
    bg_rate  = (bg_mis / bg_cov) if bg_cov else 0.0       # raw per‐base
    bg_rate_per_kb = bg_rate * 1e3

    logging.info(
        f"Background (plasmid excluding target): total_mismatches={bg_mis}, "
        f"covered_bases={bg_cov}, mapped_reads={bg_reads}, "
        f"rate_per_kb={bg_rate_per_kb:.4f}"
    )

    # Compute mismatch stats for hit (target)
    mismatch_hit = compute_mismatch_stats_sam(sam_hit, {hit_id: hit_seq})
    hit_info     = mismatch_hit[hit_id]
    hit_mis      = hit_info["total_mismatches"]
    hit_cov      = hit_info["total_covered_bases"]
    hit_reads    = hit_info["mapped_reads"]
    hit_rate     = hit_info["avg_mismatch_rate"]     # raw per‐base
    hit_rate_per_kb  = hit_rate * 1e3

    logging.info(
        f"Target ({hit_id}): total_mismatches={hit_mis}, "
        f"covered_bases={hit_cov}, mapped_reads={hit_reads}, "
        f"rate_per_kb={hit_rate_per_kb:.4f}"
    )

    # Two‐proportion Z‐test: is target rate > background rate?
    z_stat, p_val = z_test_two_proportions(hit_mis, hit_cov, bg_mis, bg_cov)

    # Compute "Estimated mutations per target copy (basepairs)" (float)
    length_of_target = len(hit_seq)
    true_diff_rate   = hit_rate - bg_rate
    est_mut_per_copy = max(true_diff_rate * length_of_target, 0.0)

    # Determine if ROI is a valid protein‐coding sequence (updated definition)
    is_protein = True
    reasons    = []
    seq_upper  = hit_seq.upper()
    # Must be multiple of 3
    if len(seq_upper) % 3 != 0:
        is_protein = False
        reasons.append(f"length {len(seq_upper)} is not a multiple of 3")

    if is_protein:
        prot_full = str(Seq(seq_upper).translate(to_stop=False))
        # Check for premature stop codons (anything except possibly at the end)
        if "*" in prot_full[:-1]:
            is_protein = False
            reasons.append("premature stop codon detected before the last codon")
        # No requirement to start with ATG or end with stop beyond the last codon

    # If protein, simulate AA distribution per copy using Poisson sampling
    if is_protein:
        logging.info(f"Simulating amino acid distribution with λ_bp={est_mut_per_copy:.2f}")
        aa_diffs = simulate_aa_distribution(est_mut_per_copy, hit_seq, n_trials=1000)
        avg_aa_mutations = sum(aa_diffs) / len(aa_diffs)
        
        # Log simulation results for debugging
        logging.info(f"AA simulation results: min={min(aa_diffs)}, max={max(aa_diffs)}, mean={avg_aa_mutations:.3f}")
        logging.info(f"AA simulation distribution: {len([x for x in aa_diffs if x == 0])} zeros, {len([x for x in aa_diffs if x > 0])} non-zeros")
    else:
        aa_diffs = []
        avg_aa_mutations = None

    # Update Q-score info for titles
    qscore_info = f" ({qscore_desc})" if qscore_desc != "unfiltered" else ""
    
    # ----------------------------
    # SAVE CSV FOR MUTATION RATES (PANEL 1)
    # ----------------------------
    gene_mismatch_csv = os.path.join(qscore_results_dir, "gene_mismatch_rates.csv")
    with open(gene_mismatch_csv, "w", newline="") as csvfile:
        csvfile.write(f"# gene_id: {hit_id}\n")
        csvfile.write(f"# background_rate_per_kb: {bg_rate_per_kb:.6f}\n")
        csvfile.write("position_1based,mismatch_rate_per_base\n")
        for pos0, rate in enumerate(hit_info["pos_rates"]):
            csvfile.write(f"{pos0 + 1},{rate:.6e}\n")
    logging.info(f"Saved CSV for gene mismatch rates: {gene_mismatch_csv}")

    # ----------------------------
    # SAVE CSV FOR BASE DISTRIBUTION (PANEL 2)
    # ----------------------------
    base_dist_csv = os.path.join(qscore_results_dir, "base_distribution.csv")
    with open(base_dist_csv, "w", newline="") as csvfile:
        csvfile.write(f"# gene_id: {hit_id}\n")
        csvfile.write("position_1based,ref_base,A_count,C_count,G_count,T_count,N_count\n")
        for pos0, counts in enumerate(base_counts):
            ref_base = seq_upper[pos0]
            csvfile.write(f"{pos0 + 1},{ref_base},{counts['A']},{counts['C']},{counts['G']},{counts['T']},{counts['N']}\n")
    logging.info(f"Saved CSV for base distribution: {base_dist_csv}")

    # ----------------------------
    # SAVE CSV FOR AA SUBSTITUTIONS (PANEL 3) - only if protein
    # ----------------------------
    if is_protein:
        aa_subst_csv = os.path.join(qscore_results_dir, "aa_substitutions.csv")
        with open(aa_subst_csv, "w", newline="") as csvfile:
            csvfile.write(f"# gene_id: {hit_id}\n")
            csvfile.write(f"# lambda_bp_mut: {est_mut_per_copy:.6f}\n")
            csvfile.write("position_1based,ref_aa,alt_aa,count\n")
            # This would need to be implemented based on the specific requirements
            # For now, just write a header
        logging.info(f"Saved CSV for AA substitutions: {aa_subst_csv}")

    # ----------------------------
    # SAVE CSV FOR PLASMID COVERAGE (PANEL 4)
    # ----------------------------
    plasmid_cov_csv = os.path.join(qscore_results_dir, "plasmid_coverage.csv")
    with open(plasmid_cov_csv, "w", newline="") as csvfile:
        csvfile.write("position_1based,coverage\n")
        for pos0, cov in enumerate(plasmid_cov):
            csvfile.write(f"{pos0 + 1},{cov}\n")
    logging.info(f"Saved CSV for plasmid coverage: {plasmid_cov_csv}")

    # ----------------------------
    # SAVE CSV FOR AA MUTATION DISTRIBUTION (PANEL 3)
    # ----------------------------
    aa_dist_csv = os.path.join(qscore_results_dir, "aa_mutation_distribution.csv")
    with open(aa_dist_csv, "w", newline="") as csvfile:
        csvfile.write(f"# gene_id: {hit_id}\n")
        csvfile.write(f"# lambda_bp_mut: {est_mut_per_copy:.6f}\n")
        csvfile.write(f"# n_trials: 1000\n")
        if is_protein:
            csvfile.write("trial_index,aa_mutations\n")
            for idx_trial, aa_count in enumerate(aa_diffs, start=1):
                csvfile.write(f"{idx_trial},{aa_count}\n")
        else:
            csvfile.write("# No AA distribution because region is not protein-coding\n")
    logging.info(f"Saved CSV for AA mutation distribution: {aa_dist_csv}")

    # ----------------------------
    # PREPARE PANEL FIGURE WITH 4 SUBPLOTS
    # ----------------------------
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), constrained_layout=True)
    # axes[0,0]: Mutation rate over gene of interest
    # axes[0,1]: Rolling mutation rate across plasmid (20 bp window)
    # axes[1,0]: Coverage of plasmid with ROI shaded
    # axes[1,1]: KDE of AA mutations per copy

    # --- Panel 1: Mutation rate over gene of interest ---
    ax0 = axes[0, 0]
    positions_gene = np.arange(1, len(hit_info["pos_rates"]) + 1)
    ax0.axhspan(0, bg_rate, color='gray', alpha=0.3, label="Background rate")
    ax0.plot(positions_gene, hit_info["pos_rates"],
             color="#2E86AB", linestyle='-', linewidth=1.5, alpha=0.8,
             label="Mutation rate")
    ax0.set_title(f"Mismatch Rate per Position: Gene of Interest{qscore_info}", fontsize=14, fontweight='bold')
    ax0.set_xlabel("Position in Gene (bp)", fontsize=12)
    ax0.set_ylabel("Mismatch Rate", fontsize=12)
    ax0.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
    ax0.tick_params(axis='both', which='minor', direction='in', length=3)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.legend(loc="upper right", frameon=False, fontsize=10)

    # --- Panel 2: Rolling mutation rate across plasmid ---
    ax1 = axes[0, 1]
    # Calculate rolling mutation rate across plasmid
    window_size = 20
    rolling_positions = []
    rolling_rates = []
    
    # Calculate mismatches per position across the plasmid
    plasmid_mismatches = [0] * len(plasmid_seq)
    samfile_rolling = pysam.AlignmentFile(sam_plasmid, "r")
    for read in samfile_rolling.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if ref_pos is not None and 0 <= ref_pos < len(plasmid_seq):
                if read_pos is not None:
                    read_base = read.query_sequence[read_pos].upper()
                    ref_base = plasmid_seq[ref_pos].upper()
                    if read_base != ref_base and read_base in "ACGT" and ref_base in "ACGT":
                        plasmid_mismatches[ref_pos] += 1
    samfile_rolling.close()
    
    # Calculate rolling mutation rate
    for i in range(len(plasmid_cov) - window_size + 1):
        window_cov = plasmid_cov[i:i + window_size]
        window_mismatches = plasmid_mismatches[i:i + window_size]
        
        total_coverage = sum(window_cov)
        total_mismatches = sum(window_mismatches)
        
        if total_coverage > 0:  # Only include windows with coverage
            rolling_positions.append(i + window_size // 2)
            mutation_rate = total_mismatches / total_coverage
            rolling_rates.append(mutation_rate)
    
    ax1.plot(rolling_positions, rolling_rates,
             color="#FF6B6B", linestyle='-', linewidth=2, alpha=0.8,
             label="Rolling average (20 bp)")
    ax1.set_title(f"Rolling Mutation Rate Across Plasmid (20 bp Window){qscore_info}", fontsize=14, fontweight='bold')
    ax1.set_xlabel("Position on Plasmid (bp)", fontsize=12)
    ax1.set_ylabel("Mismatch Rate", fontsize=12)
    ax1.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
    ax1.tick_params(axis='both', which='minor', direction='in', length=3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.legend(loc="upper right", frameon=False, fontsize=10)
    # Shade the ROI region
    start_roi = idx + 1
    end_roi = idx + len(hit_seq)
    ax1.axvspan(start_roi, end_roi, color='gray', alpha=0.3, label=f"ROI: {start_roi}–{end_roi}")

    # --- Panel 3: Coverage of plasmid with ROI shaded ---
    ax2 = axes[1, 0]
    plasmid_positions = np.arange(1, len(plasmid_cov) + 1)
    ax2.plot(plasmid_positions, plasmid_cov,
             linestyle='-', color='black', linewidth=1.0, alpha=0.8, label="Coverage")
    ax2.set_title(f"Full Plasmid Coverage with ROI Shaded{qscore_info}", fontsize=14, fontweight='bold')
    ax2.set_xlabel("Position on Plasmid", fontsize=12)
    ax2.set_ylabel("Coverage (# reads)", fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
    ax2.tick_params(axis='both', which='minor', direction='in', length=3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    start_roi = idx + 1
    end_roi   = idx + len(hit_seq)
    ax2.axvspan(start_roi, end_roi, color='gray', alpha=0.3, label=f"ROI: {start_roi}–{end_roi}")

    # --- Panel 4: KDE of AA mutations per copy ---
    ax3 = axes[1, 1]
    if is_protein and aa_diffs and len(aa_diffs) > 0:
        x_vals = np.array(aa_diffs)
        unique_vals = np.unique(x_vals)
        
        if len(unique_vals) > 1:
            # Multiple unique values - use KDE or histogram
            if HAVE_SCIPY:
                try:
                    kde = gaussian_kde(x_vals)
                    x_grid = np.linspace(0, max(x_vals), 200)
                    kde_values = kde(x_grid)
                    ax3.plot(x_grid, kde_values,
                             color="#C44E52", linewidth=2.0, alpha=0.8, label="KDE")
                    ax3.fill_between(x_grid, kde_values, color="#C44E52", alpha=0.3)
                    ax3.set_ylim(bottom=0)
                except Exception as e:
                    logging.warning(f"KDE failed: {e}, falling back to histogram")
                    ax3.hist(x_vals, bins=min(20, len(unique_vals)), 
                            color="#C44E52", alpha=0.7, density=True, edgecolor='black')
            else:
                ax3.hist(x_vals, bins=min(20, len(unique_vals)), 
                        color="#C44E52", alpha=0.7, density=True, edgecolor='black')
        else:
            # Single unique value - just show a bar
            ax3.bar(unique_vals, [1.0], color="#C44E52", alpha=0.7, width=0.1)
            ax3.set_xlim(unique_vals[0] - 0.5, unique_vals[0] + 0.5)
    else:
        # Not protein or no AA differences
        ax3.text(0.5, 0.5, "Not a protein‐coding region", 
                 horizontalalignment='center', verticalalignment='center', 
                 fontsize=12, color='gray', transform=ax3.transAxes)
        
        ax3.set_title("AA Mutation Distribution", fontsize=14, fontweight='bold')
        ax3.set_xlabel("Number of AA Mutations", fontsize=12)
        ax3.set_ylabel("Density", fontsize=12)
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.set_xticks([])
        ax3.set_yticks([])

    # Save the combined figure as both PNG and PDF
    panel_path_png = os.path.join(qscore_results_dir, "summary_panels.png")
    panel_path_pdf = os.path.join(qscore_results_dir, "summary_panels.pdf")
    fig.savefig(panel_path_png, dpi=150, transparent=False)
    fig.savefig(panel_path_pdf)  # vector format
    plt.close(fig)
    logging.info(f"Saved combined panel figure as PNG: {panel_path_png}")
    logging.info(f"Saved combined panel figure as PDF: {panel_path_pdf}")

    # ----------------------------
    # COMPUTE MUTATION SPECTRUM FOR ABOVE-BACKGROUND POSITIONS
    # ----------------------------
    # Define categories and reference percentages
    categories = {
        "A→G, T→C":    {"pairs": [("A", "G"), ("T", "C")], "ref": 17.5},
        "G→A, C→T":    {"pairs": [("G", "A"), ("C", "T")], "ref": 25.5},
        "A→T, T→A":    {"pairs": [("A", "T"), ("T", "A")], "ref": 28.5},
        "A→C, T→G":    {"pairs": [("A", "C"), ("T", "G")], "ref": 4.7},
        "G→C, C→G":    {"pairs": [("G", "C"), ("C", "G")], "ref": 4.1},
        "G→T, C→A":    {"pairs": [("G", "T"), ("C", "A")], "ref": 14.1},
    }

    # Tally observed counts at above-background positions
    category_counts = {cat: 0 for cat in categories}
    total_alt_counts = 0

    for pos0, rate in enumerate(hit_info["pos_rates"]):
        if rate <= bg_rate:
            continue
        ref_base = seq_upper[pos0]
        counts   = base_counts[pos0]
        for alt_base in ("A", "C", "G", "T"):
            if alt_base == ref_base:
                continue
            cnt = counts.get(alt_base, 0)
            if cnt == 0:
                continue
            # Determine which category this (ref→alt) belongs to
            for cat, info in categories.items():
                if (ref_base, alt_base) in info["pairs"]:
                    category_counts[cat] += cnt
                    total_alt_counts += cnt
                    break

    # Compute sample percentages
    sample_percent = {}
    if total_alt_counts > 0:
        for cat, cnt in category_counts.items():
            sample_percent[cat] = 100.0 * cnt / total_alt_counts
    else:
        for cat in categories:
            sample_percent[cat] = 0.0

    # ----------------------------
    # GENERATE PDF TABLE (MUTATION SPECTRUM)
    # ----------------------------
    pdf_path = os.path.join(qscore_results_dir, f"{sample_name}_mutation_spectrum.pdf")
    # Prepare table data
    table_rows = []
    for cat in categories:
        ref_pct = categories[cat]["ref"]
        samp_pct = sample_percent[cat]
        table_rows.append([cat, f"{ref_pct:.1f}%", f"{samp_pct:.1f}%"])

    # Create a matplotlib figure for the table
    fig, ax = plt.subplots(figsize=(6, 3))  # adjust size as needed
    ax.axis("off")

    col_labels = ["Mutation Type", "Mutazyme II reference", sample_name]
    tbl = ax.table(
        cellText=table_rows,
        colLabels=col_labels,
        cellLoc="center",
        colLoc="center",
        loc="center"
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1, 1.5)  # stretch rows

    # Add a title
    ax.set_title("Mutation Spectrum (Above-Background Sites)", fontsize=12, fontweight="bold", pad=20)

    # Save as PDF
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved mutation spectrum table as PDF: {pdf_path}")

    # ----------------------------
    # WRITE PER-SAMPLE SUMMARY TXT
    # ----------------------------
    sample_summary_path = os.path.join(qscore_results_dir, "summary.txt")
    with open(sample_summary_path, "w") as txtf:
        txtf.write(f"Sample: {sample_name}\n")
        txtf.write(f"{'=' * (8 + len(sample_name))}\n\n")
        txtf.write("1) Background (plasmid excluding target):\n")
        txtf.write(f"   • Total mismatches:    {bg_mis}\n")
        txtf.write(f"   • Total covered bases: {bg_cov}\n")
        txtf.write(f"   • Mapped reads:        {bg_reads}\n")
        txtf.write(f"   • Rate (per base):     {bg_rate:.6e}\n")
        txtf.write(f"   • Rate (per kb):       {bg_rate_per_kb:.6e}\n\n")

        txtf.write("2) Target (ROI) stats:\n")
        txtf.write(f"   • Gene ID:            {hit_id}\n")
        txtf.write(f"   • Total mismatches:   {hit_mis}\n")
        txtf.write(f"   • Total covered bases:{hit_cov}\n")
        txtf.write(f"   • Mapped reads:       {hit_reads}\n")
        txtf.write(f"   • Rate (per base):    {hit_rate:.6e}\n")
        txtf.write(f"   • Rate (per kb):      {hit_rate_per_kb:.6e}\n")
        txtf.write(f"   • Z‐statistic:        {z_stat:.4f}\n")
        txtf.write(f"   • p‐value:            {p_val if p_val is not None else 'N/A'}\n")
        txtf.write(f"   • Estimated mutations per copy: {est_mut_per_copy:.6e}\n\n")

        txtf.write("3) Protein‐coding evaluation:\n")
        txtf.write(f"   • Is protein: {is_protein}\n")
        if is_protein:
            txtf.write(f"   • Average AA mutations per copy (simulated): {avg_aa_mutations:.3f}\n")
        else:
            txtf.write(f"   • Reason(s): {('; '.join(reasons) if reasons else 'N/A')}\n")
        txtf.write("\n4) Mutation spectrum (above-background sites):\n")
        for cat in categories:
            txtf.write(f"   • {cat}: {sample_percent[cat]:.1f}%  (Ref: {categories[cat]['ref']:.1f}%)\n")
        txtf.write("\n5) Output files written to:\n")
        txtf.write(f"   • {gene_mismatch_csv}\n")
        txtf.write(f"   • {base_dist_csv}\n")
        if is_protein:
            txtf.write(f"   • {aa_subst_csv}\n")
        txtf.write(f"   • {plasmid_cov_csv}\n")
        txtf.write(f"   • {aa_dist_csv}\n")
        txtf.write(f"   • {panel_path_png} (figure)\n")
        txtf.write(f"   • {panel_path_pdf} (figure)\n")
        txtf.write(f"   • {pdf_path} (mutation spectrum table)\n")

    logging.info(f"Wrote per-sample summary to: {sample_summary_path}")
    
    return {
        'qscore': qscore,
        'qscore_desc': qscore_desc,
        'qscore_results_dir': qscore_results_dir,
        'bg_mis': bg_mis,
        'bg_cov': bg_cov,
        'bg_reads': bg_reads,
        'bg_rate': bg_rate,
        'bg_rate_per_kb': bg_rate_per_kb,
        'hit_mis': hit_mis,
        'hit_cov': hit_cov,
        'hit_reads': hit_reads,
        'hit_rate': hit_rate,
        'hit_rate_per_kb': hit_rate_per_kb,
        'hit_info': hit_info,
        'z_stat': z_stat,
        'p_val': p_val,
        'est_mut_per_copy': est_mut_per_copy,
        'is_protein': is_protein,
        'reasons': reasons,
        'aa_diffs': aa_diffs,
        'avg_aa_mutations': avg_aa_mutations,
        'base_counts': base_counts,
        'qscore_info': qscore_info,
        'sam_plasmid': sam_plasmid
    }


def process_single_fastq(fastq_path, master_summary_path):
    """
    Run the entire mutation‐rate analysis pipeline on a single FASTQ file.
    Results and intermediate files are written into a sample‐specific folder.
    A per‐sample summary.txt is created, and its key stats are appended to master_summary.txt.
    Also produces a PDF table of mutation spectrum (sample vs. Mutazyme II reference).
    """
    sample_name = os.path.basename(fastq_path).replace(".fastq.gz", "")
    # Define work and results directories for this sample
    work_dir    = os.path.join(BASE_WORK_DIR, sample_name)
    results_dir = os.path.join(BASE_RESULTS_DIR, sample_name)

    # Clean up old directories if they exist
    if os.path.isdir(work_dir):
        shutil.rmtree(work_dir)
    if os.path.isdir(results_dir):
        shutil.rmtree(results_dir)

    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    # Configure logging into this sample's results folder
    setup_logging(results_dir)
    logging.info(f"--- Starting analysis for sample: {sample_name} ---")

    # Paths to reference FASTA files (shared across samples)
    ref_hit_fasta = os.path.join(DATA_DIR, "region_of_interest.fasta")
    plasmid_fasta = os.path.join(DATA_DIR, "plasmid.fasta")

    # Load sequences (same for every sample)
    hit_seq, hit_id         = load_single_sequence(ref_hit_fasta)
    plasmid_seq, plasmid_id = load_single_sequence(plasmid_fasta)

    logging.info(f"Plasmid length: {len(plasmid_seq)} bp")
    logging.info(f"Gene of interest length: {len(hit_seq)} bp")

    # Remove hit region from plasmid → "background plasmid"
    # Use case-insensitive search to find the hit sequence in the plasmid
    idx = plasmid_seq.upper().find(hit_seq.upper())
    if idx == -1:
        logging.error("Gene region not found in plasmid")
        return
    plasmid_no_gene = plasmid_seq[:idx] + plasmid_seq[idx + len(hit_seq):]
    
    logging.info(f"Gene found at position {idx+1}-{idx+len(hit_seq)} (1-based)")
    logging.info(f"Background region length: {len(plasmid_no_gene)} bp")

    # Chunk the background plasmid into N pieces
    n_chunks = 10
    length   = len(plasmid_no_gene)
    size     = length // n_chunks
    
    # Check if chunks would be too small and adjust accordingly
    # Minimum chunk size ensures sufficient data for reliable alignment and statistics
    # 50 bp minimum provides better alignment performance with minimap2
    min_chunk_size = 50
    if size < min_chunk_size:
        logging.warning(f"Background region ({length} bp) would create chunks smaller than {min_chunk_size} bp. Adjusting chunk count.")
        n_chunks = max(1, length // min_chunk_size)
        size = length // n_chunks
        logging.info(f"Adjusted to {n_chunks} chunks of approximately {size} bp each")
    
    chunks   = [
        plasmid_no_gene[i * size : (length if i == n_chunks - 1 else (i + 1) * size)]
        for i in range(n_chunks)
    ]
    
    # Log chunk sizes for debugging
    chunk_sizes = [len(chunk) for chunk in chunks]
    logging.info(f"Chunk sizes: {chunk_sizes} bp")

    # ----------------------------
    # RUN QC ANALYSIS TO GET Q-SCORE RESULTS
    # ----------------------------
    logging.info("Running QC analysis to get Q-score results...")
    qc_results = None
    optimal_qscore = None
    try:
        qc_results, optimal_qscore = run_qc_analysis(fastq_path, results_dir)
        if qc_results is not None:
            logging.info(f"QC analysis completed successfully. Found {len(qc_results)} Q-score results.")
            if optimal_qscore is not None:
                logging.info(f"Optimal Q-score determined: {optimal_qscore}")
        else:
            logging.warning("QC analysis completed but no Q-score results found.")
    except Exception as e:
        logging.error(f"QC analysis failed: {e}")
        logging.warning("Proceeding with unfiltered data only.")
        # Don't fail the entire analysis if QC fails

    # ----------------------------
    # RUN MAIN ANALYSIS FOR EACH Q-SCORE (INCLUDING UNFILTERED)
    # ----------------------------
    qscores_to_analyze = []
    
    # Add unfiltered analysis
    qscores_to_analyze.append((None, fastq_path, "unfiltered"))
    
    # Add Q-score filtered analyses if QC results are available
    if qc_results is not None:
        for result in qc_results:
            qscore = result['quality_threshold']
            filtered_fastq_path = os.path.join(work_dir, f"{sample_name}_q{qscore}.fastq.gz")
            if run_nanofilt_filtering(fastq_path, qscore, filtered_fastq_path):
                qscores_to_analyze.append((qscore, filtered_fastq_path, f"Q{qscore}"))
                logging.info(f"Successfully created Q{qscore} filtered data for analysis.")
            else:
                logging.warning(f"Failed to create Q{qscore} filtered data.")
    
    logging.info(f"Will run main analysis for {len(qscores_to_analyze)} conditions: {[desc for _, _, desc in qscores_to_analyze]}")
    
    # Run main analysis for each Q-score condition
    analysis_results = []
    for qscore, analysis_fastq_path, qscore_desc in qscores_to_analyze:
        result = run_main_analysis_for_qscore(
            analysis_fastq_path, qscore, qscore_desc, sample_name, work_dir, results_dir,
            chunks, ref_hit_fasta, plasmid_fasta, hit_seq, hit_id, plasmid_seq, idx
        )
        analysis_results.append(result)

    # ----------------------------
    # CLEAN UP: remove this sample's work_dir entirely
    # ----------------------------
    if os.path.isdir(work_dir):
        shutil.rmtree(work_dir)
        logging.info(f"Removed temporary work directory: {work_dir}")

    logging.info(f"--- Finished analysis for sample: {sample_name} ---\n")

def main():
    # Find all FASTQ.GZ files in DATA_DIR
    fastq_files = glob.glob(os.path.join(DATA_DIR, "*.fastq.gz"))
    if not fastq_files:
        print(f"No FASTQ files found in {DATA_DIR}. Exiting.")
        sys.exit(1)

    # Ensure base work/results directories exist
    os.makedirs(BASE_WORK_DIR, exist_ok=True)
    os.makedirs(BASE_RESULTS_DIR, exist_ok=True)

    # Prepare (or overwrite) master summary
    master_summary_path = os.path.join(BASE_RESULTS_DIR, "master_summary.txt")
    with open(master_summary_path, "w") as masterf:
        masterf.write("Sample\tBackground_Rate\tTarget_Rate\tZ_stat\tP_value\tEst_Mut_per_Copy\tIs_Protein\n")

    # Process each FASTQ file separately
    for fastq_path in sorted(fastq_files):
        process_single_fastq(fastq_path, master_summary_path)

if __name__ == "__main__":
    main()
