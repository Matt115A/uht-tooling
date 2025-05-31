#!/usr/bin/env python3

import subprocess
import sys
import glob
import os
import shutil
import logging
import pysam
import csv
import random
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import math

# Use a built-in Matplotlib style ("ggplot") for consistency
plt.style.use("ggplot")

# Try to import scipy for Z‐test p‐value and KDE
try:
    from scipy.stats import norm, gaussian_kde
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
    print("Warning: scipy is not installed. Z-test p-values and KDE may be limited.")

# Define directories
DATA_DIR     = os.path.join("data", "mutation_rate_calculator")
WORK_DIR     = os.path.join(DATA_DIR, "tmp")
RESULTS_DIR  = os.path.join("results", "mutation_rate_calculator")

# Ensure directories exist
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

def setup_logging():
    """
    Configure logging so that INFO and above go into run.log only.
    """
    log_file = os.path.join(RESULTS_DIR, "run.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            # No StreamHandler: INFO logs go into file only
        ]
    )

def run_minimap2(reads_pattern, ref_fasta, out_prefix):
    """
    Runs minimap2 on FASTQ files matching reads_pattern.
    Saves SAM to WORK_DIR/out_prefix.sam.
    """
    fastq_files = glob.glob(reads_pattern)
    if not fastq_files:
        logging.error(f"No FASTQ files found matching {reads_pattern}")
        raise FileNotFoundError(f"No FASTQ files found matching {reads_pattern}")

    sam_output = os.path.join(WORK_DIR, f"{out_prefix}.sam")
    cmd = ["minimap2", "-ax", "map-ont", ref_fasta, *fastq_files]
    logging.info(f"Running minimap2: {' '.join(cmd)} → {sam_output}")
    with open(sam_output, "w") as out_sam:
        subprocess.run(cmd, stdout=out_sam, check=True)
    return sam_output

def load_single_sequence(fasta_path):
    """
    Loads exactly one sequence from a FASTA. Raises if not exactly one record.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) != 1:
        logging.error(f"Expected exactly 1 record in {fasta_path}, found {len(records)}.")
        raise ValueError(f"Expected exactly 1 record in {fasta_path}, found {len(records)}.")
    seq_id = records[0].id
    seq_str = str(records[0].seq)
    logging.info(f"Loaded sequence {seq_id} from {fasta_path}")
    return seq_str, seq_id

def create_multi_fasta(chunks, out_fasta_prefix="plasmid_chunks"):
    """
    Writes each chunk as a separate FASTA entry:
      >chunk_1
      ACTG...
      >chunk_2
      TTAG...
    """
    out_fasta = os.path.join(WORK_DIR, f"{out_fasta_prefix}.fasta")
    logging.info(f"Writing {len(chunks)} chunks to {out_fasta}")
    with open(out_fasta, "w") as f:
        for i, seq in enumerate(chunks):
            f.write(f">chunk_{i+1}\n")
            f.write(str(seq) + "\n")
    return out_fasta

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

def simulate_aa_distribution(lambda_bp, cds_seq, n_trials=1000):
    """
    Monte Carlo: each trial draws n_bp_mut ~ Poisson(lambda_bp),
    introduces those random single‐base substitutions, translates,
    and returns a list of amino acid differences per trial.
    """
    prot_orig = str(Seq(cds_seq).translate(to_stop=False))
    prot_len  = len(prot_orig)
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

def main():
    setup_logging()

    # Input paths
    reads_pattern   = os.path.join(DATA_DIR, "*.fastq.gz")
    ref_hit_fasta   = os.path.join(DATA_DIR, "region_of_interest.fasta")
    plasmid_fasta   = os.path.join(DATA_DIR, "plasmid.fasta")

    logging.info("Starting mutation rate analysis run")
    logging.info(f"Configuration: reads_pattern={reads_pattern}, ref_hit={ref_hit_fasta}, plasmid={plasmid_fasta}")

    # Load sequences
    hit_seq, hit_id         = load_single_sequence(ref_hit_fasta)
    plasmid_seq, plasmid_id = load_single_sequence(plasmid_fasta)

    # Remove hit region from plasmid → “background plasmid”
    idx = plasmid_seq.find(hit_seq)
    if idx == -1:
        logging.error("Gene region not found in plasmid")
        sys.exit(1)
    plasmid_no_gene = plasmid_seq[:idx] + plasmid_seq[idx + len(hit_seq):]

    # Chunk the background plasmid into N pieces
    n_chunks = 10
    length   = len(plasmid_no_gene)
    size     = length // n_chunks
    chunks   = [
        plasmid_no_gene[i * size : (length if i == n_chunks - 1 else (i + 1) * size)]
        for i in range(n_chunks)
    ]

    # Write chunks FASTA & align to background‐chunks
    chunks_fasta     = create_multi_fasta(chunks)
    sam_chunks       = run_minimap2(reads_pattern, chunks_fasta, "plasmid_chunks_alignment")

    # Align to hit (target) alone
    sam_hit          = run_minimap2(reads_pattern, ref_hit_fasta,   "hit_alignment")

    # Compute mismatch stats for background chunks
    chunk_refs       = { f"chunk_{i+1}": seq for i, seq in enumerate(chunks) }
    mismatch_chunks  = compute_mismatch_stats_sam(sam_chunks, chunk_refs)

    # Summarize across all chunks → background totals
    bg_mis   = sum(info["total_mismatches"]    for info in mismatch_chunks.values())
    bg_cov   = sum(info["total_covered_bases"]  for info in mismatch_chunks.values())
    bg_reads = sum(info["mapped_reads"]         for info in mismatch_chunks.values())
    bg_rate  = (bg_mis / bg_cov) if bg_cov else 0.0       # raw per‐base
    bg_rate_per_kb = bg_rate * 1e3

    logging.info(
        f"Background (all chunks): total_mismatches={bg_mis}, "
        f"covered_bases={bg_cov}, mapped_reads={bg_reads}, "
        f"rate_per_kb={bg_rate_per_kb:.4f}"
    )

    # Compute mismatch stats for hit (target)
    mismatch_hit     = compute_mismatch_stats_sam(sam_hit, {hit_id: hit_seq})
    hit_info         = mismatch_hit[hit_id]
    hit_mis          = hit_info["total_mismatches"]
    hit_cov          = hit_info["total_covered_bases"]
    hit_reads        = hit_info["mapped_reads"]
    hit_rate         = hit_info["avg_mismatch_rate"]     # raw per‐base
    hit_rate_per_kb  = hit_rate * 1e3

    logging.info(
        f"Target ({hit_id}): total_mismatches={hit_mis}, "
        f"covered_bases={hit_cov}, mapped_reads={hit_reads}, "
        f"rate_per_kb={hit_rate_per_kb:.4f}"
    )

    # Two‐proportion Z‐test: is target rate > background rate?
    z_stat, p_val = z_test_two_proportions(hit_mis, hit_cov, bg_mis, bg_cov)

    # Compute “Estimated mutations per target copy (basepairs)” (float)
    length_of_target = len(hit_seq)
    true_diff_rate   = hit_rate - bg_rate
    est_mut_per_copy = max(true_diff_rate * length_of_target, 0.0)

    # Determine if ROI is a valid protein‐coding sequence
    is_protein = True
    reasons    = []
    seq_upper  = hit_seq.upper()
    if len(seq_upper) % 3 != 0:
        is_protein = False
        reasons.append(f"length {len(seq_upper)} is not a multiple of 3")
    if not seq_upper.startswith("ATG"):
        is_protein = False
        reasons.append("does not start with ATG")
    if seq_upper[-3:] not in {"TAA", "TAG", "TGA"}:
        is_protein = False
        reasons.append(f"does not end with a stop codon (found '{seq_upper[-3:]}')")

    # If protein, simulate AA distribution per copy using Poisson sampling
    if is_protein:
        logging.info(f"Simulating amino acid distribution with λ_bp={est_mut_per_copy:.2f}")
        aa_diffs = simulate_aa_distribution(est_mut_per_copy, hit_seq, n_trials=1000)
        avg_aa_mutations = sum(aa_diffs) / len(aa_diffs)
    else:
        aa_diffs = []
        avg_aa_mutations = None

    # ----------------------------
    # SAVE CSV FOR PANEL 1 DATA:
    # Mutation rates per position with background rate metadata
    # ----------------------------
    gene_mismatch_csv = os.path.join(RESULTS_DIR, "gene_mismatch_rate.csv")
    with open(gene_mismatch_csv, "w", newline="") as csvfile:
        # Write metadata as commented lines
        csvfile.write(f"# gene_id: {hit_id}\n")
        csvfile.write(f"# background_rate_per_base: {bg_rate:.6f}\n")
        csvfile.write(f"# background_rate_per_kb: {bg_rate_per_kb:.6f}\n")
        if is_protein:
            csvfile.write(f"# is_protein: True\n")
        else:
            csvfile.write(f"# is_protein: False\n")
            csvfile.write(f"# protein_invalid_reasons: {'; '.join(reasons)}\n")
        csvfile.write("position,mismatch_rate\n")
        for pos, rate in enumerate(hit_info["pos_rates"], start=1):
            csvfile.write(f"{pos},{rate:.6f}\n")
    logging.info(f"Saved CSV for gene mismatch rates: {gene_mismatch_csv}")

    # ----------------------------
    # ALIGN TO FULL PLASMID TO GET COVERAGE
    # ----------------------------
    sam_plasmid = run_minimap2(reads_pattern, plasmid_fasta, "plasmid_full_alignment")

    plasmid_cov = [0] * len(plasmid_seq)
    samfile_full = pysam.AlignmentFile(sam_plasmid, "r")
    for read in samfile_full.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if ref_pos is not None and 0 <= ref_pos < len(plasmid_seq):
                plasmid_cov[ref_pos] += 1
    samfile_full.close()

    # ----------------------------
    # SAVE CSV FOR PANEL 2 DATA:
    # Plasmid coverage with ROI metadata
    # ----------------------------
    plasmid_cov_csv = os.path.join(RESULTS_DIR, "plasmid_coverage.csv")
    with open(plasmid_cov_csv, "w", newline="") as csvfile:
        csvfile.write(f"# plasmid_id: {plasmid_id}\n")
        csvfile.write(f"# roi_gene_id: {hit_id}\n")
        csvfile.write(f"# roi_start (1-based): {idx+1}\n")
        csvfile.write(f"# roi_end (1-based): {idx + len(hit_seq)}\n")
        csvfile.write("position,coverage\n")
        for pos, cov in enumerate(plasmid_cov, start=1):
            csvfile.write(f"{pos},{cov}\n")
    logging.info(f"Saved CSV for plasmid coverage: {plasmid_cov_csv}")

    # ----------------------------
    # SAVE CSV FOR PANEL 3 DATA:
    # AA mutation distribution with metadata
    # ----------------------------
    aa_dist_csv = os.path.join(RESULTS_DIR, "aa_mutation_distribution.csv")
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
    # PREPARE PANEL FIGURE WITH 3 SUBPLOTS
    # ----------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    # axes[0]: Mutation rate over gene of interest
    # axes[1]: Coverage of plasmid with ROI shaded
    # axes[2]: KDE of AA mutations per copy

    # --- Panel 1: Mutation rate over gene of interest ---
    ax0 = axes[0]
    positions_gene = np.arange(1, len(hit_info["pos_rates"]) + 1)
    # Shade area under background rate
    ax0.axhspan(0, bg_rate, color='gray', alpha=0.3, label="Background rate")
    # Plot per-base mismatch rate for the hit
    ax0.plot(positions_gene, hit_info["pos_rates"],
             color="#2E8B57", linestyle='-', marker='o', markersize=4, alpha=0.8,
             label="Hit mismatch rate")
    ax0.set_title("Mismatch Rate per Position: Gene of Interest", fontsize=14, fontweight='bold')
    ax0.set_xlabel("Position (bp)", fontsize=12)
    ax0.set_ylabel("Mismatch Rate", fontsize=12)
    ax0.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
    ax0.tick_params(axis='both', which='minor', direction='in', length=3)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.legend(loc="upper right", frameon=False, fontsize=10)

    # --- Panel 2: Coverage of plasmid with ROI shaded ---
    ax1 = axes[1]
    plasmid_positions = np.arange(1, len(plasmid_cov) + 1)
    ax1.plot(plasmid_positions, plasmid_cov,
             linestyle='-', color='black', linewidth=1.0, alpha=0.8, label="Coverage")
    ax1.set_title("Full Plasmid Coverage with ROI Shaded", fontsize=14, fontweight='bold')
    ax1.set_xlabel("Position on Plasmid", fontsize=12)
    ax1.set_ylabel("Coverage (# reads)", fontsize=12)
    ax1.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
    ax1.tick_params(axis='both', which='minor', direction='in', length=3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    # Shade the ROI
    start_roi = idx + 1
    end_roi   = idx + len(hit_seq)
    ax1.axvspan(start_roi, end_roi, color='gray', alpha=0.3, label=f"ROI: {start_roi}–{end_roi}")
    #ax1.legend(loc="upper right", frameon=False, fontsize=10)

    # --- Panel 3: KDE of AA mutations per copy ---
    ax2 = axes[2]
    if is_protein and aa_diffs:
        x_vals = np.array(aa_diffs)
        # Use KDE if scipy is available and there are enough unique points
        if HAVE_SCIPY and len(np.unique(x_vals)) > 1:
            kde = gaussian_kde(x_vals)
            x_grid = np.linspace(0, max(x_vals), 200)
            ax2.plot(x_grid, kde(x_grid),
                     color="#C44E52", linewidth=2.0, alpha=0.8, label="KDE")
            ax2.fill_between(x_grid, kde(x_grid), color="#C44E52", alpha=0.3)
            ax2.set_ylim(bottom=0)
        else:
            # Fall back to histogram density
            ax2.hist(x_vals, bins=range(0, max(x_vals) + 2), density=True,
                     color="#C44E52", alpha=0.6, edgecolor='black', label="Histogram")
        ax2.set_title("Distribution of Amino Acid Mutations per gene", fontsize=14, fontweight='bold')
        ax2.set_xlabel("Number of AA Mutations", fontsize=12)
        ax2.set_ylabel("Density", fontsize=12)
        ax2.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
        ax2.tick_params(axis='both', which='minor', direction='in', length=3)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        #ax2.legend(loc="upper right", frameon=False, fontsize=10)
    else:
        ax2.text(0.5, 0.5, "Not a protein‐coding region", horizontalalignment='center',
                 verticalalignment='center', fontsize=12, color='gray', transform=ax2.transAxes)
        ax2.set_title("AA Mutation Distribution", fontsize=14, fontweight='bold')
        ax2.set_xlabel("Number of AA Mutations", fontsize=12)
        ax2.set_ylabel("Density", fontsize=12)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.set_xticks([])
        ax2.set_yticks([])

    # Save the combined figure as both PNG and PDF
    panel_path_png = os.path.join(RESULTS_DIR, "summary_panels.png")
    panel_path_pdf = os.path.join(RESULTS_DIR, "summary_panels.pdf")
    fig.savefig(panel_path_png, dpi=150, transparent=False)
    fig.savefig(panel_path_pdf)  # vector format
    plt.close(fig)
    logging.info(f"Saved combined panel figure as PNG: {panel_path_png}")
    logging.info(f"Saved combined panel figure as PDF: {panel_path_pdf}")

    # Clean up intermediate files
    if os.path.isdir(WORK_DIR):
        shutil.rmtree(WORK_DIR)

if __name__ == "__main__":
    main()
