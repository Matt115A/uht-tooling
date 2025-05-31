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

# Use a built-in Matplotlib style ("ggplot") for consistency
plt.style.use("ggplot")

# Try to import scipy for Z‐test p‐value and KDE
try:
    from scipy.stats import norm, gaussian_kde
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

    # Remove hit region from plasmid → “background plasmid”
    idx = plasmid_seq.find(hit_seq)
    if idx == -1:
        logging.error("Gene region not found in plasmid")
        return
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
    chunks_fasta = create_multi_fasta(chunks, work_dir)
    sam_chunks   = run_minimap2(fastq_path, chunks_fasta, "plasmid_chunks_alignment", work_dir)

    # Align to hit (target) alone
    sam_hit = run_minimap2(fastq_path, ref_hit_fasta, "hit_alignment", work_dir)

    # Compute mismatch stats for background chunks
    chunk_refs      = { f"chunk_{i+1}": seq for i, seq in enumerate(chunks) }
    mismatch_chunks = compute_mismatch_stats_sam(sam_chunks, chunk_refs)

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

    # Compute “Estimated mutations per target copy (basepairs)” (float)
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
        aa_diffs         = simulate_aa_distribution(est_mut_per_copy, hit_seq, n_trials=1000)
        avg_aa_mutations = sum(aa_diffs) / len(aa_diffs)
    else:
        aa_diffs         = []
        avg_aa_mutations = None

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
    # SAVE CSV FOR POSITIONS ABOVE BACKGROUND: BASE DISTRIBUTION
    # ----------------------------
    base_dist_csv = os.path.join(results_dir, "gene_base_distribution.csv")
    with open(base_dist_csv, "w", newline="") as csvfile:
        csvfile.write(f"# gene_id: {hit_id}\n")
        csvfile.write(f"# background_rate_per_base: {bg_rate:.6f}\n")
        csvfile.write("position,A_count,C_count,G_count,T_count,N_count,coverage\n")
        for pos0, rate in enumerate(hit_info["pos_rates"]):
            if rate > bg_rate:
                counts   = base_counts[pos0]
                coverage = hit_info["cov"][pos0]
                csvfile.write(
                    f"{pos0+1},"
                    f"{counts['A']},{counts['C']},{counts['G']},{counts['T']},{counts['N']},"
                    f"{coverage}\n"
                )
    logging.info(f"Saved CSV for base distribution at above-background positions: {base_dist_csv}")

    # ----------------------------
    # IF PROTEIN: SAVE CSV FOR AA SUBSTITUTION RATES AT ABOVE-BACKGROUND POSITIONS
    # ----------------------------
    if is_protein:
        aa_subst_csv = os.path.join(results_dir, "gene_aa_substitution_rates.csv")
        aa_rates     = {}
        for pos0, rate in enumerate(hit_info["pos_rates"]):
            if rate <= bg_rate:
                continue
            ref_base = seq_upper[pos0]
            cov      = hit_info["cov"][pos0]
            if cov == 0:
                continue
            counts = base_counts[pos0]
            codon_start         = (pos0 // 3) * 3
            codon_pos_in_codon  = pos0 % 3
            codon_seq           = seq_upper[codon_start : codon_start + 3]
            ref_aa              = str(Seq(codon_seq).translate(to_stop=False))
            aa_position         = (codon_start // 3) + 1  # 1-based AA index

            for base_allele, cnt in counts.items():
                if base_allele == ref_base or cnt == 0 or base_allele == "N":
                    continue
                mutated_codon_list = list(codon_seq)
                mutated_codon_list[codon_pos_in_codon] = base_allele
                mutated_codon = "".join(mutated_codon_list)
                mutated_aa    = str(Seq(mutated_codon).translate(to_stop=False))
                if mutated_aa == ref_aa:
                    continue
                rate_contribution = cnt / cov
                key = (pos0 + 1, aa_position, ref_aa, mutated_aa)
                aa_rates[key] = aa_rates.get(key, 0.0) + rate_contribution

        with open(aa_subst_csv, "w", newline="") as csvfile:
            csvfile.write(f"# gene_id: {hit_id}\n")
            csvfile.write(f"# background_rate_per_base: {bg_rate:.6f}\n")
            csvfile.write("nt_position,aa_position,ref_aa,alt_aa,rate\n")
            for (nt_pos, aa_pos, ref_aa, alt_aa), rate_val in sorted(aa_rates.items()):
                csvfile.write(f"{nt_pos},{aa_pos},{ref_aa},{alt_aa},{rate_val:.6f}\n")
        logging.info(f"Saved CSV for amino acid substitution rates: {aa_subst_csv}")

    # ----------------------------
    # SAVE CSV FOR MUTATION RATES (PANEL 1)
    # ----------------------------
    gene_mismatch_csv = os.path.join(results_dir, "gene_mismatch_rate.csv")
    with open(gene_mismatch_csv, "w", newline="") as csvfile:
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
    sam_plasmid = run_minimap2(fastq_path, plasmid_fasta, "plasmid_full_alignment", work_dir)

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
    # SAVE CSV FOR PLASMID COVERAGE (PANEL 2)
    # ----------------------------
    plasmid_cov_csv = os.path.join(results_dir, "plasmid_coverage.csv")
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
    # SAVE CSV FOR AA MUTATION DISTRIBUTION (PANEL 3)
    # ----------------------------
    aa_dist_csv = os.path.join(results_dir, "aa_mutation_distribution.csv")
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
    ax0.axhspan(0, bg_rate, color='gray', alpha=0.3, label="Background rate")
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
    start_roi = idx + 1
    end_roi   = idx + len(hit_seq)
    ax1.axvspan(start_roi, end_roi, color='gray', alpha=0.3, label=f"ROI: {start_roi}–{end_roi}")

    # --- Panel 3: KDE of AA mutations per copy ---
    ax2 = axes[2]
    if is_protein and aa_diffs:
        x_vals = np.array(aa_diffs)
        if HAVE_SCIPY and len(np.unique(x_vals)) > 1:
            kde = gaussian_kde(x_vals)
            x_grid = np.linspace(0, max(x_vals), 200)
            ax2.plot(x_grid, kde(x_grid),
                     color="#C44E52", linewidth=2.0, alpha=0.8, label="KDE")
            ax2.fill_between(x_grid, kde(x_grid), color="#C44E52", alpha=0.3)
            ax2.set_ylim(bottom=0)
        else:
            ax2.hist(x_vals, bins=range(0, max(x_vals) + 2), density=True,
                     color="#C44E52", alpha=0.6, edgecolor='black', label="Histogram")
        ax2.set_title("Distribution of Amino Acid Mutations per Gene", fontsize=14, fontweight='bold')
        ax2.set_xlabel("Number of AA Mutations", fontsize=12)
        ax2.set_ylabel("Density", fontsize=12)
        ax2.tick_params(axis='both', which='major', labelsize=10, direction='in', length=6)
        ax2.tick_params(axis='both', which='minor', direction='in', length=3)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
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
    panel_path_png = os.path.join(results_dir, "summary_panels.png")
    panel_path_pdf = os.path.join(results_dir, "summary_panels.pdf")
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
    pdf_path = os.path.join(results_dir, f"{sample_name}_mutation_spectrum.pdf")
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
    sample_summary_path = os.path.join(results_dir, "summary.txt")
    with open(sample_summary_path, "w") as txtf:
        txtf.write(f"Sample: {sample_name}\n")
        txtf.write(f"{'=' * (8 + len(sample_name))}\n\n")
        txtf.write("1) Background (all chunks):\n")
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

    # ----------------------------
    # APPEND TO MASTER SUMMARY TXT
    # ----------------------------
    # Each sample's key stats go into a single line in master_summary.txt
    with open(master_summary_path, "a") as masterf:
        masterf.write(
            f"{sample_name}\t"
            f"{bg_rate:.6e}\t"
            f"{hit_rate:.6e}\t"
            f"{z_stat:.4f}\t"
            f"{p_val if p_val is not None else 'NA'}\t"
            f"{est_mut_per_copy:.6e}\t"
            f"{is_protein}\n"
        )

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
