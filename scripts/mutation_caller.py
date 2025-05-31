#!/usr/bin/env python3

import gzip
import glob
import re
import os
import tempfile
import subprocess
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline

from scipy.stats import fisher_exact, gaussian_kde

# ----------------------------------------------------
# Load flanking information
flanks_df = pd.read_csv("data/mutation_caller/mutation_caller.csv")
gene_start = flanks_df.loc[0, 'gene_flanks']
gene_end   = flanks_df.loc[1, 'gene_flanks']
gene_min   = int(flanks_df.loc[0, 'gene_min_max'])
gene_max   = int(flanks_df.loc[1, 'gene_min_max'])
pattern_gene = re.compile(
    rf'{gene_start}([ACGTNacgtn]{{{gene_min},{gene_max}}}){gene_end}',
    re.IGNORECASE
)

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def extract_gene(seq: str) -> str or None:
    """
    Given a raw nucleotide sequence string, search for the region between
    'gene_start' and 'gene_end'. If found in the forward strand, return
    the intervening sequence; otherwise, try the reverse complement. If
    still not found, return None.
    """
    match = pattern_gene.search(seq)
    if match:
        return match.group(1)
    match = pattern_gene.search(reverse_complement(seq))
    if match:
        return match.group(1)
    return None

def process_fastq(file_path: str) -> dict:
    """
    Read a single .fastq.gz, extract all reads whose sequence contains
    the gene region (via extract_gene), and return a dict mapping
    read_id -> gene_sequence.
    """
    gene_reads = {}
    with gzip.open(file_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # skip '+'
            f.readline()  # skip quality
            gene = extract_gene(seq)
            if gene and gene_min <= len(gene) <= gene_max:
                read_id = header.strip()[1:]  # drop '@'
                gene_reads[read_id] = gene
    return gene_reads

def align_to_reference(gene_seqs: dict, reference: str) -> (str, dict):
    """
    Align the reference sequence plus every gene_seq (dict: read_id -> seq)
    using MAFFT. Returns (aligned_reference_seq, aligned_reads_dict).
    'aligned_reads_dict' maps read_id -> aligned_sequence (with gaps).
    """
    # 1) Write reference + reads into a temp FASTA
    with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="tmp_align_") as tmp_fasta:
        ref_record = SeqRecord(Seq(reference), id="REF", description="")
        SeqIO.write(ref_record, tmp_fasta, "fasta")
        for rid, seq in gene_seqs.items():
            SeqIO.write(SeqRecord(Seq(seq), id=rid, description=""), tmp_fasta, "fasta")
        input_fasta = tmp_fasta.name

    # 2) Run MAFFT (capture stdout)
    mafft_cline = MafftCommandline(input=input_fasta)
    process = subprocess.Popen(
        str(mafft_cline),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    stdout, stderr = process.communicate()

    # 3) Write MAFFT's alignment output to another temp file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="tmp_align_out_") as tmp_out:
        tmp_out.write(stdout)
        output_fasta = tmp_out.name

    # 4) Parse the aligned FASTA
    alignments = AlignIO.read(output_fasta, "fasta")
    aligned_reads = {}
    aligned_ref = None
    for record in alignments:
        if record.id == "REF":
            aligned_ref = str(record.seq)
        else:
            aligned_reads[record.id] = str(record.seq)

    # 5) Clean up temp files
    os.remove(input_fasta)
    os.remove(output_fasta)

    return aligned_ref, aligned_reads

def identify_substitutions(ref: str, aligned_reads: dict) -> dict:
    """
    Given:
      - ref: aligned reference sequence (with gaps)
      - aligned_reads: dict mapping read_id -> aligned_read_seq (with gaps)
    Return:
      subs_by_read: a dict mapping read_id -> list of strings like "a76c (S26R)".
    """
    subs_by_read = defaultdict(list)

    # Build mapping: alignment index -> ungapped reference index
    alignment_to_ref_idx = {}
    ref_seq_clean_chars = []
    ref_nt_index = 0
    for ai, base in enumerate(ref):
        if base != "-":
            alignment_to_ref_idx[ai] = ref_nt_index
            ref_seq_clean_chars.append(base)
            ref_nt_index += 1
        else:
            alignment_to_ref_idx[ai] = None
    ref_seq_clean = "".join(ref_seq_clean_chars)

    # For each read, compare aligned positions
    for rid, aln_seq in aligned_reads.items():
        for ai, query_base in enumerate(aln_seq):
            ref_idx = alignment_to_ref_idx.get(ai)
            if ref_idx is None:
                continue  # gap in ref
            ref_base = ref_seq_clean[ref_idx]
            if query_base == "-" or query_base == ref_base:
                continue  # no substitution
            # Nucleotide-level mutation (1-based indexing)
            nt_mut = f"{ref_base}{ref_idx+1}{query_base}"
            # Codon index (0-based)
            codon_index = ref_idx // 3
            codon_start = codon_index * 3
            codon_end = codon_start + 3
            if codon_end > len(ref_seq_clean):
                continue  # partial codon at end
            ref_codon = ref_seq_clean[codon_start:codon_end]
            if len(ref_codon) != 3:
                continue
            # Build the read's codon
            rc_list = list(ref_codon)
            pos_in_codon = ref_idx % 3
            rc_list[pos_in_codon] = query_base
            try:
                aa_from = str(Seq(ref_codon).translate())
                aa_to   = str(Seq("".join(rc_list)).translate())
            except Exception:
                aa_from, aa_to = "?", "?"
            aa_mut = f"{aa_from}{codon_index+1}{aa_to}"
            subs_by_read[rid].append(f"{nt_mut} ({aa_mut})")

    return subs_by_read

def find_cooccurring_aa(subs_by_read_aa: dict, frequent_aa: set, output_dir: str):
    """
    Build co-occurrence on amino-acid–level substitutions.
    subs_by_read_aa: { read_id: [ "S26R", "L35V", ... ], ... }
    frequent_aa: set of those AA changes that pass the threshold
    output_dir: path to write two CSVs:
      - cooccurring_AA_baseline.csv
      - cooccurring_AA_fisher.csv
    """
    # 1) Sort frequent_aa to fix column order
    aa_list = sorted(frequent_aa)
    aa_idx = { aa: i for i, aa in enumerate(aa_list) }

    # 2) Build a binary matrix: rows=reads, cols=frequent AA changes
    matrix = []
    for aa_calls in subs_by_read_aa.values():
        row = [0] * len(aa_list)
        any_frequent = False
        for aa in aa_calls:
            if aa in aa_idx:
                row[aa_idx[aa]] = 1
                any_frequent = True
        if any_frequent:
            matrix.append(row)
    if not matrix:
        # No reads with frequent AA changes → write empty CSVs
        pd.DataFrame(columns=["AA1","AA2","Both_Count","AA1_Count","AA2_Count"])\
          .to_csv(os.path.join(output_dir, "cooccurring_AA_baseline.csv"), index=False)
        pd.DataFrame(columns=["AA1","AA2","p-value"])\
          .to_csv(os.path.join(output_dir, "cooccurring_AA_fisher.csv"), index=False)
        return

    df = pd.DataFrame(matrix, columns=aa_list)

    # 3) Baseline co-occurrence: both >=2 OR perfect co-presence
    cooccur_simple = []
    for i in range(len(aa_list)):
        for j in range(i+1, len(aa_list)):
            A = df.iloc[:, i]
            B = df.iloc[:, j]
            both     = int(((A == 1) & (B == 1)).sum())
            A_total  = int((A == 1).sum())
            B_total  = int((B == 1).sum())
            if (both >= 2) or (both > 0 and both == A_total == B_total):
                cooccur_simple.append((aa_list[i], aa_list[j], both, A_total, B_total))
    simple_df = pd.DataFrame(
        cooccur_simple,
        columns=["AA1", "AA2", "Both_Count", "AA1_Count", "AA2_Count"]
    )
    simple_df.to_csv(
        os.path.join(output_dir, "cooccurring_AA_baseline.csv"),
        index=False
    )

    # 4) Fisher's exact test for p < 0.05
    cooccur_fisher = []
    for i in range(len(aa_list)):
        for j in range(i+1, len(aa_list)):
            A = df.iloc[:, i]
            B = df.iloc[:, j]
            both    = int(((A == 1) & (B == 1)).sum())
            only_A  = int(((A == 1) & (B == 0)).sum())
            only_B  = int(((A == 0) & (B == 1)).sum())
            neither = int(((A == 0) & (B == 0)).sum())
            table = [[both, only_A], [only_B, neither]]
            try:
                _, p = fisher_exact(table)
                if p < 0.05:
                    cooccur_fisher.append((aa_list[i], aa_list[j], p))
            except Exception:
                continue
    fisher_df = pd.DataFrame(
        cooccur_fisher,
        columns=["AA1", "AA2", "p-value"]
    )
    fisher_df.to_csv(
        os.path.join(output_dir, "cooccurring_AA_fisher.csv"),
        index=False
    )

def main():
    # 0) Load the reference template (once)
    template_path = "data/mutation_caller/mutation_caller_template.fasta"
    template = str(next(SeqIO.parse(template_path, "fasta")).seq)

    # 1) Find all .fastq.gz files
    fastq_files = glob.glob("data/mutation_caller/*.fastq.gz")
    if not fastq_files:
        print("No .fastq.gz files found in data/mutation_caller/. Exiting.")
        return

    # 2) Process each FASTQ file individually
    for fastq_path in fastq_files:
        base_name = os.path.basename(fastq_path)    # e.g. "SAMPLE1.fastq.gz"
        sample_name = os.path.splitext(base_name)[0]
        if sample_name.endswith(".fastq"):
            sample_name = sample_name[:-6]  # strip ".fastq"

        out_dir = os.path.join("results", "mutation_caller", base_name)
        os.makedirs(out_dir, exist_ok=True)

        print(f"\nProcessing {fastq_path} …")
        reads_this_file = process_fastq(fastq_path)
        print(f"  → {len(reads_this_file)} valid gene reads found in {base_name}")

        if not reads_this_file:
            print("  → No valid reads; skipping this sample.")
            continue

        # 3) Align these reads to the reference
        ref_aligned, aligned_reads = align_to_reference(reads_this_file, template)

        # 4) Identify nucleotide + AA substitutions per read
        subs_by_read = identify_substitutions(ref_aligned, aligned_reads)

        # 5) Build an AA-only dictionary: read_id -> [ "S26R", "L35V", ... ]
        subs_by_read_aa = {}
        for rid, nt_aa_list in subs_by_read.items():
            aa_list = []
            for nt_aa in nt_aa_list:
                parts = nt_aa.split()
                if len(parts) >= 2 and parts[1].startswith("(") and parts[1].endswith(")"):
                    aa_change = parts[1][1:-1]  # strip parentheses
                    aa_list.append(aa_change)
            if aa_list:
                subs_by_read_aa[rid] = aa_list

        # 6) Count all AA substitutions, plot bar + KDE (with log x-axis)
        all_aa_counts = Counter(
            aa for aa_list in subs_by_read_aa.values() for aa in aa_list
        )
        if not all_aa_counts:
            print("  → No amino-acid substitutions found; skipping plotting/co-occurrence.")
            continue

        # Convert counts to arrays for KDE
        counts_array = np.array(list(all_aa_counts.values()), dtype=float)

        # Set up a figure with two side-by-side subplots
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # 6a) Bar plot of AA substitution frequencies
        ax_bar = axes[0]
        substitutions = list(all_aa_counts.keys())
        frequencies   = list(all_aa_counts.values())
        ax_bar.bar(substitutions, frequencies)
        ax_bar.set_xticklabels(substitutions, rotation=90, fontsize=8)
        ax_bar.set_ylabel("Count")
        ax_bar.set_title("Amino-Acid Substitution Frequencies")

        # 6b) KDE plot of the frequency distribution (log-scale x-axis)
        ax_kde = axes[1]
        # Compute Gaussian KDE on raw counts
        try:
            kde = gaussian_kde(counts_array)
            x_min = max(1, counts_array.min())
            x_max = counts_array.max()
            # Generate points evenly spaced on log scale
            x_vals = np.logspace(np.log10(x_min), np.log10(x_max), 200)
            y_vals = kde(x_vals)
            ax_kde.plot(x_vals, y_vals, linewidth=2)
        except Exception:
            # Fall back: if KDE fails (e.g., only one unique count), plot a histogram instead
            ax_kde.hist(counts_array, bins='auto', density=True, alpha=0.6)
        ax_kde.set_xscale("log")
        ax_kde.set_xlabel("Substitution Count (log scale)")
        ax_kde.set_ylabel("Density")
        ax_kde.set_title("KDE of AA Substitution Frequencies")

        plt.tight_layout()
        hist_path = os.path.join(out_dir, "aa_substitution_frequency.png")
        plt.savefig(hist_path)
        print(f"  → Saved combined bar + KDE plot to {hist_path}")

        # Show the plot so user can choose a threshold
        plt.show()
        plt.close()

        # 7) Ask the user for an amino-acid–level threshold
        threshold_aa = int(input(f"Enter amino-acid substitution threshold for {base_name}: "))

        # 8) Determine which AA changes meet the threshold
        frequent_aa = { aa for aa, cnt in all_aa_counts.items() if cnt >= threshold_aa }
        print(f"  → {len(frequent_aa)} amino-acid substitutions meet threshold ≥ {threshold_aa}")

        # 9) Run AA-level co-occurrence and save results
        find_cooccurring_aa(subs_by_read_aa, frequent_aa, out_dir)
        print(f"  → AA co-occurrence CSVs written into {out_dir}/")

    print("\nAll FASTQ files processed.")

if __name__ == "__main__":
    main()
