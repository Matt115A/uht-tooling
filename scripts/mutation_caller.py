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
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def extract_gene(seq: str) -> str or None:
    match = pattern_gene.search(seq)
    if match:
        return match.group(1)
    match = pattern_gene.search(reverse_complement(seq))
    if match:
        return match.group(1)
    return None

# ----------------------------------------------------
def process_fastq(file_path: str) -> dict:
    gene_reads = {}
    with gzip.open(file_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # '+'
            f.readline()  # quality
            gene = extract_gene(seq)
            if gene and gene_min <= len(gene) <= gene_max:
                read_id = header.strip()[1:]
                gene_reads[read_id] = gene
    return gene_reads

# ----------------------------------------------------
def align_to_reference(gene_seqs: dict, reference: str) -> (str, dict):
    print(f"DEBUG: Aligning {len(gene_seqs)} reads to reference of length {len(reference)}")  # debug
    with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="tmp_align_") as tmp:
        SeqIO.write(SeqRecord(Seq(reference), id="REF", description=""), tmp, "fasta")
        for rid, seq in gene_seqs.items():
            SeqIO.write(SeqRecord(Seq(seq), id=rid, description=""), tmp, "fasta")
        inp = tmp.name

    mafft = MafftCommandline(input=inp)
    proc = subprocess.Popen(str(mafft), shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    stdout, stderr = proc.communicate()
    print(f"DEBUG: MAFFT stderr:\n{stderr}")  # debug output

    with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="tmp_align_out_") as tmp2:
        tmp2.write(stdout)
        outp = tmp2.name

    aln = AlignIO.read(outp, "fasta")
    aligned_ref = None
    aligned_reads = {}
    for rec in aln:
        if rec.id == "REF":
            aligned_ref = str(rec.seq)
            print(f"DEBUG: Aligned REF (first 50 bp): {aligned_ref[:50]}")
        else:
            aligned_reads[rec.id] = str(rec.seq)
    print(f"DEBUG: Sample aligned read lengths (first 5): {[len(s) for s in list(aligned_reads.values())[:5]]}")

    os.remove(inp)
    os.remove(outp)
    return aligned_ref, aligned_reads

# ----------------------------------------------------
def identify_substitutions(ref: str, aligned_reads: dict) -> dict:
    """
    Identify amino-acid substitutions by collapsing per-base mismatches within each codon
    into a single codon-level mutation.
    """
    subs_by_read = defaultdict(list)

    # build ungapped reference and mapping from alignment idx to ref idx
    aln2ref = {}
    ref_clean = []
    ri = 0
    for ai, b in enumerate(ref):
        if b != '-':
            aln2ref[ai] = ri
            ref_clean.append(b)
            ri += 1
        else:
            aln2ref[ai] = None
    ref_clean = ''.join(ref_clean)

    # invert mapping: ref idx -> alignment idx (first occurrence)
    ref2aln = {}
    for ai, r_idx in aln2ref.items():
        if r_idx is not None and r_idx not in ref2aln:
            ref2aln[r_idx] = ai

    codon_count = len(ref_clean) // 3
    for rid, seq in aligned_reads.items():
        # for each codon in the reference
        for codon_i in range(codon_count):
            start_r = codon_i * 3
            end_r   = start_r + 3
            codon_ref = ref_clean[start_r:end_r]
            codon_read = []
            valid = True
            diff = False
            # collect read bases for each codon position
            for offset in range(3):
                r_idx = start_r + offset
                ai = ref2aln.get(r_idx)
                if ai is None:
                    valid = False
                    break
                base_r = codon_ref[offset]
                base_q = seq[ai]
                if base_q == '-' or base_q is None:
                    valid = False
                    break
                codon_read.append(base_q)
                if base_q != base_r:
                    diff = True
            if not valid or not diff:
                continue
            # translate codons
            try:
                aa_from = str(Seq(codon_ref).translate())
                aa_to   = str(Seq(''.join(codon_read)).translate())
            except Exception:
                aa_from, aa_to = '?', '?'
            aa_mut = f"{aa_from}{codon_i+1}{aa_to}"
            nt_mut = f"{codon_ref}->{''.join(codon_read)}"
            subs_by_read[rid].append(f"{nt_mut} ({aa_mut})")

    return subs_by_read

# ----------------------------------------------------
def find_cooccurring_aa(subs_by_read_aa: dict, frequent_aa: set, output_dir: str, sample_name: str):
    aa_list = sorted(frequent_aa)
    aa_idx  = {aa:i for i,aa in enumerate(aa_list)}

    matrix = []
    for calls in subs_by_read_aa.values():
        row = [0]*len(aa_list)
        any_f = False
        for aa in calls:
            if aa in aa_idx:
                row[aa_idx[aa]] = 1
                any_f = True
        if any_f:
            matrix.append(row)

    base1 = f"{sample_name}_cooccurring_AA_baseline.csv"
    base2 = f"{sample_name}_cooccurring_AA_fisher.csv"

    if not matrix:
        pd.DataFrame(columns=["AA1","AA2","Both_Count","AA1_Count","AA2_Count"])\
          .to_csv(os.path.join(output_dir, base1), index=False)
        pd.DataFrame(columns=["AA1","AA2","p-value"])\
          .to_csv(os.path.join(output_dir, base2), index=False)
        return base1, base2

    df = pd.DataFrame(matrix, columns=aa_list)
    simple, fisher = [], []

    for i in range(len(aa_list)):
        for j in range(i+1, len(aa_list)):
            A, B = df.iloc[:,i], df.iloc[:,j]
            both   = int(((A==1)&(B==1)).sum())
            a_tot  = int((A==1).sum())
            b_tot  = int((B==1).sum())
            if (both>=2) or (both>0 and both==a_tot==b_tot):
                simple.append((aa_list[i], aa_list[j], both, a_tot, b_tot))
            table = [[both, int(((A==1)&(B==0)).sum())],
                     [int(((A==0)&(B==1)).sum()), int(((A==0)&(B==0)).sum())]]
            try:
                _, p = fisher_exact(table)
                if p < 0.05:
                    fisher.append((aa_list[i], aa_list[j], p))
            except:
                pass

    pd.DataFrame(simple, columns=["AA1","AA2","Both_Count","AA1_Count","AA2_Count"])\
      .to_csv(os.path.join(output_dir, base1), index=False)
    pd.DataFrame(fisher, columns=["AA1","AA2","p-value"])\
      .to_csv(os.path.join(output_dir, base2), index=False)
    return base1, base2

# ----------------------------------------------------
def main():
    # load template (which is the ORF only)
    tmpl = "data/mutation_caller/mutation_caller_template.fasta"
    full_ref_record = next(SeqIO.parse(tmpl, "fasta"))
    full_ref        = str(full_ref_record.seq)

    # if primers flank the ORF on the template, trim; otherwise use full ORF
    if full_ref.startswith(gene_start) and full_ref.endswith(gene_end):
        reference = full_ref[len(gene_start):-len(gene_end)]
    else:
        reference = full_ref

    fastqs = glob.glob("data/mutation_caller/*.fastq.gz")
    if not fastqs:
        print("No FASTQ found; exiting.")
        return

    for fq in fastqs:
        base   = os.path.basename(fq)
        sample = os.path.splitext(base)[0]
        if sample.endswith(".fastq"):
            sample = sample[:-6]

        out_dir = os.path.join("results","mutation_caller", base)
        os.makedirs(out_dir, exist_ok=True)

        print(f"\nProcessing {base} …")
        gene_reads = process_fastq(fq)
        n_reads    = len(gene_reads)
        print(f"  → {n_reads} valid gene reads")

        if n_reads == 0:
            continue

        aligned_ref, aligned_reads = align_to_reference(gene_reads, reference)
        subs = identify_substitutions(aligned_ref, aligned_reads)

        # AA-only dict
        subs_aa = {rid: [it.split()[1][1:-1]
                         for it in items if "(" in it and it.endswith(")")]
                   for rid, items in subs.items()}

        counts = Counter(aa for aas in subs_aa.values() for aa in aas)
        if not counts:
            print("  → No AA substitutions; skipping.")
            continue

        # plot frequencies with proper tick alignment
        keys = list(counts.keys())
        values = [counts[k] for k in keys]
        idx = np.arange(len(keys))
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,6))
        ax1.bar(idx, values)
        ax1.set_xticks(idx)
        ax1.set_xticklabels(keys, rotation=90, fontsize=8)
        ax1.set_ylabel("Count")
        ax1.set_title("Amino-Acid Substitution Frequencies")

        # density/KDE
        arr = np.array(values, dtype=float)
        try:
            kde = gaussian_kde(arr)
            xmin, xmax = max(1, arr.min()), arr.max()
            xs = np.logspace(np.log10(xmin), np.log10(xmax), 200)
            ax2.plot(xs, kde(xs), linewidth=2)
        except:
            ax2.hist(arr, bins="auto", density=True, alpha=0.6)
        ax2.set_xscale("log")
        ax2.set_xlabel("Substitution Count (log scale)")
        ax2.set_ylabel("Density")
        ax2.set_title("KDE of AA Substitution Frequencies")

        plt.tight_layout()
        hist_png = os.path.join(out_dir, f"{sample}_aa_substitution_frequency.png")
        plt.savefig(hist_png)
        plt.show()
        plt.close()
        print(f"  → Saved plot to {hist_png}")

        # prompt threshold after plot
        thr = int(input(f"Threshold for {base}: "))
        frequent = {aa for aa,c in counts.items() if c >= thr}
        print(f"  → {len(frequent)} substitutions ≥ {thr}")

        freq_csv = os.path.join(out_dir, f"{sample}_frequent_aa_counts.csv")
        pd.DataFrame(
            [(aa, counts[aa]) for aa in sorted(frequent)],
            columns=["AA","Count"]
        ).to_csv(freq_csv, index=False)
        print(f"  → Saved frequent counts to {freq_csv}")

        co_base, co_fish = find_cooccurring_aa(subs_aa, frequent, out_dir, sample)
        print(f"  → Co-occurrence files: {co_base}, {co_fish}")

        report_path = os.path.join(out_dir, f"{sample}_report.txt")
        with open(report_path, "w") as rpt:
            rpt.write(f"Sample: {sample}\n")
            rpt.write(f"Valid gene reads: {n_reads}\n")
            rpt.write(f"Unique AA substitutions: {len(counts)}\n")
            rpt.write(f"User threshold: {thr}\n")
            rpt.write(f"Frequent AA substitutions (≥ {thr}): {len(frequent)}\n\n")
            rpt.write("Frequent AA counts:\n")
            rpt.write("AA\tCount\n")
            for aa in sorted(frequent):
                rpt.write(f"{aa}\t{counts[aa]}\n")
            rpt.write("\nGenerated files:\n")
            rpt.write(f"- Plot: {os.path.basename(hist_png)}\n")
            rpt.write(f"- Frequent counts: {os.path.basename(freq_csv)}\n")
            rpt.write(f"- Co-occurrence baseline: {co_base}\n")
            rpt.write(f"- Co-occurrence fisher:   {co_fish}\n")
        print(f"  → Report written to {report_path}")

    print("\nAll FASTQ files processed.")

if __name__ == "__main__":
    main()
