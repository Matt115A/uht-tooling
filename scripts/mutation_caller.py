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

def align_to_reference(gene_seqs: dict, reference: str) -> (str, dict):
    with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="tmp_align_") as tmp:
        SeqIO.write(SeqRecord(Seq(reference), id="REF", description=""), tmp, "fasta")
        for rid, seq in gene_seqs.items():
            SeqIO.write(SeqRecord(Seq(seq), id=rid, description=""), tmp, "fasta")
        inp = tmp.name

    mafft = MafftCommandline(input=inp)
    proc = subprocess.Popen(str(mafft), shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    stdout, _ = proc.communicate()

    with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix="tmp_align_out_") as tmp2:
        tmp2.write(stdout)
        outp = tmp2.name

    aln = AlignIO.read(outp, "fasta")
    aligned_ref = None
    aligned_reads = {}
    for rec in aln:
        if rec.id == "REF":
            aligned_ref = str(rec.seq)
        else:
            aligned_reads[rec.id] = str(rec.seq)

    os.remove(inp)
    os.remove(outp)
    return aligned_ref, aligned_reads

def identify_substitutions(ref: str, aligned_reads: dict) -> dict:
    subs_by_read = defaultdict(list)

    # map alignment idx -> ungapped ref idx
    aln2ref = {}
    ref_clean = []
    ri = 0
    for ai, b in enumerate(ref):
        if b != "-":
            aln2ref[ai] = ri
            ref_clean.append(b)
            ri += 1
        else:
            aln2ref[ai] = None
    ref_clean = "".join(ref_clean)

    for rid, seq in aligned_reads.items():
        for ai, qb in enumerate(seq):
            r_idx = aln2ref.get(ai)
            if r_idx is None:
                continue
            rb = ref_clean[r_idx]
            if qb == "-" or qb == rb:
                continue
            nt_mut = f"{rb}{r_idx+1}{qb}"
            codon_i = r_idx // 3
            start = codon_i*3
            end   = start+3
            if end > len(ref_clean):
                continue
            codon = ref_clean[start:end]
            if len(codon) != 3:
                continue
            rc = list(codon)
            rc[r_idx%3] = qb
            try:
                aa_from = str(Seq(codon).translate())
                aa_to   = str(Seq("".join(rc)).translate())
            except:
                aa_from, aa_to = "?", "?"
            aa_mut = f"{aa_from}{codon_i+1}{aa_to}"
            subs_by_read[rid].append(f"{nt_mut} ({aa_mut})")

    return subs_by_read

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

    # prepare filenames
    base1 = f"{sample_name}_cooccurring_AA_baseline.csv"
    base2 = f"{sample_name}_cooccurring_AA_fisher.csv"

    if not matrix:
        pd.DataFrame(columns=["AA1","AA2","Both_Count","AA1_Count","AA2_Count"])\
          .to_csv(os.path.join(output_dir, base1), index=False)
        pd.DataFrame(columns=["AA1","AA2","p-value"])\
          .to_csv(os.path.join(output_dir, base2), index=False)
        return base1, base2

    df = pd.DataFrame(matrix, columns=aa_list)

    # baseline
    simple = []
    for i in range(len(aa_list)):
        for j in range(i+1, len(aa_list)):
            A, B = df.iloc[:,i], df.iloc[:,j]
            both   = int(((A==1)&(B==1)).sum())
            a_tot  = int((A==1).sum())
            b_tot  = int((B==1).sum())
            if (both>=2) or (both>0 and both==a_tot==b_tot):
                simple.append((aa_list[i], aa_list[j], both, a_tot, b_tot))
    pd.DataFrame(simple, columns=["AA1","AA2","Both_Count","AA1_Count","AA2_Count"])\
      .to_csv(os.path.join(output_dir, base1), index=False)

    # fisher
    fisher = []
    for i in range(len(aa_list)):
        for j in range(i+1, len(aa_list)):
            A, B = df.iloc[:,i], df.iloc[:,j]
            both   = int(((A==1)&(B==1)).sum())
            onlyA  = int(((A==1)&(B==0)).sum())
            onlyB  = int(((A==0)&(B==1)).sum())
            neither= int(((A==0)&(B==0)).sum())
            table = [[both, onlyA], [onlyB, neither]]
            try:
                _, p = fisher_exact(table)
                if p < 0.05:
                    fisher.append((aa_list[i], aa_list[j], p))
            except:
                pass
    pd.DataFrame(fisher, columns=["AA1","AA2","p-value"])\
      .to_csv(os.path.join(output_dir, base2), index=False)

    return base1, base2

def main():
    # load template
    tmpl = "data/mutation_caller/mutation_caller_template.fasta"
    reference = str(next(SeqIO.parse(tmpl, "fasta")).seq)

    fastqs = glob.glob("data/mutation_caller/*.fastq.gz")
    if not fastqs:
        print("No FASTQ found; exiting.")
        return

    for fq in fastqs:
        base = os.path.basename(fq)             # e.g. "S1.fastq.gz"
        sample = os.path.splitext(base)[0]      # "S1"
        if sample.endswith(".fastq"):
            sample = sample[:-6]

        out_dir = os.path.join("results","mutation_caller", base)
        os.makedirs(out_dir, exist_ok=True)

        print(f"\nProcessing {base} …")
        gene_reads = process_fastq(fq)
        n_reads = len(gene_reads)
        print(f"  → {n_reads} valid gene reads")

        if n_reads == 0:
            continue

        aligned_ref, aligned_reads = align_to_reference(gene_reads, reference)
        subs = identify_substitutions(aligned_ref, aligned_reads)

        # AA-only dict
        subs_aa = {}
        for rid, items in subs.items():
            aas = [it.split()[1][1:-1]
                   for it in items
                   if "(" in it and it.endswith(")")]
            if aas:
                subs_aa[rid] = aas

        counts = Counter(aa for aas in subs_aa.values() for aa in aas)
        if not counts:
            print("  → No AA substitutions; skipping.")
            continue

        # plot
        arr = np.array(list(counts.values()), dtype=float)
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,6))
        ax1.bar(list(counts.keys()), list(counts.values()))
        ax1.set_xticklabels(list(counts.keys()), rotation=90, fontsize=8)
        ax1.set_ylabel("Count")
        ax1.set_title("Amino-Acid Substitution Frequencies")

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
        print(f"  → Saved plot to {hist_png}")
        plt.show()
        plt.close()

        # threshold
        thr = int(input(f"Threshold for {base}: "))
        frequent = {aa for aa,c in counts.items() if c >= thr}
        print(f"  → {len(frequent)} substitutions ≥ {thr}")

        # save frequent counts
        freq_csv = os.path.join(out_dir, f"{sample}_frequent_aa_counts.csv")
        pd.DataFrame(
            [(aa, counts[aa]) for aa in sorted(frequent)],
            columns=["AA","Count"]
        ).to_csv(freq_csv, index=False)
        print(f"  → Saved frequent counts to {freq_csv}")

        # co-occurrence
        co_base, co_fish = find_cooccurring_aa(subs_aa, frequent, out_dir, sample)
        print(f"  → Co-occurrence files: {co_base}, {co_fish}")

        # write report
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
