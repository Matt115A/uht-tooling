import gzip
import glob
import re
import os
import tempfile
import subprocess
from collections import defaultdict, Counter
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from scipy.stats import fisher_exact
from tqdm import tqdm

# Load flanking data
flanks_df = pd.read_csv("data/mutation_caller/mutation_caller.csv")
gene_start = flanks_df.loc[0, 'gene_flanks']
gene_end = flanks_df.loc[1, 'gene_flanks']
gene_min = int(flanks_df.loc[0, 'gene_min_max'])
gene_max = int(flanks_df.loc[1, 'gene_min_max'])

pattern_gene = re.compile(rf'{gene_start}([ACGTNacgtn]{{{gene_min},{gene_max}}}){gene_end}', re.IGNORECASE)

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def extract_gene(seq):
    match = pattern_gene.search(seq)
    if match:
        return match.group(1)
    match = pattern_gene.search(reverse_complement(seq))
    if match:
        return match.group(1)
    return None

def process_fastq(file_path):
    gene_reads = {}
    with gzip.open(file_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            gene = extract_gene(seq)
            if gene and gene_min <= len(gene) <= gene_max:
                gene_reads[header.strip()[1:]] = gene
    return gene_reads

def align_to_reference(gene_seqs, reference):
    aligned = {}
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        ref_record = SeqRecord(Seq(reference), id="REF", description="")
        SeqIO.write(ref_record, tmp, "fasta")
        for rid, seq in gene_seqs.items():
            SeqIO.write(SeqRecord(Seq(seq), id=rid, description=""), tmp, "fasta")
        input_path = tmp.name

    mafft_cline = MafftCommandline(input=input_path)
    stdout, _ = mafft_cline()
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_out:
        tmp_out.write(stdout)
        output_path = tmp_out.name

    alignments = AlignIO.read(output_path, "fasta")
    ref_seq = None
    for record in alignments:
        if record.id == "REF":
            ref_seq = str(record.seq)
            break

    for record in alignments:
        if record.id == "REF":
            continue
        aligned[record.id] = str(record.seq)

    os.remove(input_path)
    os.remove(output_path)
    return ref_seq, aligned

def identify_substitutions(ref, aligned_reads):
    subs_by_read = defaultdict(list)

    # Step 1: Build map from alignment index to ungapped reference index
    alignment_to_ref_index = {}
    ref_index_to_alignment = {}
    ref_seq_clean = []
    ref_nt_index = 0
    for i, base in enumerate(ref):
        if base != "-":
            alignment_to_ref_index[i] = ref_nt_index
            ref_index_to_alignment[ref_nt_index] = i
            ref_seq_clean.append(base)
            ref_nt_index += 1
        else:
            alignment_to_ref_index[i] = None

    ref_seq_clean = "".join(ref_seq_clean)

    for rid, aligned_seq in aligned_reads.items():
        for align_pos, query_base in enumerate(aligned_seq):
            ref_nt_pos = alignment_to_ref_index.get(align_pos)
            if ref_nt_pos is None:
                continue  # this position is a gap in the reference

            ref_base = ref_seq_clean[ref_nt_pos]
            if query_base == "-" or query_base == ref_base:
                continue  # skip gaps and matches

            # Nucleotide mutation
            nt_mut = f"{ref_base}{ref_nt_pos + 1}{query_base}"

            # Amino acid mutation
            codon_index = ref_nt_pos // 3
            codon_start_nt = codon_index * 3
            codon_end_nt = codon_start_nt + 3
            if codon_end_nt > len(ref_seq_clean):
                continue

            ref_codon = ref_seq_clean[codon_start_nt:codon_end_nt]
            if len(ref_codon) != 3:
                continue

            read_codon = list(ref_codon)
            pos_in_codon = ref_nt_pos % 3
            read_codon[pos_in_codon] = query_base

            try:
                aa_from = str(Seq(ref_codon).translate())
                aa_to = str(Seq("".join(read_codon)).translate())
            except Exception:
                aa_from, aa_to = "?", "?"

            aa_mut = f"{aa_from}{codon_index + 1}{aa_to}"
            subs_by_read[rid].append(f"{nt_mut} ({aa_mut})")

    return subs_by_read


def find_cooccurring(subs_by_read, frequent_subs):
    sub_list = sorted(frequent_subs)
    sub_idx = {s: i for i, s in enumerate(sub_list)}

    matrix = []
    for subs in subs_by_read.values():
        row = [0] * len(sub_list)
        for s in subs:
            if s in sub_idx:
                row[sub_idx[s]] = 1
        matrix.append(row)

    df = pd.DataFrame(matrix, columns=sub_list)
    cooccur = []

    for i in range(len(sub_list)):
        for j in range(i + 1, len(sub_list)):
            A = df.iloc[:, i]
            B = df.iloc[:, j]
            both = ((A == 1) & (B == 1)).sum()
            A_total = (A == 1).sum()
            B_total = (B == 1).sum()

            # Report only if they always co-occur or if they co-occur >2 times
            if both >= 2 or (both > 0 and both == A_total == B_total):
                cooccur.append((sub_list[i], sub_list[j], both, A_total, B_total))

    pd.DataFrame(cooccur, columns=["Sub1", "Sub2", "Both_Count", "Sub1_Count", "Sub2_Count"])\
        .to_csv("results/mutation_caller/cooccuring_mutations.csv", index=False)


    df = pd.DataFrame(matrix, columns=sub_list)
    cooccur = []
    for i in range(len(sub_list)):
        for j in range(i + 1, len(sub_list)):
            A = df.iloc[:, i]
            B = df.iloc[:, j]
            both = ((A == 1) & (B == 1)).sum()
            only_A = ((A == 1) & (B == 0)).sum()
            only_B = ((A == 0) & (B == 1)).sum()
            neither = ((A == 0) & (B == 0)).sum()

            table = [[both, only_A], [only_B, neither]]
            if all(v >= 0 for row in table for v in row):  # ensure all nonnegative
                try:
                    _, p = fisher_exact(table)
                    if p < 0.05:
                        cooccur.append((sub_list[i], sub_list[j], p))
                except Exception as e:
                    print(f"Skipping {sub_list[i]} vs {sub_list[j]}: {e}")

    pd.DataFrame(cooccur, columns=["Sub1", "Sub2", "p-value"]).to_csv("results/mutation_caller/cooccuring_mutations.csv", index=False)


def main():
    all_reads = {}
    for f in glob.glob("data/mutation_caller/*.fastq.gz"):
        print(f"Reading {f}...")
        reads = process_fastq(f)
        all_reads.update(reads)
    print(f"Total valid reads: {len(all_reads)}")

    # Load template
    template = str(next(SeqIO.parse("data/mutation_caller/mutation_caller_template.fasta", "fasta")).seq)
    ref_aligned, aligned = align_to_reference(all_reads, template)
    subs_by_read = identify_substitutions(ref_aligned, aligned)

    # Plot substitutions and show the histogram before asking for a threshold
    substitution_counts = Counter([s for subs in subs_by_read.values() for s in subs])
    os.makedirs("results", exist_ok=True)
    plt.figure(figsize=(12, 6))
    plt.bar(substitution_counts.keys(), substitution_counts.values())
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig("results/mutation_caller/substitution_rate.png")
    plt.show()

    # Now ask for threshold
    threshold = int(input("Enter substitution occurrence threshold: "))

    # Save frequent substitutions
    with open("results/mutation_caller/singles.csv", "w") as f:
        f.write("Substitution,Count\n")
        for sub, count in substitution_counts.items():
            if count >= threshold:
                f.write(f"{sub},{count}\n")
    frequent_subs = {s for s, c in substitution_counts.items() if c >= threshold}

    find_cooccurring(subs_by_read, frequent_subs)


if __name__ == "__main__":
    main()
