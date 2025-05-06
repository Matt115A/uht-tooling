import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
import os

# -----------------------------------------
# Basic constraints
# -----------------------------------------
MIN_LEN = 20
MAX_LEN = 55
TARGET_TM = 60.0
MAX_TM = 70.0
UPSTREAM_15 = 12

# -----------------------------------------
# Codon utilities
# -----------------------------------------
def codon_table():
    return {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
    }

def translate_codon(cd):
    return codon_table().get(cd.upper(), "?")

def pick_mutant_codon(wt_codon, target_aa):
    best_list = []
    for c, aa in codon_table().items():
        if aa == target_aa:
            diff = sum(a != b for a,b in zip(c.upper(), wt_codon.upper()))
            best_list.append((c.upper(), diff))
    if not best_list:
        return None
    best_list.sort(key=lambda x: x[1])
    return best_list[0][0]

def calc_tm(seq):
    return mt.Tm_NN(seq)

def pick_downstream_matched(dna_segment):
    best_sub, best_diff = "", float("inf")
    for length in range(1, len(dna_segment)+1):
        sub = dna_segment[:length]
        if len(sub) > MAX_LEN:
            break
        T = calc_tm(sub)
        if T <= MAX_TM:
            diff = abs(T - TARGET_TM)
            if diff < best_diff:
                best_diff, best_sub = diff, sub
    return best_sub

def pick_upstream_matched(dna_rev):
    best_sub, best_diff = "", float("inf")
    for length in range(1, len(dna_rev)+1):
        sub = dna_rev[:length]
        if len(sub) > MAX_LEN:
            break
        T = calc_tm(sub)
        if T <= MAX_TM:
            diff = abs(T - TARGET_TM)
            if diff < best_diff:
                best_diff, best_sub = diff, sub
    return best_sub

def design_long_forward(full_seq, codon_start, new_codon):
    if codon_start < UPSTREAM_15:
        raise ValueError("Not enough upstream for Lf (15 bp).")
    lf_up = full_seq[codon_start - UPSTREAM_15 : codon_start]
    down_seg = full_seq[codon_start+3 :]
    matched_down = pick_downstream_matched(down_seg)
    lf_seq = lf_up + new_codon + matched_down
    L = len(lf_seq)
    if L < MIN_LEN or L > MAX_LEN:
        raise ValueError(f"Lf length {L} not in [20..50].")
    lf_start = codon_start - UPSTREAM_15
    lf_end = lf_start + L - 1
    return lf_seq.upper(), lf_start, lf_end

def design_long_reverse(full_seq, codon_start, new_codon):
    up_rev = full_seq[:codon_start][::-1]
    matched_up = pick_upstream_matched(up_rev)
    up_fwd = matched_up[::-1]
    lr_fwd = up_fwd + new_codon
    L = len(lr_fwd)
    if L < MIN_LEN or L > MAX_LEN:
        raise ValueError(f"Lr length {L} not in [20..50].")
    lr_start = codon_start - len(up_fwd)
    lr_end = codon_start + 2
    if lr_start < 0 or lr_end >= len(full_seq):
        raise ValueError("Lr coverage out of range.")
    lr_final = str(Seq(lr_fwd).reverse_complement()).upper()
    return lr_final, lr_start, lr_end

def design_short_reverse(full_seq, lf_start):
    sr_end = lf_start - 1
    if sr_end < 0:
        raise ValueError("No space for Sr (lf_start=0).")
    best_seg, best_diff = "", float("inf")
    for length in range(MIN_LEN, MAX_LEN+1):
        start = sr_end - length + 1
        if start < 0:
            break
        region = full_seq[start : sr_end+1]
        T = calc_tm(region)
        if T <= MAX_TM:
            diff = abs(T - TARGET_TM)
            if diff < best_diff:
                best_diff, best_seg = diff, region
    if not best_seg:
        raise ValueError("Cannot find Sr [20..50] with Tm <=70 upstream.")
    return str(Seq(best_seg).reverse_complement()).upper()

def design_short_forward(full_seq, lr_end):
    sf_start = lr_end + 1
    if sf_start >= len(full_seq):
        raise ValueError("No space for Sf (lr_end near end).")
    best_seg, best_diff = "", float("inf")
    max_possible = len(full_seq) - sf_start
    for length in range(MIN_LEN, min(MAX_LEN, max_possible)+1):
        sub = full_seq[sf_start : sf_start + length]
        T = calc_tm(sub)
        if T <= MAX_TM:
            diff = abs(T - TARGET_TM)
            if diff < best_diff:
                best_diff, best_seg = diff, sub
    if not best_seg:
        raise ValueError("Cannot find Sf [20..50] with Tm <=70.")
    return best_seg.upper()

def main():
    # Load gene and context sequences from FASTA
    gene = str(next(SeqIO.parse("data/slim_template_gene.fasta", "fasta")).seq).upper()
    context = str(next(SeqIO.parse("data/slim_context.fasta", "fasta")).seq).upper()

    # Load mutations from CSV column
    df = pd.read_csv("data/slim_target_mutations.csv")
    mutations = df["mutations"].dropna().tolist()

    # Align gene in context
    try:
        gene_offset = context.index(gene)
    except ValueError:
        raise ValueError("Could not align gene within context. No perfect substring match found.")

    full_seq = context
    gene_len = len(gene)

    os.makedirs("results", exist_ok=True)
    with open("results/SLIM_primers.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Primer Name", "Sequence"])

        for mutation in mutations:
            wt_aa = mutation[0]
            pos = int(mutation[1:-1])
            mut_aa = mutation[-1]

            codon_start_in_gene = (pos - 1) * 3
            if codon_start_in_gene + 3 > gene_len:
                raise ValueError(f"Mutation {mutation} out of range for gene length {gene_len}.")

            codon_start = gene_offset + codon_start_in_gene
            if codon_start + 3 > len(full_seq):
                raise ValueError(f"Mutation {mutation} out of range in context.")

            wt_codon = full_seq[codon_start : codon_start + 3]
            found_aa = translate_codon(wt_codon)
            if found_aa != wt_aa:
                raise ValueError(f"For mutation {mutation}, expected {wt_aa}, found {found_aa} at codon {wt_codon}")

            new_codon = pick_mutant_codon(wt_codon, mut_aa)
            if not new_codon:
                raise ValueError(f"No minimal-change codon for {wt_aa}->{mut_aa}.")

            lf_seq, lf_start, lf_end = design_long_forward(full_seq, codon_start, new_codon)
            lr_seq, lr_start, lr_end = design_long_reverse(full_seq, codon_start, new_codon)
            sr_seq = design_short_reverse(full_seq, lf_start)
            sf_seq = design_short_forward(full_seq, lr_end)

            writer.writerow([mutation + "_Lf", lf_seq])
            writer.writerow([mutation + "_Sr", sr_seq])
            writer.writerow([mutation + "_Lr", lr_seq])
            writer.writerow([mutation + "_Sf", sf_seq])

    print("results/SLIM_primers.csv created successfully.")

if __name__ == "__main__":
    main()
