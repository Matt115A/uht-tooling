import re
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
    """Translate a 3-nt codon to its amino acid."""
    return codon_table().get(cd.upper(), "?")

def pick_mutant_codon(wt_codon, target_aa):
    """
    Choose the codon encoding 'target_aa' that has minimal nucleotide differences
    from the wild-type codon 'wt_codon'.
    """
    best_list = []
    for codon, aa in codon_table().items():
        if aa == target_aa:
            diff = sum(a != b for a, b in zip(codon.upper(), wt_codon.upper()))
            best_list.append((codon.upper(), diff))
    if not best_list:
        return None
    best_list.sort(key=lambda x: x[1])
    return best_list[0][0]

def calc_tm(seq):
    return mt.Tm_NN(seq)

# -----------------------------------------
# Primer design core functions
# -----------------------------------------
def pick_downstream_matched(dna_segment):
    best_sub, best_diff = "", float("inf")
    for length in range(1, len(dna_segment)+1):
        sub = dna_segment[:length]
        if len(sub) > MAX_LEN:
            break
        Tm = calc_tm(sub)
        if Tm <= MAX_TM:
            diff = abs(Tm - TARGET_TM)
            if diff < best_diff:
                best_diff, best_sub = diff, sub
    return best_sub

def pick_upstream_matched(dna_rev):
    best_sub, best_diff = "", float("inf")
    for length in range(1, len(dna_rev)+1):
        sub = dna_rev[:length]
        if len(sub) > MAX_LEN:
            break
        Tm = calc_tm(sub)
        if Tm <= MAX_TM:
            diff = abs(Tm - TARGET_TM)
            if diff < best_diff:
                best_diff, best_sub = diff, sub
    return best_sub

def design_long_forward(full_seq, region_start, old_len, new_codon):
    """
    Design the Long Forward primer (Lf).
    - region_start: 0-based index where original region begins
    - old_len: number of nucleotides deleted in original
    - new_codon: nucleotide string to insert
    """
    # ensure sufficient upstream
    if region_start < UPSTREAM_15:
        raise ValueError("Not enough upstream for Lf (15 bp).")
    # take 15 nt upstream
    lf_up = full_seq[region_start - UPSTREAM_15 : region_start]
    # downstream region starts after the deleted region
    down_seg = full_seq[region_start + old_len :]
    # pick segment giving Tm closest to TARGET_TM
    matched_down = pick_downstream_matched(down_seg)
    # concatenate upstream + mutation + downstream
    lf_seq = lf_up + new_codon + matched_down
    L = len(lf_seq)
    if L < MIN_LEN or L > MAX_LEN:
        raise ValueError(f"Lf length {L} not in [{MIN_LEN}..{MAX_LEN}].")
    lf_start = region_start - UPSTREAM_15
    lf_end = lf_start + L - 1
    return lf_seq.upper(), lf_start, lf_end

def design_long_reverse(full_seq, region_start, old_len, new_codon):
    """
    Design the Long Reverse primer (Lr).
    - region_start: 0-based index where original region begins
    - old_len: number of nucleotides deleted in original
    - new_codon: nucleotide string to insert
    """
    # reverse of sequence upstream of region_start
    up_rev = full_seq[:region_start][::-1]
    matched_up = pick_upstream_matched(up_rev)
    up_fwd = matched_up[::-1]
    # primer covers inserted mutation at 3' end
    lr_fwd = up_fwd + new_codon
    L = len(lr_fwd)
    if L < MIN_LEN or L > MAX_LEN:
        raise ValueError(f"Lr length {L} not in [{MIN_LEN}..{MAX_LEN}].")
    # start index of primer
    lr_start = region_start - len(up_fwd)
    # primer ends at original region end (inclusive)
    lr_end = region_start + old_len - 1
    if lr_start < 0 or lr_end >= len(full_seq):
        raise ValueError("Lr coverage out of range.")
    # reverse complement to get actual primer sequence
    lr_final = str(Seq(lr_fwd).reverse_complement()).upper()
    return lr_final, lr_start, lr_end

def design_short_reverse(full_seq, lf_start):
    """Design the Short Reverse primer (Sr) ending just before Lf."""
    sr_end = lf_start - 1
    if sr_end < 0:
        raise ValueError("No space for Sr (lf_start=0).")
    best_seg, best_diff = "", float("inf")
    for length in range(MIN_LEN, MAX_LEN+1):
        start = sr_end - length + 1
        if start < 0:
            break
        region = full_seq[start : sr_end+1]
        Tm = calc_tm(region)
        if Tm <= MAX_TM:
            diff = abs(Tm - TARGET_TM)
            if diff < best_diff:
                best_diff, best_seg = diff, region
    if not best_seg:
        raise ValueError("Cannot find Sr [20..50] with Tm <=70 upstream.")
    return str(Seq(best_seg).reverse_complement()).upper()

def design_short_forward(full_seq, lr_end):
    """Design the Short Forward primer (Sf) starting just after Lr."""
    sf_start = lr_end + 1
    if sf_start >= len(full_seq):
        raise ValueError("No space for Sf (lr_end near end).")
    best_seg, best_diff = "", float("inf")
    max_possible = len(full_seq) - sf_start
    for length in range(MIN_LEN, min(MAX_LEN, max_possible)+1):
        sub = full_seq[sf_start : sf_start + length]
        Tm = calc_tm(sub)
        if Tm <= MAX_TM:
            diff = abs(Tm - TARGET_TM)
            if diff < best_diff:
                best_diff, best_seg = diff, sub
    if not best_seg:
        raise ValueError("Cannot find Sf [20..50] with Tm <=70.")
    return best_seg.upper()

# -----------------------------------------
# Main workflow
# -----------------------------------------
def main():
    # Load gene and context sequences from FASTA
    gene = str(next(SeqIO.parse("data/design_slim/slim_template_gene.fasta", "fasta")).seq).upper()
    context = str(next(SeqIO.parse("data/design_slim/slim_context.fasta", "fasta")).seq).upper()

    # Load mutations from CSV
    df = pd.read_csv("data/design_slim/slim_target_mutations.csv")
    mutations = df["mutations"].dropna().tolist()

    # Align gene within context
    try:
        gene_offset = context.index(gene)
    except ValueError:
        raise ValueError("Could not align gene within context. No perfect substring match found.")
    full_seq = context
    gene_len = len(gene)

    os.makedirs("results", exist_ok=True)
    with open("results/design_slim/SLIM_primers.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Primer Name", "Sequence"])

        for mutation in mutations:
            # Parse mutation according to nomenclature:
            # Substitution:        A123G
            # Deletion:            T241Del
            # InDel (inter-codon): T241InDelA242S
            # Insertion after codon: T241TS (insert Ser after Thr241)
            # Codon replacement insertion: L46GP (replace Leu46 with Gly-Pro)
            m = mutation
            m_del   = re.match(r'^([A-Z])(\d+)Del$', m)
            m_indel = re.match(r'^([A-Z])(\d+)InDel([A-Z])(\d+)([A-Z]+)$', m)
            m_sub   = re.match(r'^([A-Z])(\d+)([A-Z])$', m)
            m_ins   = re.match(r'^([A-Z])(\d+)([A-Z]{2,})$', m)

            if m_del:
                wt_aa, pos1 = m_del.group(1), int(m_del.group(2))
                region_start = gene_offset + (pos1 - 1) * 3
                old_len = 3
                new_seq = ''
            elif m_indel:
                wt1, pos1, wt2, pos2, ins_aa = m_indel.groups()
                pos1, pos2 = int(pos1), int(pos2)
                region_start = gene_offset + (pos1 - 1) * 3
                old_len = (pos2 - pos1 + 1) * 3
                wt_codon = full_seq[region_start:region_start+3]
                new_seq = ''.join(pick_mutant_codon(wt_codon, aa) or '' for aa in ins_aa)
            elif m_ins:
                wt_aa, pos1_s, ins_str = m_ins.groups()
                pos1 = int(pos1_s)
                codon_start_old = gene_offset + (pos1 - 1) * 3
                # determine if this is insertion after codon or replacement
                if ins_str[0] == wt_aa:
                    # insertion after codon: first letter = original aa
                    inserted_aas = ins_str[1:]
                    region_start = codon_start_old + 3
                    old_len = 0
                else:
                    # codon replacement insertion: replace original codon
                    inserted_aas = ins_str
                    region_start = codon_start_old
                    old_len = 3
                wt_codon = full_seq[codon_start_old:codon_start_old+3]
                new_seq = ''.join(
                    pick_mutant_codon(wt_codon, aa) or '' for aa in inserted_aas
                )
            elif m_sub:
                wt_aa, pos1, mut_aa = m_sub.group(1), int(m_sub.group(2)), m_sub.group(3)
                region_start = gene_offset + (pos1 - 1) * 3
                old_len = 3
                wt_codon = full_seq[region_start:region_start+3]
                if translate_codon(wt_codon) != wt_aa:
                    raise ValueError(f"For {mutation}: expected {wt_aa}, found {translate_codon(wt_codon)} at {wt_codon}")
                new_seq = pick_mutant_codon(wt_codon, mut_aa)
                if not new_seq:
                    raise ValueError(f"No minimal-change codon for {wt_aa}->{mut_aa}")
            else:
                raise ValueError(f"Unknown mutation format: {mutation}")

            # Design primers
            lf_seq, lf_start, lf_end = design_long_forward(full_seq, region_start, old_len, new_seq)
            lr_seq, lr_start, lr_end = design_long_reverse(full_seq, region_start, old_len, new_seq)
            sr_seq = design_short_reverse(full_seq, lf_start)
            sf_seq = design_short_forward(full_seq, lr_end)

            # Output
            writer.writerow([mutation + '_Lf', lf_seq])
            writer.writerow([mutation + '_Sr', sr_seq])
            writer.writerow([mutation + '_Lr', lr_seq])
            writer.writerow([mutation + '_Sf', sf_seq])

    print("results/design_slim/SLIM_primers.csv created successfully.")

if __name__ == "__main__":
    main()
