import re
import csv
import os
from math import floor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
from collections import defaultdict

# -----------------------------------------
# CONSTANTS for Gibson overlap–assembly primers
# -----------------------------------------
OVERHANG_LEN = 20   # total nt of homology (must contain the entire new_seq)
ANNEAL_LEN   = 20   # nt that actually anneal to template beyond that overlap

# -----------------------------------------
# CODON–translation utilities (unchanged)
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

def translate_codon(cd: str) -> str:
    """Translate a 3‐nt codon to its amino acid (or '?' if invalid)."""
    return codon_table().get(cd.upper(), "?")

def pick_mutant_codon(wt_codon: str, target_aa: str) -> str:
    """
    Among all codons encoding 'target_aa', pick the one with the fewest
    nucleotide differences from wt_codon. Returns a 3‐nt string.
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

# -----------------------------------------
# CIRCULAR‐sequence helper (wrap‐around)
# -----------------------------------------
def get_subseq_circ(seq: str, start: int, length: int) -> str:
    """
    Return `length` bases of `seq`, starting at index `start` (0-based),
    but wrap around if start < 0 or (start + length) > len(seq).  
    Raises ValueError if length >= len(seq).
    """
    N = len(seq)
    if length >= N:
        raise ValueError(f"Requested length {length} ≥ sequence length {N}")
    s_mod = start % N
    end = s_mod + length
    if end <= N:
        return seq[s_mod:end]
    else:
        return seq[s_mod:] + seq[: (end % N)]

# -----------------------------------------
# GIBSON–primer design for a single mutation (circular)
# -----------------------------------------
def design_gibson_primers(full_seq: str, region_start: int, old_len: int, new_seq: str):
    """
    Given a circular `full_seq` and a mutation defined by:
      • region_start (0‐based index of the first base to remove/replace),
      • old_len      (# of nucleotides to delete at region_start),
      • new_seq      (inserted/replacement bases; e.g. 3 nt for a codon, or '' for deletion),

    produce two Gibson primers (5′→3′) that introduce exactly that mutation.

    Each primer is:
      – 5′ overhang of length OVERHANG_LEN (containing new_seq),
      – 3′ annealing region of length ANNEAL_LEN immediately outside that overhang.

    Returns a tuple:
      (gibson_fwd_seq, gibson_rev_seq, fwd_start, rev_start)

    where fwd_start is the plus‐strand index of the first base that the forward primer’s
    3′‐annealing region binds, and rev_start is the plus‐strand index of the first base
    that the reverse primer’s 3′‐annealing region binds.
    """
    m = len(new_seq)
    if OVERHANG_LEN < m:
        raise ValueError(f"Mutation length {m} > OVERHANG_LEN ({OVERHANG_LEN})")

    N = len(full_seq)
    flank_left  = floor((OVERHANG_LEN - m)/2)
    flank_right = OVERHANG_LEN - m - flank_left

    if flank_left >= N or flank_right >= N or ANNEAL_LEN >= N:
        raise ValueError("Sequence too short for requested OVERHANG/ANNEAL lengths")

    # 1) Build the 5′ overhang (exactly OVERHANG_LEN):
    oh_left   = get_subseq_circ(full_seq, region_start - flank_left,  flank_left)
    oh_right  = get_subseq_circ(full_seq, region_start + old_len,    flank_right)
    overhang  = oh_left + new_seq + oh_right  # length = OVERHANG_LEN

    # 2) Compute 3′‐annealing start positions:
    fwd_start   = (region_start + old_len + flank_right) % N
    rev_start   = (region_start - flank_left - ANNEAL_LEN) % N

    fwd_anneal  = get_subseq_circ(full_seq, fwd_start, ANNEAL_LEN)
    rev_anneal  = get_subseq_circ(full_seq, rev_start, ANNEAL_LEN)

    # 3) Construct primers (5′→3′):
    gibson_fwd  = (overhang + fwd_anneal).upper()
    overhang_rc = str(Seq(overhang).reverse_complement()).upper()
    rev_anneal_rc = str(Seq(rev_anneal).reverse_complement()).upper()
    gibson_rev  = (overhang_rc + rev_anneal_rc).upper()

    return gibson_fwd, gibson_rev, fwd_start, rev_start

# -----------------------------------------
# MAIN workflow
# -----------------------------------------
def main():
    # 1) Load “template gene” and “context” (circular) from FASTA
    gene_fasta    = "data/design_gibson/gibson_template_gene.fasta"
    context_fasta = "data/design_gibson/gibson_context.fasta"

    gene_seq    = str(next(SeqIO.parse(gene_fasta,    "fasta")).seq).upper()
    context_seq = str(next(SeqIO.parse(context_fasta, "fasta")).seq).upper()

    # 2) Read mutations from CSV (column “mutations”)
    df = pd.read_csv("data/design_gibson/gibson_target_mutations.csv")
    entries = df["mutations"].dropna().tolist()
    # e.g. ["A123G", "T241Del", "A123G+T241Del", ...]

    # 3) Locate gene_seq in context_seq circularly by searching doubled string
    double_seq = context_seq + context_seq
    idx = double_seq.find(gene_seq)
    if idx == -1 or idx >= len(context_seq):
        raise ValueError("Could not align gene within circular context.")
    gene_offset = idx % len(context_seq)
    full_seq = context_seq  # now treat as circular

    # 4) Prepare output folder + filenames
    os.makedirs("results/design_gibson", exist_ok=True)
    primers_csv = "results/design_gibson/Gibson_primers.csv"
    plan_csv    = "results/design_gibson/Gibson_assembly_plan.csv"

    # For each group, accumulate a list of dicts with:
    #   sub, fwd_name, rev_name, fwd_pos, rev_pos, fwd_seq, rev_seq
    group_entries = defaultdict(list)

    # 5) Write Gibson_primers.csv (one primer per row)
    with open(primers_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Group", "Submutation", "Primer Name", "Sequence"])

        # 6) Loop through each CSV entry (may contain multiple sub‐mutations with “+”)
        for entry in entries:
            submuts    = entry.split('+')
            group_name = entry.replace('+', '_')  # e.g. "A123G_T241Del"

            for sub in submuts:
                # Parse the submutation string
                m_del   = re.match(r'^([A-Z])(\d+)Del$',               sub)
                m_indel = re.match(r'^([A-Z])(\d+)InDel([A-Z])(\d+)([A-Z]+)$', sub)
                m_sub   = re.match(r'^([A-Z])(\d+)([A-Z])$',             sub)
                m_ins   = re.match(r'^([A-Z])(\d+)([A-Z]{2,})$',         sub)

                if m_del:
                    wt_aa, pos1 = m_del.group(1), int(m_del.group(2))
                    region_start = gene_offset + (pos1 - 1)*3
                    old_len  = 3
                    new_seq  = ""
                elif m_indel:
                    wt1, pos1_s, wt2, pos2_s, ins_aa = m_indel.groups()
                    pos1, pos2 = int(pos1_s), int(pos2_s)
                    region_start = gene_offset + (pos1 - 1)*3
                    old_len  = (pos2 - pos1 + 1)*3
                    wt_codon = get_subseq_circ(full_seq, region_start, 3)
                    new_seq  = ""
                    for aa in ins_aa:
                        codon = pick_mutant_codon(wt_codon, aa)
                        if not codon:
                            raise ValueError(f"No codon found for {wt1}->{ins_aa}")
                        new_seq += codon
                elif m_ins:
                    wt_aa, pos1_s, ins_str = m_ins.groups()
                    pos1 = int(pos1_s)
                    codon_start_old = gene_offset + (pos1 - 1)*3
                    wt_codon = get_subseq_circ(full_seq, codon_start_old, 3)
                    if ins_str[0] == wt_aa:
                        # insertion _after_ that codon
                        inserted_aas = ins_str[1:]
                        region_start = codon_start_old + 3
                        old_len  = 0
                    else:
                        # codon‐replacement insertion
                        inserted_aas = ins_str
                        region_start = codon_start_old
                        old_len  = 3
                    new_seq = ""
                    for aa in inserted_aas:
                        codon = pick_mutant_codon(wt_codon, aa)
                        if not codon:
                            raise ValueError(f"No codon for insertion AA {aa}")
                        new_seq += codon
                elif m_sub:
                    wt_aa, pos1_s, mut_aa = m_sub.group(1), int(m_sub.group(2)), m_sub.group(3)
                    region_start = gene_offset + (pos1_s - 1)*3
                    old_len  = 3
                    wt_codon = get_subseq_circ(full_seq, region_start, 3)
                    if translate_codon(wt_codon) != wt_aa:
                        raise ValueError(f"For {sub}: expected {wt_aa}, found {translate_codon(wt_codon)}")
                    new_codon = pick_mutant_codon(wt_codon, mut_aa)
                    if not new_codon:
                        raise ValueError(f"No minimal‐change codon for {wt_aa}->{mut_aa}")
                    new_seq = new_codon
                else:
                    raise ValueError(f"Unknown mutation format: {sub}")

                # 7) Design Gibson primers (circularly) and get binding start positions:
                try:
                    gibson_fwd, gibson_rev, fwd_start, rev_start = design_gibson_primers(
                        full_seq, region_start, old_len, new_seq
                    )
                except ValueError as e:
                    raise ValueError(f"Error designing {sub} in group {entry}: {e}")

                # 8) Name them and write to Gibson_primers.csv
                primer_fwd_name = f"{group_name}__{sub}_Gibson_F"
                primer_rev_name = f"{group_name}__{sub}_Gibson_R"

                writer.writerow([group_name, sub, primer_fwd_name, gibson_fwd])
                writer.writerow([group_name, sub, primer_rev_name, gibson_rev])

                # 9) Record all relevant info under this group
                group_entries[group_name].append({
                    "sub":       sub,
                    "fwd_name":  primer_fwd_name,
                    "rev_name":  primer_rev_name,
                    "fwd_pos":   fwd_start % len(full_seq),
                    "rev_pos":   rev_start % len(full_seq),
                    "fwd_seq":   gibson_fwd,
                    "rev_seq":   gibson_rev
                })

    # 10) Write Gibson_assembly_plan.csv with one line per *PCR primer pair*
    #     Pairing rule: sort forwards by fwd_pos ascending, sort reverses by rev_pos ascending,
    #     then pair forward[i] with reverse[(i+1) % n].
    with open(plan_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Group",
            "Submutation",
            "PCR_Primer_Forward",
            "PCR_Primer_Reverse",
            "Tm (celsius)",              # lower Tm of the two primers in the pair
            "Amplicon Size (bp)"    # number of bases to be amplified
        ])

        for group_name, entries in group_entries.items():
            # Sort by forward‐binding position (smallest→largest)
            sorted_forwards = sorted(entries, key=lambda e: e["fwd_pos"])
            # Sort by reverse‐binding position (smallest→largest)
            sorted_reverses = sorted(entries, key=lambda e: e["rev_pos"])

            n = len(sorted_forwards)
            N = len(full_seq)
            for i in range(n):
                f_entry = sorted_forwards[i]
                r_entry = sorted_reverses[(i + 1) % n]

                fwd_name = f_entry["fwd_name"]
                rev_name = r_entry["rev_name"]
                submut   = f_entry["sub"]  # which submutation is associated with this forward

                # Compute Tm of each primer, then take the minimum:
                Tm_fwd = mt.Tm_NN(f_entry["fwd_seq"])
                Tm_rev = mt.Tm_NN(r_entry["rev_seq"])
                Tm_pair = min(Tm_fwd, Tm_rev)

                # Compute amplicon size: from fwd_start to rev_end (inclusive), circularly
                fwd_start = f_entry["fwd_pos"]
                rev_start = r_entry["rev_pos"]
                rev_end   = (rev_start + ANNEAL_LEN - 1) % N

                if rev_end >= fwd_start:
                    amp_size = rev_end - fwd_start + 1
                else:
                    amp_size = (N - fwd_start) + (rev_end + 1)

                writer.writerow([
                    group_name,
                    submut,
                    fwd_name,
                    rev_name,
                    f"{Tm_pair:.1f}",
                    amp_size
                ])

    print(f"Created Gibson primers CSV: {primers_csv}")
    print(f"Created Gibson assembly plan CSV: {plan_csv}")


if __name__ == "__main__":
    main()
