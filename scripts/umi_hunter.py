#### Tried and tested script for identifying indels and substitutions from nanopore data using UMIs, taking all the .fastq.gz files in the wd. It also requires a reference sequence. ####


import gzip
import glob
import re
import csv
import os
import tempfile
import subprocess
from collections import Counter
from tqdm import tqdm
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline

# Load flanking data
flanks_df = pd.read_csv("data/umi_hunter.csv")

# UMI flanks and limits
umi_start = flanks_df.loc[0, 'umi_flanks']
umi_end = flanks_df.loc[1, 'umi_flanks']
umi_min = int(flanks_df.loc[0, 'umi_min_max'])
umi_max = int(flanks_df.loc[1, 'umi_min_max'])

# Gene flanks
gene_start = flanks_df.loc[0, 'gene_flanks']
gene_end = flanks_df.loc[1, 'gene_flanks']

# Compile regex patterns
pattern_umi = re.compile(rf'{umi_start}([ACGT]{{{umi_min},{umi_max}}}){umi_end}', re.IGNORECASE)
pattern_gene = re.compile(rf'{gene_start}(.*?){gene_end}', re.IGNORECASE)

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]

def extract_read_info(seq):
    """
    Attempt to extract both UMI and gene region from a sequence.
    First, try the forward orientation; if either is not found,
    try the reverse complement.
    Returns a tuple (umi, gene) if both are found, otherwise (None, None).
    """
    umi_match = pattern_umi.search(seq)
    gene_match = pattern_gene.search(seq)
    if umi_match and gene_match:
        return umi_match.group(1), gene_match.group(1)
    # Try reverse complement if forward fails
    rev_seq = reverse_complement(seq)
    umi_match = pattern_umi.search(rev_seq)
    gene_match = pattern_gene.search(rev_seq)
    if umi_match and gene_match:
        return umi_match.group(1), gene_match.group(1)
    return None, None

def process_fastq(file_path):
    """
    Process a FASTQ file and extract UMI and gene sequences from each read.
    Returns the total read count and a dictionary mapping each extracted UMI
    to a list of gene sequences observed in reads with that UMI.
    """
    read_count = 0
    umi_info = {}  # key: UMI string, value: list of gene sequences
    with gzip.open(file_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                break  # end of file
            seq = f.readline().strip()
            f.readline()  # skip '+' line
            f.readline()  # skip quality scores line
            read_count += 1

            umi, gene = extract_read_info(seq)
            if umi and gene:
                umi_info.setdefault(umi, []).append(gene)
    return read_count, umi_info

def levenshtein(s1, s2):
    """Compute the Levenshtein edit distance between two strings."""
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def percent_identity(seq1, seq2):
    """Calculate percent identity based on Levenshtein distance."""
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return 1.0
    dist = levenshtein(seq1, seq2)
    return (max_len - dist) / max_len

def cluster_umis(umi_info, threshold=0.9):
    """
    Cluster UMIs using a greedy algorithm.
    Each cluster is a dictionary with keys:
      'rep': representative UMI,
      'total_count': total number of reads (gene sequences) in the cluster,
      'members': dictionary of original UMIs and their counts,
      'gene_seqs': concatenated list of gene sequences from all members.
    """
    sorted_umis = sorted(umi_info.items(), key=lambda x: len(x[1]), reverse=True)
    clusters = []
    for umi, gene_list in sorted_umis:
        count = len(gene_list)
        added = False
        for cluster in clusters:
            rep = cluster['rep']
            if percent_identity(umi, rep) >= threshold:
                cluster['total_count'] += count
                cluster['members'][umi] = count
                cluster['gene_seqs'].extend(gene_list)
                added = True
                break
        if not added:
            clusters.append({
                'rep': umi,
                'total_count': count,
                'members': {umi: count},
                'gene_seqs': list(gene_list)
            })
    return clusters

def generate_custom_consensus_mafft(gene_seqs, verbose=True, mutation_threshold=0.5):
    """
    Generate a consensus sequence using MAFFT with a custom, gap-aware consensus caller.
    
    The function:
      1. Loads a reference from template.fasta.
      2. Writes the reference and the gene sequences to a temporary FASTA file.
      3. Runs MAFFT to align the sequences.
      4. For each column, compares the non-reference bases to the reference base.
         If a base (or gap) different from the reference is observed in at least
         'mutation_threshold' fraction of the non-reference sequences, that base is called;
         otherwise, the reference base is used.
    
    Returns a tuple (consensus_sequence, _).
    (The second value was previously the median gene length, but is no longer used.)
    """
    if verbose:
        print("  [MAFFT] Aligning {} gene sequences plus reference...".format(len(gene_seqs)))
        
    # Load the reference sequence from template.fasta.
    try:
        template_record = next(SeqIO.parse("data/template.fasta", "fasta"))
    except Exception as e:
        raise RuntimeError("Failed to load template.fasta: " + str(e))
        
    # Create a temporary FASTA file with the reference and gene sequences.
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_in:
        fasta_in = tmp_in.name
        template_record.id = "REF_template"
        template_record.description = ""
        SeqIO.write(template_record, tmp_in, "fasta")
        for i, seq in enumerate(gene_seqs):
            record = SeqRecord(Seq(seq), id="seq{}".format(i), description="")
            SeqIO.write(record, tmp_in, "fasta")
    
    # Run MAFFT on the temporary FASTA file.
    mafft_cline = MafftCommandline(input=fasta_in)
    try:
        stdout, stderr = mafft_cline()
    except Exception as e:
        os.remove(fasta_in)
        raise RuntimeError("MAFFT failed: {}".format(e))
    
    # Write the alignment output to a temporary file.
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as align_tmp:
        align_filename = align_tmp.name
        align_tmp.write(stdout)
    
    # Parse the alignment.
    alignment = AlignIO.read(align_filename, "fasta")
    # Find the reference record.
    ref_record = None
    other_records = []
    for record in alignment:
        if record.id == "REF_template":
            ref_record = record
        else:
            other_records.append(record)
    if ref_record is None or not other_records:
        consensus = ""
    else:
        consensus_chars = []
        num_reads = len(other_records)
        alignment_length = alignment.get_alignment_length()
        for i in range(alignment_length):
            ref_base = ref_record.seq[i]
            col_bases = [record.seq[i] for record in other_records]
            counts = Counter(col_bases)
            most_common, count = counts.most_common(1)[0]
            freq = count / num_reads
            if most_common != ref_base and freq >= mutation_threshold:
                consensus_chars.append(most_common)
            else:
                consensus_chars.append(ref_base)
        consensus = "".join(consensus_chars)
    
    os.remove(fasta_in)
    os.remove(align_filename)
    if verbose:
        print("  [MAFFT] Alignment and custom consensus complete.")
    return consensus, None

def write_umi_csv(output_file, clusters):
    """
    Write UMI cluster information to CSV.
    Columns:
      - Cluster Representative (UMI)
      - Total Count
      - Members (UMI:count pairs)
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Cluster Representative", "Total Count", "Members"])
        for cluster in clusters:
            members_str = "; ".join(["{}:{}".format(umi, cnt) for umi, cnt in cluster['members'].items()])
            writer.writerow([cluster['rep'], cluster['total_count'], members_str])

def write_gene_csv(output_file, clusters, verbose=True):
    """
    Write gene consensus information to CSV for clusters with >10 reads.
    Columns:
      - Cluster Representative (UMI)
      - Total Count
      - Consensus Gene (gaps removed)
      - Length Difference (consensus length - reference length)
      - Members (UMI:count pairs)
    Also returns a list of consensus FASTA records.
    """
    # Load the reference to get its length (remove any gaps)
    try:
        ref_record = next(SeqIO.parse("data/template.fasta", "fasta"))
        ref_seq = str(ref_record.seq).replace("-", "")
        ref_length = len(ref_seq)
    except Exception as e:
        raise RuntimeError("Failed to load template.fasta: " + str(e))
    
    consensus_records = []
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Cluster Representative", "Total Count", "Consensus Gene", "Length Difference", "Members"])
        clusters_to_align = [c for c in clusters if c['total_count'] > 10]
        for i, cluster in enumerate(tqdm(clusters_to_align, desc="Processing UMI clusters", unit="cluster")):
            if verbose:
                print("Processing cluster {}/{} (Rep: {}, Count: {})".format(
                    i+1, len(clusters_to_align), cluster['rep'], cluster['total_count']))
            consensus_gene, _ = generate_custom_consensus_mafft(cluster['gene_seqs'], verbose=verbose, mutation_threshold=0.7)
            # Remove gap characters from the consensus
            consensus_clean = consensus_gene.replace("-", "")
            length_diff = len(consensus_clean) - ref_length
            members_str = "; ".join(["{}:{}".format(umi, cnt) for umi, cnt in cluster['members'].items()])
            writer.writerow([cluster['rep'], cluster['total_count'], consensus_clean, length_diff, members_str])
            record_id = "{}_cluster{}".format(cluster['rep'], i+1)
            consensus_records.append(SeqRecord(Seq(consensus_clean), id=record_id,
                                               description="Length diff: {}".format(length_diff)))
    return consensus_records

def main():
    fastq_files = glob.glob("data/*.fastq.gz")
    if not fastq_files:
        print("No .fastq.gz files found in data/.")
        return

    # Process each FASTQ file
    for fastq_file in fastq_files:
        print("\nProcessing file: {}".format(fastq_file))
        reads, umi_info = process_fastq(fastq_file)
        print("Total reads processed: {}".format(reads))
        clusters = cluster_umis(umi_info, threshold=0.9)
        print("Total UMI clusters found: {}".format(len(clusters)))

        base = os.path.basename(fastq_file)
        if base.endswith(".gz"):
            base = os.path.splitext(base)[0]
        if base.endswith(".fastq"):
            base = os.path.splitext(base)[0]
        
        umi_csv = "results/{}_UMI_clusters.csv".format(base)
        write_umi_csv(umi_csv, clusters)
        print("UMI clusters written to {}".format(umi_csv))

        gene_csv = "results/{}_gene_consensus.csv".format(base)
        consensus_records = write_gene_csv(gene_csv, clusters, verbose=True)
        print("Gene consensus data (for clusters >10 reads) written to {}".format(gene_csv))

        # Write a FASTA file containing all consensus sequences for this FASTQ file.
        fasta_out = "results/{}_consensuses.fasta".format(base)
        SeqIO.write(consensus_records, fasta_out, "fasta")
        print("Consensus sequences saved to {}\n".format(fasta_out))

if __name__ == "__main__":
    main()
