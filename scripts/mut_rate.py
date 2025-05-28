#!/usr/bin/env python3

import subprocess
import sys
import glob
import os
import shutil
import logging
import pysam
from Bio import SeqIO
import matplotlib.pyplot as plt
import math
from math import sqrt

# Try to import scipy for Z-test p-value calculation
try:
    from scipy.stats import norm  # for Z-test p-value
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
    print("Warning: scipy is not installed. Z-test p-values cannot be computed.")

# Define directories
DATA_DIR = os.path.join("data", "mutation_rate_calculator")
WORK_DIR = os.path.join(DATA_DIR, "tmp")
RESULTS_DIR = os.path.join("results", "mutation_rate_calculator")

# Ensure directories exist
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Configure logging
def setup_logging():
    log_file = os.path.join(RESULTS_DIR, "run.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
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
    logging.info(f"Running minimap2: {' '.join(cmd)} -> {sam_output}")
    with open(sam_output, "w") as out_sam:
        subprocess.run(cmd, stdout=out_sam, check=True)
    return sam_output


def load_single_sequence(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) != 1:
        logging.error(f"Expected exactly 1 record in {fasta_path}, found {len(records)}.")
        raise ValueError(f"Expected exactly 1 record in {fasta_path}, found {len(records)}.")
    logging.info(f"Loaded sequence {records[0].id} from {fasta_path}")
    return str(records[0].seq), records[0].id


def create_multi_fasta(chunks, out_fasta_prefix="plasmid_chunks"):
    out_fasta = os.path.join(WORK_DIR, f"{out_fasta_prefix}.fasta")
    logging.info(f"Writing {len(chunks)} chunks to {out_fasta}")
    with open(out_fasta, "w") as f:
        for i, seq in enumerate(chunks):
            f.write(f">chunk_{i+1}\n")
            f.write(str(seq) + "\n")
    return out_fasta


def compute_mismatch_stats_sam(sam_file, refs_dict):
    mismatch_data = { name: {"pos_rates":[0.0]*len(seq),"cov":[0]*len(seq),"mismatch":[0]*len(seq),"avg_mismatch_rate":0.0,"total_mismatches":0,"total_covered_bases":0} for name,seq in refs_dict.items() }
    logging.info(f"Computing mismatch stats for {sam_file}")
    samfile = pysam.AlignmentFile(sam_file, "r")
    for read in samfile.fetch():
        if read.is_unmapped or read.query_sequence is None:
            continue
        ref_name = samfile.get_reference_name(read.reference_id)
        if ref_name not in mismatch_data:
            continue
        data = mismatch_data[ref_name]
        ref_seq = refs_dict[ref_name]
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if read_pos is not None and ref_pos is not None and 0 <= ref_pos < len(ref_seq):
                data["cov"][ref_pos] += 1
                if read.query_sequence[read_pos].upper() != ref_seq[ref_pos].upper():
                    data["mismatch"][ref_pos] += 1
    samfile.close()
    for name, info in mismatch_data.items():
        total_mis = total_cov = 0
        pos_rates=[]
        for cov, mis in zip(info["cov"], info["mismatch"]):
            if cov>0:
                rate = mis/cov
                pos_rates.append(rate)
                total_mis += mis
                total_cov += cov
            else:
                pos_rates.append(0.0)
        info["pos_rates"] = pos_rates
        info["total_mismatches"] = total_mis
        info["total_covered_bases"] = total_cov
        info["avg_mismatch_rate"] = (total_mis/total_cov) if total_cov>0 else 0.0
        logging.info(f"{name}: mismatches={total_mis}, covered={total_cov}, avg_rate={info['avg_mismatch_rate']:.6f}")
    return mismatch_data


def z_test_two_proportions(mis1, cov1, mis2, cov2):
    if cov1==0 or cov2==0:
        return 0.0,1.0
    p1,p2=mis1/cov1,mis2/cov2
    p=(mis1+mis2)/(cov1+cov2)
    denom=math.sqrt(p*(1-p)*(1/cov1+1/cov2))
    if denom==0:
        return 0.0,1.0
    z=(p1-p2)/denom
    p_val=2*(1-norm.cdf(abs(z))) if HAVE_SCIPY else None
    logging.info(f"Z-test: z={z:.3f}, p-value={p_val}")
    return z,p_val


def main():
    setup_logging()
    # Input paths
    reads_pattern = os.path.join(DATA_DIR, "*.fastq.gz")
    ref_hit_fasta = os.path.join(DATA_DIR, "region_of_interest.fasta")
    plasmid_fasta = os.path.join(DATA_DIR, "plasmid.fasta")

    logging.info("Starting mutation rate analysis run")
    logging.info(f"Configuration: reads_pattern={reads_pattern}, ref_hit={ref_hit_fasta}, plasmid={plasmid_fasta}")

    # Load sequences
    hit_seq, hit_id = load_single_sequence(ref_hit_fasta)
    plasmid_seq, plasmid_id = load_single_sequence(plasmid_fasta)

    # Remove hit region
    idx = plasmid_seq.find(hit_seq)
    if idx == -1:
        logging.error("Gene region not found in plasmid")
        sys.exit(1)
    plasmid_no_gene = plasmid_seq[:idx] + plasmid_seq[idx+len(hit_seq):]

    # Chunk plasmid
    n_chunks=10
    length=len(plasmid_no_gene)
    size=length//n_chunks
    chunks=[plasmid_no_gene[i*size: (length if i==n_chunks-1 else (i+1)*size)] for i in range(n_chunks)]

    # Write and align
    chunks_fasta=create_multi_fasta(chunks)
    sam_chunks=run_minimap2(reads_pattern, chunks_fasta, "plasmid_chunks_alignment")
    sam_hit=run_minimap2(reads_pattern, ref_hit_fasta, "hit_alignment")

    # Stats
    chunk_refs={f"chunk_{i+1}":seq for i,seq in enumerate(chunks)}
    mismatch_chunks=compute_mismatch_stats_sam(sam_chunks, chunk_refs)
    bg_mis=sum(info['total_mismatches'] for info in mismatch_chunks.values())
    bg_cov=sum(info['total_covered_bases'] for info in mismatch_chunks.values())
    bg_rate=(bg_mis/bg_cov) if bg_cov else 0.0
    logging.info(f"Background total mismatches={bg_mis}, covered={bg_cov}, rate={bg_rate:.6f}")

    mismatch_hit=compute_mismatch_stats_sam(sam_hit, {hit_id:hit_seq})
    hit_info=mismatch_hit[hit_id]
    hit_mis=hit_info['total_mismatches']
    hit_cov=hit_info['total_covered_bases']
    hit_rate=hit_info['avg_mismatch_rate']
    logging.info(f"Hit total mismatches={hit_mis}, covered={hit_cov}, rate={hit_rate:.6f}")

    # Z-test & significance
    z,p=z_test_two_proportions(hit_mis, hit_cov, bg_mis, bg_cov)
    significant = (p is not None and p < 0.05 and hit_rate > bg_rate)
    if significant:
        diff=hit_rate-bg_rate
        excess_kb=diff*1e3
        avg_mut=diff*len(hit_seq)
        logging.info(f"Hit error rate significantly higher (p<{0.05}), excess_kb={excess_kb:.4f}, avg_mutations={avg_mut:.4f}")
    else:
        logging.info("No significant increase in error rate or insufficient coverage")

    # Print summary and plot
    print("Results logged to ", os.path.join(RESULTS_DIR, "run.log"))
    plt.figure()
    plt.plot(hit_info["pos_rates"], marker='.', linestyle='-')
    plt.title("Mismatch Rate per Position: Hit Gene")
    plt.xlabel("Position")
    plt.ylabel("Mismatch Rate")
    plt.show()

    labels, rates = zip(*[(name, info['avg_mismatch_rate']) for name, info in mismatch_chunks.items()])
    plt.figure()
    plt.bar(labels, rates)
    plt.title("Average Mismatch Rate per Plasmid Chunk")
    plt.xlabel("Chunk")
    plt.ylabel("Mismatch Rate")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    # Clean up
    if os.path.isdir(WORK_DIR):
        shutil.rmtree(WORK_DIR)

if __name__ == "__main__":
    main()
