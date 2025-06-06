import csv
import yaml
import pandas as pd


##### INPUT: In /data/nextera_designer/, place a file called nextera_designer.csv, that has the template binding sequences in the first column, which should be named binding_region. These should be the PCR primers you intend to use. Then run this script ######

# === LOAD DATA ===
df = pd.read_csv("data/nextera_designer/nextera_designer.csv")

# === USER SETTINGS ===
TEMPLATE_BINDING_SEQUENCE_i7 = df['binding_region'].iloc[0]  # First entry
TEMPLATE_BINDING_SEQUENCE_i5 = df['binding_region'].iloc[1]  # Second entry

OUTPUT_FILE = "results/nextera_designer/nextera_xt_primers.csv"

# === Index definitions ===
i7_indexes = {
    "701": "TCGCCTTA",
    "702": "CTAGTACG",
    "703": "TTCTGCCT",
    "704": "GCTCAGGA",
    "705": "AGGAGTCC",
    "706": "CATGCCTA",
    "707": "GTAGAGAG",
    "708": "CCTCTCTG",
    "709": "AGCGTAGC",
    "710": "CAGCCTCG",
    "711": "TGCCTCTT",
    "712": "TCCTCTAC"
}

i5_indexes = {
    "501": "TAGATCGC",
    "502": "CTCTCTAT",
    "503": "TATCCTCT",
    "504": "AGAGTAGA",
    "505": "GTAAGGAG",
    "506": "ACTGCATA",
    "507": "AAGGAGTA",
    "508": "CTAAGCCT",
    "510": "CGTCTAAT",
    "511": "TCTCTCCG",
    "513": "TCGACTAG",
    "515": "TTCTAGCT"
}

# === Constant primer structures ===
i7_prefix = "CAAGCAGAAGACGGCATACGAGAT"
i7_suffix = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

i5_prefix = "AATGATACGGCGACCACCGAGATCTACAC"
i5_suffix = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"

# === Primer generation ===
def generate_primers():
    primers = []

    for idx, seq in i7_indexes.items():
        name = f"i7_{idx}"
        full_seq = f"{i7_prefix}{seq}{i7_suffix}{TEMPLATE_BINDING_SEQUENCE_i7}"
        primers.append((name, full_seq))

    for idx, seq in i5_indexes.items():
        name = f"i5_{idx}"
        full_seq = f"{i5_prefix}{seq}{i5_suffix}{TEMPLATE_BINDING_SEQUENCE_i5}"
        primers.append((name, full_seq))

    return primers

def main():
    with open("configs/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    primers = generate_primers()

    with open(OUTPUT_FILE, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["primer_name", "sequence"])  # header
        for name, seq in primers:
            writer.writerow([name, seq])

    print(f"Saved {len(primers)} primers to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
