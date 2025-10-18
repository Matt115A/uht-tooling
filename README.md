# uht-tooling: Automating the Molecular Biology Accessory to Ultra-High Throughput Screening

---

## Getting Started

You will need **conda** installed. Follow the instructions here:
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Clone the repo into a directory of your choice:

```bash
	git clone https://github.com/Matt115A/uht-tooling.git
```

Navigate to `/uht-tooling/` and run:

```bash
	chmod +x setup.sh
	./setup.sh
```

*Note: This setup.sh file is verified for Mac only. Results may vary on other OS.*

---

## General Structure

- Configs and data are found in `data/.../`, where `...` depends on the script you wish to run.
- Results are saved to `results/.../`.
- All major scripts now generate detailed log files in their respective `results/.../` directories. These logs capture all major steps, errors, and outputs, making it easier to debug and trace your runs.

### Available Commands

- `make nextera_primers`
- `make umi_hunter`
- `make design_slim`
- `make mutation_caller`
- `make ep-library-profile`
- `make design_gibson`
- `make profile_inserts`

---

## Automated Nextera XT Primer Design Software

**One-PCR-to-flowcell workflow**

To use this software:
1. Upload a `.csv` file in `data/nextera_designer/` called `nextera_designer.csv`. The first column should be called `binding_region`, and the first two rows should be the forward and reverse primers you want to use to generate your Illumina amplicon.
2. Use normal primer design rules, inputting the sequence in the 5'→3' order (same as if you were ordering primers).
3. The script is pre-loaded with twelve i5 and twelve i7 indices, allowing up to 144 unique amplicons.
4. After setting up your conda environment, run `make nextera_primers`. A new `.csv` will appear in `results/nextera_designer/`, ready for ordering.

**Pipeline for kit-free Illumina sample prep:**
- PCR using an i5/i7 primer pair, monitor a portion of the reaction using qPCR. Cap the cycle number to get only 10% of the final yield (minimises PCR amplification bias).
- Purify the DNA product, removing primers (important to avoid primer excess on the chip). Use SPRIselect beads from Beckman Coulter. Start with a 0.65:1 bead:DNA volume ratio.
- Verify primer removal using electrophoresis (e.g., 2100 BioAnalyser and a DNA detection chip).

---

## Profile Inserts: Analyzing Insert Sequences from FASTQ Data

Extract and analyze insert sequences from FASTQ data using upstream and downstream probe sequences:

- Place your `.fastq.gz` file and a `.csv` file with `upstream` and `downstream` columns in `data/profile_inserts/`.
- Run `make profile_inserts`.
- The script uses fuzzy matching to find probe sequences (forward and reverse complement).
- Outputs include extracted inserts in FASTA format, comprehensive QC plots, and detailed metrics.
- QC metrics include length distribution, GC content, sequence composition, probe performance, and duplicate analysis.

**Features:**
- Fuzzy matching for probe annealing (configurable threshold)
- Progress bars for long-running analyses
- Comprehensive QC plots (12 different visualizations)
- Detailed logging and error handling
- Support for multiple probe pairs per CSV

---

## Identifying, Counting, and Generating Consensus UMI-Gene Pairings

For long-read sequencing libraries tagged with UMIs:
- Place your `.fastq.gz` files in `data/umi_hunter/` along with a gene template (`template.fasta`) and a configuration file `umi_hunter.csv`.
- The script outputs UMI-gene clusters, their counts, and consensus genes for clusters with >10 representatives to `results/umi_hunter/`.
- By default, the script clusters at 90% UMI identity.

---

## Identifying and Counting Mutants from Long-Read Data Without UMIs

- Save your gene reference (coding sequence only) in `data/mutation_caller/mutation_caller_template.fasta`.
- Place your `.fastq.gz` file in `data/mutation_caller/` and run `make mutation_caller`.
- The script outputs mutation counts and co-occurrence data to `results/mutation_caller/`.

---

## Making Mutant Proteins with SLIM


- This is a quick tool to design primers for SLIM cloning, which can add mutations to specific spots in a protein with overnight ease.
- Add your gene template (coding sequence only) to `data/design_slim/slim_template_gene.fasta` and your whole plasmid to `data/design_slim/slim_context.fasta`.
- Specify mutants in `data/design_slim/slim_target_mutations.csv` (first column: `mutations`).
- The script outputs designed primers to `results/design_slim/SLIM_primers.csv`.

**Mutation Nomenclature Examples:**
- Substitution: `A123G`
- Deletion: `T241Del`
- InDel (inter-codon): `T241InDelA242S`
- Insertion after codon: `T241TS` (insert Ser after Thr241)
- Codon replacement insertion: `L46GP` (replace Leu46 with Gly-Pro)

**Experimental:**
Total time:
~ 3h hands-on (not inc. protein purification), 72h DNA -> pure mutant protein

Contributions:
2x PCRs (~2h + 10 mins setup), SLIM thermocycling (50 mins), transformation (30 mins setup + overnight incubation), colony growth (5 mins hands-on, 12 h growth), DNA recovery (important for validation, 30 mins hands-on), protein expression and purification (2x overnight)

SLIM protocol:
Once you have the primers, run two normal PCRs using (A) long fwd + short rvs and (B) long rvs + short fwd. You can then add 10 ul of each PCR product to 10 ul of H-buffer, composed of 150 mM Tris pH 8, 400 mM NaCl and 60 mM EDTA. Incubate this (total volume 30 ul) in a thermocycler using the following protocol: 99 oC, 3:00 -> 2x [65 oC, 5:00 -> 30 oC, 15:00] -> Hold at 4 oC. You may then transform either NEB 5a or BL21 (DE3) with this mixture without further purification.

This code has been experimentally validated for simultaneous cloning of tens of mutants.
---

## Making Mutant Proteins with Gibson Assembly

- Add your gene template (coding sequence only) to `data/design_gibson/gibson_template_gene.fasta` and your whole plasmid to `data/design_gibson/gibson_context.fasta`.
- Specify mutants in `data/design_gibson/gibson_target_mutations.csv` (first column: `mutations`).
- The script outputs designed primers and an assembly plan to `results/design_gibson/`.
- Multi-mutants can be specified by adding a `+` between mutations in one cell of the `.csv` file.

**Mutation Nomenclature Examples:**
- Substitution: `A123G`
- Deletion: `T241Del`
- InDel (inter-codon): `T241InDelA242S`
- Insertion after codon: `T241TS` (insert Ser after Thr241)
- Codon replacement insertion: `L46GP` (replace Leu46 with Gly-Pro)

*Note: If target mutations are too close together and primer regions would overlap, run the reaction sequentially to incorporate multi-mutations.*

---

## ep-library-profile: Profiling a DNA Library Without UMI Dependence

- Place `.fastq.gz` files in `/data/ep-library-profile/`.
- Provide the mutational region of interest (`/data/ep-library-profile/region_of_interest.fasta`) and the whole plasmid (`/data/ep-library-profile/plasmid.fasta`).
- Run `make ep-library-profile`.
- Outputs include coverage, mutation rate, spectrum, and occurrence rates of specific mutations, all saved in `results/ep-library-profile/` along with a log file.

---

## Logging

All major scripts generate detailed log files in their respective `results/.../` directories. These logs capture all major steps, errors, and outputs, making it easier to debug and trace your runs.

---

## Coming Soon

- Primer design for KLD cloning of mutants
- Expansion of all primer design software for multi-mutations
- Expansion of mutation caller to handle indels


