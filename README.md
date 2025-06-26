# uht-tooling: automating the molecular biology accessory to ultra-high throughput screening

To use this repo:

You will need conda installed - follow these instructions;

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Then, clone the repo into a directory of your choice by navigating there via command line and running: 

	git clone https://github.com/Matt115A/uht-tooling.git

Navigate to /uht-tooling/ and run:

	chmod +x setup.sh

	./setup.sh

Please note, this setup.sh file is verified for Mac only. Results may vary on other OS.

##  General structure

This code is set up do that configs and data are found in data/.../, where ... depends on the script that you care to run. Likewise, results are saved to results/.../

The currently available commands are: 

	make nextera_primers
	make umi_hunter
	make design_slim
	make mutation_caller
	make ep-library-profile
	make design_gibson

-----------------------------------------------------

## Automated Nextera XT primer deșign software

One-PCR-to-flowcell workflow

Illumina will try to sell you expensive kits to use their sequencers, when this often isn't needed if you can order primers! To use this software, upload a .csv file in data/nextera_designer/ called nextera_designer.csv. The first column should be called binding_region, and the first two rows should the forward and reverse primers that you want to use to generate your illumina amplicon. Use normal primer design rules here, inputting the sequence in the 5'->3' order (same as if you were ordering primers). I have pre-loaded twelve i5 and twelve i7 indices here for your enjoyment, and every one will be used to give you a set of forward and reverse primers. You can use the primers in any combination to create up to 144 unique amplicons. 

Once you navigate to the root and setup your conda environment using the commands at the top of the readme, you can run make nextera_primers and a new .csv will appear in results/nextera_designer/, ready for ordering! Once you have run a PCR with these primers and purified the products, they can be loaded directly on to the flow-cell. Enjoy!

Pipeline for kit-free illumina sample prep:

	1. PCR using a i5/i7 primer pair, monitor a portion of the reaction using qPCR. Cap the cycle number so as to get only 10% of the final yield. This minimises PCR amplification bias.
	2. Purify the DNA product, removing primers. This is important as the primers will be in great molar excess due to their small size, and will flood the chip with their binding regions. Purification is easiest using SPRIselect beads from Beckman Coulter. You will have to empirically determine the bead:DNA volume ratio that removes primers while retaining amplicon for your sample - a good starting point is 0.65 parts bead volume to 1 part DNA solution.
	3. Verify primer removal using electrophoresis. This is best performed using a 2100 BioAnalyser and a DNA detection chip. Ensure that the molar ratio is satisfactory - you will loose read depth in a manner proportional to the primer impurities in the sample. 

-----------------------------------------------------

## Identifying, counting and generating consensus UMI-gene pairings

If you send a library for long-read sequencing tagged with UMIs, you will get back a .fastq.gz file that can be complex to interpret. This script runs in a dictionary-free way, finding UMIs and genes based on user-provided flanking regions. Simply place your .fastq.gz files in data/umi_hunter/ along with a gene template for your coding sequence (template.fasta) and a configuration umi_hunter.csv. This file should have columns umi_flanks,gene_flanks,umi_min_max. These columns should be populated with the first row being the upstream flank and the second row being the downstream flank. The umi_min_max should be populated with the first row being minimum expected UMI length, and the second row being the maximum expected umi length. The script will output a set of files into results/umi_hunter/, where you will find the UMI-gene clusters, their counts, and the consensus genes for clusters with >10 representatives. It will also tell you the length difference of the consensus compared to the reference, which will allow you to identify any indels that may be present. By default, the script clusters at 90% UMI identity.

-----------------------------------------------------

## Identifying and counting mutants from long-read data without UMIs

Sometimes, you won't have UMIs, but still want to call mutations from a pool of long-reads. This script allows you to do just that. You need to save your gene reference (coding sequence only) in data/mutation_caller/mutation_caller_template.fasta. Then, simply save your .fastq.gz file in data/mutation_caller, and run 'make mutation_caller'. It will show you a bar plot of the number of mutations at each position - use this to inform the threshold for a reliable call. Enter this in the command line, and then the list of single mutations will be saved to results/mutation_caller/singles, and single mutants that co-occur (i.e. are likely present as double mutants) will be saved in results/mutation_caller/coocuring_mutations.csv.

-----------------------------------------------------

## Making mutant proteins with SLIM

This is a quick tool to design primers for SLIM cloning, which can add mutations to specific spots in a protein with overnight ease.

Simply add your gene template (coding sequence only) to data/design_slim/slim_template_gene.fasta and your whole plasmid to data/design_slim/slim_context.fasta. Ensure that your gene is contiguous in the context file and not split over the start and end. Then, specify mutants in the file data/design_slim/slim_target_mutations.csv, naming the first column 'mutations' and then specifying mutations in the form 'M1P', or wild-type, position, mutant.

If you have an insertion, the form is e.g. G56GT (read: insertion of T after G). If the insertion is inside a codon (read: insertion of a triplet in G56 that results in TP) the nomclamenture is G56TP.

If you have a clean deletion, the form is e.g. G56Del (read: deletion of G56). If the deletion in inter-codon, (read: deletion of a triple between G56 and L57 resulting in P at these positions) the form is G56InDelL57P.

In summary, some examples:

            Substitution:        A123G
            Deletion:            T241Del
            InDel (inter-codon): T241InDelA242S
            Insertion after codon: T241TS (insert Ser after Thr241)
            Codon replacement insertion: L46GP (replace Leu46 with Gly-Pro)


Experimental:

Total time:
~ 3h hands-on (not inc. protein purification), 72h DNA -> pure mutant protein

Contributions:
2x PCRs (~2h + 10 mins setup), SLIM thermocycling (50 mins), transformation (30 mins setup + overnight incubation), colony growth (5 mins hands-on, 12 h growth), DNA recovery (important for validation, 30 mins hands-on), protein expression and purification (2x overnight)

SLIM protocol:
Once you have the primers, run two normal PCRs using (A) long fwd + short rvs and (B) long rvs + short fwd. You can then add 10 ul of each PCR product to 10 ul of H-buffer, composed of 150 mM Tris pH 8, 400 mM NaCl and 60 mM EDTA. Incubate this (total volume 30 ul) in a thermocycler using the following protocol: 99 oC, 3:00 -> 2x [65 oC, 5:00 -> 30 oC, 15:00] -> Hold at 4 oC. You may then transform either NEB 5a or BL21 (DE3) with this mixture without further purification.

This code has been experimentally validated for designing 10s of mutations into a WT sequence.

-----------------------------------------------------

##  Making mutant proteins with Gibson assembly

This is a quick tool to design primers for Gibson cloning, which can add mutations to specific spots in a protein with overnight ease.

Simply add your gene template (coding sequence only) to data/design_gibson/gibson_template_gene.fasta and your whole plasmid to data/design_gibson/gibson_context.fasta. Then, specify mutants in the file data/design_gibson/gibson_target_mutations.csv, naming the first column 'mutations' and then specifying mutations in the form 'M1P', or wild-type, position, mutant.

If you have an insertion, the form is e.g. G56GT (read: insertion of T after G). If the insertion is inside a codon (read: insertion of a triplet in G56 that results in TP) the nomclamenture is G56TP.

If you have a clean deletion, the form is e.g. G56Del (read: deletion of G56). If the deletion in inter-codon, (read: deletion of a triple between G56 and L57 resulting in P at these positions) the form is G56InDelL57P.

In summary, some examples:

            Substitution:        A123G
            Deletion:            T241Del
            InDel (inter-codon): T241InDelA242S
            Insertion after codon: T241TS (insert Ser after Thr241)
            Codon replacement insertion: L46GP (replace Leu46 with Gly-Pro)

You can specify a multi-mutant by adding a '+' between mutations, in one cell of the .csv file. You will get an accompanying csv in the output that describes the order in which to use the primers. It will also suggest the melting temperatures to use for each PCR, as well as the length of the PCR amplicon, to inform the elongation time.

Please be aware that the code currently fails if the target mutations are too close together, and the primer regions would overlap. In this case, you should run the reaction sequentially to incorporate multi-mutations.

-----------------------------------------------------

## ep-library-profile: Profiling a DNA library without UMI-dependance 

If you have cloned a DNA library that contains mutations, you will be interested in:

- Cloning efficiency
- Mutational rate
- How mutational rate transfers to amino acid mutational rate
- Mutational spectrum
- Over-representation of specific mutations

Using this code, you can get all of this information from a single whole-plasmid sequencing run. You must provide: 

- .fastq.gz file(s), placed into /data/ep-library-profile/
- mutational region of interest, /data/ep-library-profile/region_of_interest.fasta
- whole plasmid containing the region of interest /data/ep-library-profile/plasmid.fasta

Now run:

make ep-library-profile

You will receive:
- The coverage of each region of the plasmid in the sequencing run, to assess whether the plasmid has been evenly sampled (this may be used to assess library cloning efficiency)
 - - This information is saved as a .csv file with coverage at each position calculated as the number of reads mapping to that position
 - - This information is also saved to the panel .png and .pdf files, as 'Full Plasmid Coverage with ROI (region of interest) Shaded'
- The mutation rate in your area of interest normalised by the background rate (calculated based on the plasmid) in basepairs
 - - This information is found in the summary.txt file, saved in many ways and forms, which should be self-evident from the names.
 - - It is also visualised in the summary plots as 'Mismatch Rate per Position: Gene of Interest'
- The estimated rate of amino acid mutations in the gene of interest, modelled under a Poisson distribution
 - - This information is saved in the summary.txt file, and also shown as the KDE plot in the summary plot - this plot estimates the odds of obtaining a mutant of a particular mutation order given a random pick from the library.
 - - The amino acid mutation rate is based on the simulation of 1000 mutation events of your specified gene of interest, under the calculated Poisson distribution of amino acid mutations from your data.
- The mutational spectrum is calculated from above-background mutational rate positions
 - - This information can be found in the summary.txt file and also the _mutaion_spectrum.pdf file. The mutational spectrum reported in the literature for the commonly-used Mutazyme II kit is also reported as a reference.
- The occurrence rates of specific mutations at sites with above-background mutational rates. If you have over-sampled your library, you can use this to estimate the relative abundance of mutations in each library.
 - - This information is saved in gene_base_distribution.csv and gene_aa_substitution_rates.csv, for base mutations and amino acid substitutions respectively. The rate is saved as odds of occurrence per gene.

These will be saved in results/ep-library-profile/, along with the log file. Note that the amino acid mutation rate is based on the simulation of 1000 mutation events, with the calculated Poisson distribution of amino acid mutations.

-----------------------------------------------------

Coming soon: 

- Primer design for KLD cloning of mutants
- Detailed log files for all scripts
- Expansion of all primer design software for multi-mutations
- Expansion of mutation caller to handle indels


