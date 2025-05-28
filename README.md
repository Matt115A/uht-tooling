To use this repo:

chmod +x setup.sh
./setup.sh


---------- General structure ----------

This code is set up do that configs and data are found in data/.../, where ... depends on the script that you care to run. Likewise, results are saved to results/.../

The currently available commands are: 

	make nextera_primers
	make umi_hunter
	make design_slim
	make mutation_caller
	make mutation_rate_calculator

-----------------------------------------------------

---- Automated Nextera XT primer deșign software ----

One-PCR-to-flowcell workflow

Illumina will try to sell you expensive kits to use their sequencers, when this often isn't needed if you can order primers! To use this software, upload a .csv file in data/nextera_designer/ called nextera_designer.csv. The first column should be called binding_region, and the first two rows should the forward and reverse primers that you want to use to generate your illumina amplicon. Use normal primer design rules here, inputting the sequence in the 5'->3' order (same as if you were ordering primers). I have pre-loaded six i5 and six i7 indices here for your enjoyment, and every one will be used to give you a set of forward and reverse primers. You can use the primers in any combination to create up to 36 unique amplicons. 

Once you navigate to the root and setup your conda environment using the commands at the top of the readme, you can run make nextera_primers and a new .csv will appear in results/nextera_designer/, ready for ordering! Once you have run a PCR with these primers and purified the products, they can be loaded directly on to the flow-cell. Enjoy!

Pipeline for kit-free illumina sample prep:

	1. PCR using a i5/i7 primer pair, monitor a portion of the reaction using qPCR. Cap the cycle number so as to get only 10% of the final yield. This minimises PCR amplification bias.
	2. Purify the DNA product, removing primers. This is important as the primers will be in great molar excess due to their small size, and will flood the chip with their binding regions. Purification is easiest using SPRIselect beads from Beckman Coulter. You will have to empirically determine the bead:DNA volume ratio that removes primers while retaining amplicon for your sample - a good starting point is 0.65 parts bead volume to 1 part DNA solution.
	3. Verify primer removal using electrophoresis. This is best performed using a 2100 BioAnalyser and a DNA detection chip. Ensure that the molar ratio is satisfactory - you will loose read depth in a manner proportional to the primer impurities in the sample. 

-----------------------------------------------------

---- Identifying, counting and generating consensus UMI-gene pairings ----

If you send a library for long-read sequencing tagged with UMIs, you will get back a .fastq.gz file that can be complex to interpret. This script runs in a dictionary-free way, finding UMIs and genes based on user-provided flanking regions. Simply place your .fastq.gz files in data/umi_hunter/ along with a gene template for your coding sequence (template.fasta) and a configuration umi_hunter.csv. This file should have columns umi_flanks,gene_flanks,umi_min_max. These columns should be populated with the first row being the upstream flank and the second row being the downstream flank. The umi_min_max should be populated with the first row being minimum expected UMI length, and the second row being the maximum expected umi length. The script will output a set of files into results/umi_hunter/, where you will find the UMI-gene clusters, their counts, and the consensus genes for clusters with >10 representatives. It will also tell you the length difference of the consensus compared to the reference, which will allow you to identify any indels that may be present.

-----------------------------------------------------

---- Identifying and counting mutants from long-read data without UMIs ----

Sometimes, you won't have UMIs, but still want to call mutations from a pool of long-reads. This script allows you to do just that. You need to save your gene reference (coding sequence only) in data/mutation_caller/mutation_caller_template.fasta. Then, simply save your .fastq.gz file in data/mutation_caller, and run 'make mutation_caller'. It will show you a bar plot of the number of mutations at each position - use this to inform the threshold for a reliable call. Enter this in the command line, and then the list of single mutations will be saved to results/mutation_caller/singles, and single mutants that co-occur (i.e. are likely present as double mutants) will be saved in results/mutation_caller/coocuring_mutations.csv.

-----------------------------------------------------

---- Making mutant proteins with SLIM ----

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

This code has been experimentally validated for designing a set of 12 mutations into a WT sequence.

-----------------------------------------------------

---- Calculating the error rate of a library without UMI-dependance ----

If you have cloned a DNA library that contains mutations, you will naturally be interested in the rate at which they occur to ensure you have achieved a reasonable mutational load. This code allows you to assess this mutational rate using the output from a whole-plasmid sequencing run.

Place your .fastq.gz file into /data/mutation_rate_calculator/, and also place your region of interest (region_of_interest.fasta) and the whole plasmid that contains the region of interest (plasmid.fasta) in here. Then simply run make mutation_rate_calculator, and the mutation rate in your area of interest compared to the background rate will be saved in results/mutation_rate_calculator/, in the log file. It will also be displayed in the terminal, along with plots that show how the mutations are distributed across the region of interest.

-----------------------------------------------------

Coming soon: 

- Primer design for KLD cloning of mutants
- Detailed log files for all scripts
- Expansion of the primer design software for multi-mutations
- Expansion of mutation caller to handle indels


