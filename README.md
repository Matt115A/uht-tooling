To use this repo:

chmod +x setup.sh
./setup.sh

---- Automated Nextera XT primer deșign software ----

One-PCR-to-flowcell workflow

Illumina will try to sell you expensive kits to use their sequencers, when this often isn't needed if you can order primers! To use this software, upload a .csv file in data/ called nextera_designer.csv. The first column should be called binding_region, and the first two rows should the forward and reverse primers that you want to use to generate your illumina amplicon. Use normal primer design rules here, inputting the sequence in the 5'->3' order (same as if you were ordering primers). I have pre-loaded six i5 and six i7 indices here for your enjoyment, and every one will be used to give you a set of forward and reverse primers. You can use the primers in any combination to create up to 36 unique amplicons. 

Once you navigate to the root and setup your conda environment using the commands at the top of the readme, you can run make nextera_primers and a new .csv will appear in results, ready for ordering! Once you have run a PCR with these primers and purified the products, they can be loaded directly on to the flow-cell. Enjoy!

-----------------------------------------------------

---- Identifying, counting and generating consensus UMI-gene pairings ----

If you send a library for long-read sequencing tagged with UMIs, you will get back a .fastq.gz file that can be complex to interpret. This script runs in a dictionary-free way, finding UMIs and genes based on user-provided flanking regions. Simply place your .fastq.gz files in data/ along with a gene template for your coding sequence (template.fasta) and a configuration umi_hunter.csv. This file should have columns umi_flanks,gene_flanks,umi_min_max. These columns should be populated with entry one being the upstream flank and entry two being the downstream flank. The umi_min_max should be populated with the first entry being minimum expected UMI length, and the second entry being the maximum expected umi length. The script will output a set of files into results/, where you will find the UMI-gene clusters, their counts, and the consensus genes for clusters with >10 representatives. It will also tell you the length difference of the consensus compared to the reference, which will allow you to identify any indels that may be present.

-----------------------------------------------------

---- Identifying and counting mutants from long-read data without UMIs ----

Sometimes, you won't have UMIs, but still want to call mutations from a pool of long-reads. This script allows you to do just that. You need to save your gene reference (coding sequence only) in data/mutation_caller_template.fasta. Then, simply save your .fastq.gz file in /data, and run 'make mutation_caller'. It will show you a bar plot of the number of mutations at each position - use this to inform the threshold for a reliable call. Enter this in the command line, and then the list of single mutations will be saved to results/singles, and single mutants that co-occur (i.e. are likely present as double mutants) will be saved in results/coocuring_mutations.csv.

-----------------------------------------------------

---- Making mutant proteins with SLIM ----

This is a quick tool to design primers for SLIM cloning, which can add mutations to specific spots in a protein with overnight ease.

Simply add your gene template (coding sequence only) to data/slim_template_gene.fasta and your whole plasmid to data/slim_context.fasta. Ensure that your gene is contiguous in the context file and not split over the start and end. Then, specify mutants in the file data/slim_target_mutations.csv, naming the first column 'mutations' and then specifying mutations in the form 'M1P', or wild-type, position, mutant.

Experimental:

Total time:
~ 3h hands-on (not inc. protein purification), 72h DNA -> pure mutant protein

Contributions:
2x PCRs (~2h + 10 mins setup), SLIM thermocycling (50 mins), transformation (30 mins setup + overnight incubation), colony growth (5 mins hands-on, 12 h growth), DNA recovery (important for validation, 30 mins hands-on), protein expression and purification (2x overnight)

SLIM protocol:
Once you have the primers, run two normal PCRs using (A) long fwd + short rvs and (B) long rvs + short fwd. You can then add 10 ul of each PCR product to 10 ul of H-buffer, composed of 150 mM Tris pH 8, 400 mM NaCl and 60 mM EDTA. Incubate this (total volume 30 ul) in a thermocycler using the following protocol: 99 oC, 3:00 -> 2x [65 oC, 5:00 -> 30 oC, 15:00] -> Hold at 4 oC. You may then transform either NEB 5a or BL21 (DE3) with this mixture without further purification.

This code has been experimentally validated for designing a set of 12 mutations into a WT sequence.

-----------------------------------------------------