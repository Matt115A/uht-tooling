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
