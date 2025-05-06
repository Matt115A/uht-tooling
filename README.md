To use this repo:

conda create --name uht-tooling python=3.10 ; 
conda activate uht-tooling ; 
pip install -r requirements.txt ; 

---- Automated Nextera XT primer deșign software ----

One-PCR-to-flowcell workflow

Illumina will try to sell you expensive kits to use their sequencers, when this often isn't needed if you can order primers! To use this software, upload a .csv file in data/ called nextera_designer.csv. The first column should be called binding_region, and the first two rows should the forward and reverse primers that you want to use to generate your illumina amplicon. Use normal primer design rules here, inputting the sequence in the 5'->3' order (same as if you were ordering primers). I have pre-loaded six i5 and six i7 indices here for your enjoyment, and every one will be used to give you a set of forward and reverse primers. You can use the primers in any combination to create up to 36 unique amplicons. 

Once you navigate to the root and setup your conda environment using the commands at the top of the readme, you can run make nextera_primers and a new .csv will appear in results, ready for ordering! Once you have run a PCR with these primers and purified the products, they can be loaded directly on to the flow-cell. Enjoy!

-----------------------------------------------------