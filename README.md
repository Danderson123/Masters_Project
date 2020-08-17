# Masters Project
This repository houses the scripts for my MSc in Epidemiology project "A pan-genomic functional annotation framework integrating genomic epidemiology and functional genomics". 

# Motivation 
There are currently no tools capable of functionally annotating bacterial genomes with the consistency and efficiency needed for high throughput genomic epidemiology. 

# Methods 
* retrieve_sequences.py: retrieve genomic sequences in FASTA format and functional annotations in GFF format then quality control and re-format for direct input into Panaroo. This script also outputs an "all_annotations.csv" file that is later the source of functional annotations. 
* Panaroo: Used to construct the "reference" pan-genome that clusters input features (https://github.com/gtonkinhill/panaroo).
* update_descriptions.py: Transfer clustered gene names in the Panaroo output to "all_annotations.csv"
* update_annotations.py: Update functional annotation files for isolates already in the graph
* 
