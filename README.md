# Overview
This repository houses the scripts for my MSc in Epidemiology project "A pan-genomic functional annotation framework integrating genomic epidemiology and functional genomics". 

# Motivation 
There are currently no tools capable of functionally annotating bacterial genomes with the consistency and efficiency needed for high throughput genomic epidemiology. 

# Methods 
* retrieve_sequences.py: Retrieve genomic sequences in FASTA format and functional annotations in GFF format using Biopython Entrez, then quality control and re-format for direct input into Panaroo. This script also outputs an "all_annotations.csv" file that is later the source of functional annotations. 
* retrieve_samples.py: Used for high throughput sequence retrieval as Biopython Entrez has a tendency to crash in these situtations. 
* Panaroo: Used to construct the "reference" pan-genome that clusters input features (https://github.com/gtonkinhill/panaroo).
* update_descriptions.py: Transfer clustered gene names in the Panaroo output to "all_annotations.csv"
* update_annotations.py: Update functional annotation files for isolates already in the graph
* integrate.py: Integrate the functional annotations for a single isolate into a Panaroo graph to incorporate newly characterised genomic information. 
* gen_index.py: Construct a COBS index from a Panaroo reference output (https://github.com/bingmann/cobs). 
* annotate.py: Functionally annotate a genome in fasta format using a Panaroo-updated "all_annotations.csv" and COBS-constructed index. Features are annotated by consensus between Panaroo and COBS. If they disagree, the names identified by both tools are included. 
