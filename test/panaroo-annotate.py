#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will extract updated annotations. 
First I am going to see how long it takes to extract the annotations for a specific isolate from the graphs. 
Then I will change this so that python stores and reads through all of the generated alignment files to generate the annotation. 
The easy part is extracting the annotations. 
The hard part will be updating the alignments and then extracting the updated annotations.
If extracting the annotations doesn't take long, then just recreate the gff file every time.'

So now I am loading all of the alignments into one list, then i will generate the gff from searching for the specific isolate.

This will generate a gffs from the alignments but then need to figure out how you can update the alignments without having to rerun panaroo when you want to incorporate a new genome.

So panaroo doesn't store the positions of the genes. Need to search the genome for the genes. This will require the fasta file as input but this is needed for the gff file anyway.'
This can be a fasta downloaded automatically or an input fasta. Now i am just placing the fasta in the same directory.

retrieivng information for the fasta files is fine for the majority of the genes. The graphs can be used too get the annotations and all other information except for the specific sequences and positions of genes in the genomes. 
Gene searching is fine for the majority of the genes using the find function but those with any gaps in the alignements are much harder to search for. One way to solve this would be to search for the longest substrings and if the count is larger than n1 then add another string to the search. 
This would likely work for almost all but there will certainly be exceptions at some point or another. 
I am going to implement a variable gap length searching algorithm from https://pdf.sciencedirectassets.com/271538/1-s2.0-S0304397512X00252/1-s2.0-S0304397512002915/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjENj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQD4nwXwkayawwCofgI8qTSvLqlrd1G2brnyqTo1YqADsAIgBNpj9vN3cca5JzAH%2BzkYB87gsrWVuz0jcw4iBEHTPhAqtAMIQRADGgwwNTkwMDM1NDY4NjUiDHQNNt%2BsAFRYmtqdJCqRA0m8Hu1EUOvEKeEvWik3wE79%2FatY52msUJzTEXPqAORwkqs4JjiuZ91vzFm0aHZu0erPyrPrcPbFWcxrAqXGv9jDj35TygOAN%2BEgIyQE6TsmU1f65mr2w%2FEQHHrVt9K1ZI6tjmJhyfEd9RwCu7oECgxbc45uRq0PPp8fhpkUNb2pYWLESRxtO%2F8aesw%2Fb2CPelaiBN0vO1LEzOk3p5IGFY248%2BICdAz5wQ88PrjvbbpAPFdo831%2FqXr1xRz8vryyGCNByRFJoG3V5%2Fnz2BJW1fmYk1%2BcAD95TDfN7V8QrVY0X1GjrFVQHmzMnCBeW%2B9SDPgcVYercyFJx6G08rtR0OrsUYwNDMAomg5TMP85MIBh5Y3l8NvuBdpSmqPR2Du04nP%2BhS3VnGyis7TGP%2FPlncxmbhmLWqwuncZY3x7CtGyv9d0DrHL4MedcJu8Qf%2Bq7%2Ba6LAMIIm%2FSaQMFn%2B3zlyPAilfQyZlBUNOFbiRtF5Gb2WC74Ox0hcYI9CzudysjBY%2BIIf8AoFmmv4eciqj8oWST0MP3h0vYFOusBxhepkFvZCLtwlvgiq1gi5i4rVg0O50UU1GxgBsbPcRwHer3AW5zo2cjOMdTiJ88NEwEZIkdNhaycyZctkQjkgARgfAn%2BK2WghEFlgmh90AJevx%2BceIkV78kaYsdbkQzEVv4XvCBU%2BF%2BH4nFbe77J%2BVDmZfwqQSeCmBSUJYgWJU3QfIZY%2BLw%2B4ah%2FbRRqYsRj2jVXAGGOovvJsEZ9HvgsDKf1hyifHiji0NODfnje5NVJW0rQR%2BA5FLXNFsmkFG5qAvyR1Ux4PhLJOk9NnTjcon7M7w6QpC1Z04slXs9srp%2Fx5WYsa6lVHUC0VA%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20200601T082305Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYYBEQ3ZX6%2F20200601%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=8184652b7b82f306fee5db605bf39acb9ec1358e015a145788dfcdaa917b35e9&hash=194762fa19cc0315910386c6dd99b1f253d26ac6393d4b7da30830911af72ea4&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0304397512002915&tid=spdf-d57955ef-d4ba-446f-a40a-421554af2f2e&sid=7d0c1bf12fb5734f984835d10ecf7b547a60gxrqb&type=client
This will handle all possilbe exceptions and hopefullyincrease computation time.
"""



import networkx as nx
import re
import pandas as pd
import numpy as np

G = nx.read_gml("final_graph.gml")

isolate = '452723578'


graphs = []
for x in G._node:
    if isinstance(G._node[x]['members'], list):
        if 0 in G._node[x]['members']:
            annotation = ">" + G._node[x]['name'] + "; " + G._node[x]['description'] + "; " + G._node[x]['dna']
            graphs.append(annotation)
    elif isinstance(G._node[x]['members'], int):
        if G._node[x]['members'] == 0:
            annotation = ">" + G._node[x]['name'] + "; " + G._node[x]['description'] + "; " + G._node[x]['dna']
            graphs.append(annotation)
            
with open('452723578.fasta') as fasta:
    isolate_genome = fasta.read()
    
isolate_genome = isolate_genome.split('\n', 1)
title = isolate_genome[0]
isolate_genome = (isolate_genome[1]).replace('\n', '')
    
gene_name = []
isolates = []
clusters = []
genes = []

#Import and process alignments 
import glob   
path_fas = 'aligned_gene_sequences/*.aln.fas'   
path_fasta = 'aligned_gene_sequences/*.fasta'   

files_aln_fas = glob.glob(path_fas)
files_fasta = glob.glob(path_fasta)
files = files_aln_fas + files_fasta

for x in range(len(files)):
    gene_name.append(re.split('; |, |\/|\.',files[x])[1])

regex = re.compile('[^ACTG-]')

for file in files:
    isolates_single = []
    genes_single = []
    cluster_single = []
    f=open(file, 'r') 
    stripped = (f.read()).split('>')
    f.close() 
    del stripped[0]
    for x in range(len(stripped)):
        isolated = stripped[x].split(';')
        isolates_single.append(isolated[0])
        cluster_single.append((isolated[1]).split('\n')[0])
        seq_to_clean = (isolated[1]).replace('\n', '').upper()
        seq_cleaned = regex.sub('', seq_to_clean)
        genes_single.append(seq_cleaned)
    isolates.append(isolates_single)
    clusters.append(cluster_single)
    genes.append(genes_single)
        
dictonary = {'Names': gene_name, 'Isolates' : isolates, 'Clusters': clusters, 'Sequence': genes}

alignments = pd.DataFrame(gene_name, columns= ['Names'])
alignments['Isolates'] = isolates
alignments['Clusters'] = clusters
alignments['Gene'] = genes

#Indexes for genes and clusters of specific isolates.
index = []
sub_index =[]

for x in range(len(isolates)):
    for y in range(len(isolates[x])):
        if isolate in isolates[x][y]:
            index.append(int(x))
            sub_index.append(int(y))

query = pd.DataFrame(index, columns = ["Gene Index"])
query['Sub Index'] = sub_index

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

genome_reverse = reverse_complement(isolate_genome)


def position_finder(sequence):
    length = len(sequence)
    start_position = isolate_genome.find(sequence) + 1
    end_position = start_position + length - 1
    if start_position == 0:
        end_position = len(genome_reverse) - int(genome_reverse.find(sequence)) + 2
        start_position = end_position - length - 1
        strand.append('-')
    else:
        strand.append('+')
    return start_position, end_position

genes_name_list = []
cluster_list = []
sequence_list = []
start_position_list = []
end_position_list = []
strand = []

def finder(index, subindex):
    
    genes_name_list.append(gene_name[index])
    
    clus = clusters[index][subindex]
    #cluster = clus.split('_')
    #if cluster[1] == 'refound':
       #cluster = cluster[1] + '_' + cluster[2]
    #else:
        #cluster = int(cluster[2])
    cluster_list.append(clus)
    
    sequence = genes[index][subindex]
    sequence_list.append(sequence)
    
    start, end = position_finder(sequence)
    start_position_list.append(start)    
    end_position_list.append(end) 
    return

query.apply(lambda row: finder(row["Gene Index"], row["Sub Index"]), axis=1)

query["Gene Name"] = genes_name_list
query["Cluster"] = cluster_list
query["Sequence"] = sequence_list
query["Start Position"] = start_position_list
query["End Position"] = end_position_list
query['Strand'] = strand

for x in G._node:
    if 'peb1A' == G._node[x]['annotation']:
        print(G._node[x])
        comb = G._node[x]['dna']

split = comb.split(';')
for x in split:
   print(len(x))
   
   
genome_reverse.find(split[2])

sequence = "ATGA------------------------------------------------------------------------------------CAGAAAC----------------------CTTGATAAAAATTGAAAATTTACATAAATCCTTTGG--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AAAGAATG------------------------AAGTATTGAAGGGCATCAACCTCGAGATTAAAAGAGGAGAAGTTGTCGTTATCA------------TCGGTCCTTCAGGGAGCGGGAAATCTACCTTGCTTCGCTCTATGAATTTGTTGGAAGAAGCAACCAAGGGGAAGGTTATCTTTGAGGGAGTCGATATTACGGACAAGAAGAATGACCTGTTTGCCATGCGTGAGAAGATGGGCATGGTTTTTCAACAATTCAATCTCTTTCCTAATATGACTGTGATGGAAAATATCACCTTGTCCCCTATCAAGACCAAAGGTGACAGTAAGGCCGTTGCAGAGAAAAGAGCTCAGGAACTTTTGGAAAAAGTTGGTTTGCCAGATAAGGCAGACGCTTATCCACAGAGTTTGTCAGGTGGCCAGCAACAGCGGATTGCCATCGCGCGTGGGTTGGCTATGGAACCAGATGTTTTGCTCTTTGACGAGCCAACTTCAGCCCTAGATCCTGAGATGGTTGGAGAAGTTCTGGCTGTTATGCAAGATCTAGCCAAGTCAGGAATGACCATGGTTATCGTAACACATGAGATGGGATTTGCCCGTGAGGTGGCAGATCGTGTCATCTTTATGGCAGACGGTGTGGTTGTTGAAGACGGAACACCTGAGCAGATTTTTGAACAAACCCAAGGACAAAGGACTAAAGACTTCTTGAGTAAGGTTTTATAA"
sep = sequence.split('-')

sequence_strings = []

for x in range(len(sep)):
    if not sep[x] == '':
        sequence_strings.append(sep[x])

all_positions = list(sequence)
base_positions = []
for x in range(len(all_positions)):
    if not all_positions[x] == '-':
        base_positions.append(x)
    else:
        pass
        
start_positions = []
end_positions = []

start_index = 0
middle_index = 1
next_index = 2  

for x in range(len(base_positions)): 
    if x == 0:
        start_positions.append(base_positions[x])
    elif x + 1 == len(base_positions):
        end_positions.append(base_positions[x])
    else:
        start = base_positions[start_index]
        middle = base_positions[middle_index]
        next_middle = base_positions[x]
            
        if next_middle - middle == 1:
            middle_index += 1
            
        elif next_middle - middle > 1:
            start_index = x
            middle_index += 1
            start_positions.append(next_middle)
            end_positions.append(middle)
        else:
            pass
        
sequence_strings = pd.DataFrame(sequence_strings, columns = ["Sub Sequences"])
sequence_strings['Start'] = start_positions 
sequence_strings['End'] = end_positions
sequence_strings['Search length'] = len(all_positions)

for x in range(len(sequence_strings)):
    if len(sequence_strings) == 1:
        start_pos = isolate_genome.find(sequence_strings["Sub Sequences"][x]) + 1
        if start_pos == 0:
            start_pos = genome_reverse.find(sequence_strings["Sub Sequences"][x]) + 1   
        start_pos = start_pos - sequence_strings['Start'][x]
        end_pos = start_pos + sequence_strings['Search length'][x]
                
        
        
    pos = isolate_genome.find(sequence_strings["Sub Sequences"][x])
    if pos == -1:
        pos = genome_reverse.find(sequence_strings["Sub Sequences"][x])
