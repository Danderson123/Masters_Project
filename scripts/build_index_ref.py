#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The constructs a searcheable index from representative sequences associated with Panaroo-outputted pan-genomic graphs
"""

import networkx as nx
import os 
from tqdm import tqdm
import cobs_index as cobs

def long_centroid(graph):
    """extract protein sequences associated with each node and add to a single list"""
    all_proteins = []
    G = nx.read_gml(graph)
    for node in G._node:
        y = G._node[node]
        seq = y["protein"].split(";")
        for x in seq:
            all_proteins.append(">" + y["name"] + "\n" + x)
            
    all_proteins = "\n".join(all_proteins)
    return all_proteins

all_proteins = long_centroid("ATCC_INTEGRATED/final_graph.gml")

os.mkdir("INDEX")
os.mkdir("FILES")

protein_file = all_proteins.split(">")[1:]

titles = []
cleaned = []
for x in protein_file:
    cleaned.append("".join(x.splitlines()[1:]))
    titles.append(x.splitlines()[0])

hundred_thousand_indexes = []
for x in range(len(protein_file)):
    hundred_thousand_indexes.append(x)
    
proteins = []
gene_names = []
for protein in hundred_thousand_indexes:
    info = cleaned[protein]
    proteins.append(info)
    gene_names.append(titles[protein])

#Write out indeividual files for each protein sequence
duplicate_set = set()
duplicate_count = 1
for prot_index in tqdm(range(len(hundred_thousand_indexes))):
    info = proteins[prot_index]
    if not gene_names[prot_index] in duplicate_set:
        thou = open("FILES/" + gene_names[prot_index] + ".txt", "w")
        thou.write(info)
        thou.close()
        duplicate_set.add(gene_names[prot_index])
    else:
        thou = open("FILES/" + gene_names[prot_index] + ";" + str(duplicate_count) + ".txt", "w")
        thou.write(info)
        thou.close()
        duplicate_count += 1

#Construct the index
params = cobs.CompactIndexParameters()
params.term_size = 10
params.clobber = True               # overwrite output and temporary files
params.false_positive_rate = 0.01    # higher false positive rate -> smaller index
cobs.compact_construct("FILES/", "INDEX/" + "10_index_index.cobs_compact", index_params=params)
