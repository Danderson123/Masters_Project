#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:33:34 2020

@author: danielanderson
"""

#Panaroo returns errors when trying to generate a gml file from a single gff. 
#This script generates a single gml for use in the graph merger function. 

def single_network(gff_file):
    
    from prokka import process_prokka_input
    import os
    from cdhit import run_cdhit
    from generate_network import generate_network
    from merge_graphs import cluster_centroids, simple_merge_graphs
    import networkx as nx
    
    files = list(gff_file) #Supply gff filename as a list          
    
    output_dir = str(os.getcwd()) + '/' #output directory is current working directory  
    
    process_prokka_input(files,output_dir, '', 1) #run function on a single gff file with one core  
    
    # Cluster protein sequences using cdhit
    cd_hit_out = output_dir + "single_protein_cdhit_out.txt"
    run_cdhit(input_file=output_dir + "single_protein_CDS.fasta",
              output_file=cd_hit_out,
              id=0.98,
              s=0.98,
              quiet='',
              n_cpu=2)
    
    G_single, centroid_contexts_single, seqid_to_centroid_single = generate_network(
        cluster_file=cd_hit_out + ".clstr",
        data_file=output_dir + "single_gene_data.csv",
        prot_seq_file=output_dir + "single_protein_CDS.fasta",
        all_dna=False)
    
    for adj in G_single._adj:
        for x in G_single._adj[adj]:
            y = G_single._adj[adj][x]
        
            y.pop('members')
            zero = {'members': 0}
            y.update(zero)
        
            genomes = {'genomeIDs' : '0'}
            y.update(genomes)
        
    for node in G_single._node:
        y = G_single._node[node]
        y.pop('members')
        
        zero = {'members': 0}
        y.update(zero)
        
        rep = "[']"
        for char in rep:
            y['centroid'] = (str(y['centroid'])).replace(char, '')
            y['dna'] = (str(y['dna'])).replace(char, '')
            y['protein'] = (str(y['protein'])).replace(char, '')
        
        y['hasEnd'] = int(y['hasEnd'])
        y['mergedDNA'] = int(y['mergedDNA'])
        y['paralog'] = int(y['paralog'])
        
        y['longCentroidID'] = list(y['longCentroidID'])
        y['seqIDs'] = list(y['seqIDs'])
    
    nx.write_gml(G_single, output_dir + "single_graph.gml") #members are assigned as inhibitset[0]. needs to be 0
    
    return 
