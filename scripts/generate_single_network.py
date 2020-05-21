#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:33:34 2020

@author: danielanderson
"""

from prokka import process_prokka_input
import os
from cdhit import run_cdhit
from generate_network import generate_network
from merge_graphs import cluster_centroids, simple_merge_graphs
import networkx as nx

files = ['15902044.gff'] #Supply gff filename as a list          

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

G_multi = nx.read_gml('results/final_graph.gml')
G = ['G_single','G_multi']

clusters = cluster_centroids(G, output_dir)
G_combined = simple_merge_graphs(G)
