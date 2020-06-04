#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script extracts updated annotations using the panaroo graph and mafft alignment outputs.
"""

import networkx as nx
import re
import pandas as pd
import os
import glob


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def finder(index, subindex):
    genes_name_list.append(gene_name[index])

    clus = clusters[index][subindex]
    cluster_list.append(clus)

    sequence = genes[index][subindex]
    sequence = sequence.replace('-', '')
    sequence_list.append(sequence)

    start, end = position_finder(sequence)
    start_position_list.append(start)
    end_position_list.append(end)
    return

def position_finder(sequence):
    if not genome.count(sequence) == 0:
        length = len(sequence)
        start_position = isolate_genome.find(sequence) + 1
        end_position = start_position + length - 1
        if start_position == 0:
            end_position = len(genome_reverse) - int(genome_reverse.find(sequence)) + 2
            start_position = end_position - length - 1
            strand.append('-')
        else:
            strand.append('+')
    else:
        start_position = 'None'
        end_position = 'None'
        strand.append('None')
    return start_position, end_position

def line(start, end, strand, ID, name):
    gff_line = '' + '\t' + 'Panaroo' + '\t' + 'CDS' + '\t' + str(start) + '\t' + str(end) + '\t' + '.' + '\t' + str(strand) + '\t' + '.' + '\t' + ID + ';' + 'Name=' + name + ';' + 'gbkey=Gene' + ';' + 'gene_biotype=protein_coding' #+ ';' + 'locus_tag='
    gff_lines.append(gff_line)
    return

def get_options(): #options for downloading and cleaning
 
    import argparse

    description = 'Use a panaroo-generated graph and mafft alignments to extract an updated GFF file for a specific isolate'
    parser = argparse.ArgumentParser(description=description,
                                  prog='panaroo-annotate')

    io_opts = parser.add_argument_group('isolate')

    io_opts.add_argument("-i",
                     "--genome_id",
                     dest="genome_id",
                     required=True,
                     help='the genome ID for the isolate of interest',
                     type=str)
    
    io_opts.add_argument("-f",
                        "--genome_file",
                        dest="fasta_file",
                        required=True,
                        help='.fasta file of the isolate genome. The filename must be consistent with the isolate name.',
                        type=str)
                        
    io_opts.add_argument("-g",
                     "--graph",
                     dest="graph_dir",
                     required=True,
                     help="the directory where the graph file is located",
                     type=str)
                     
    io_opts.add_argument("-a",
                     "--alignements",
                     dest="aln_dir",
                     required=True,
                     help="the directory where the mafft alignements are located",
                     type=str)
                     
    io_opts.add_argument("-o",
                     "--output",
                     dest="output_dir",
                     required=True,
                     help="output directory for updated gff file",
                     type=str)
                    
    args = parser.parse_args()

    return (args)
    
def main():
    
    args = get_options()
    graph_location = args.graph_dir + "/" + "final_graph.gml"
    G = nx.read_gml(graph_location)

    isolate = args.isolate
    
    fasta_name = args.fasta_file + ".fasta"
    
    with open(fasta) as fasta:
        isolate_genome = fasta.read()
        
    isolate_genome = isolate_genome.split('\n', 1)
    title = isolate_genome[0]
    isolate_genome = (isolate_genome[1]).replace('\n', '')
        
    gene_name = []
    isolates = []
    clusters = []
    genes = []

    #Import and process alignments
    path_fas = args.aln_dir + "/*.aln.fas"
    path_fasta = aln_dir + "/*.fasta"

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

    gff = pd.DataFrame(index, columns = ["Gene Index"])
    gff['Sub Index'] = sub_index

    genome_reverse = reverse_complement(isolate_genome)

    genome = isolate_genome + genome_reverse

    genes_name_list = []
    cluster_list = []
    sequence_list = []
    start_position_list = []
    end_position_list = []
    strand = []
                
    gff.apply(lambda row: finder(row["Gene Index"], row["Sub Index"]), axis=1)

    gff["Name"] = genes_name_list
    gff["Cluster"] = cluster_list
    gff["Sequence"] = sequence_list
    gff["Start Position"] = start_position_list
    gff["End Position"] = end_position_list
    gff['Strand'] = strand
    gff['ID'] = 'gene-'

    gff_final = gff[['Start Position', 'End Position', 'Strand', 'ID', 'Name']]

    gff_final = gff_final.sort_values(by=['Start Position'])

    gff_lines = []

    gff_final.apply(lambda row: line(row["Start Position"], row["End Position"], row['Strand'], row['ID'], row['Name']), axis=1)

    filename = args.output_dir + "/" + isolate + ".gff"

    with open(filename, 'w') as f:
        for item in gff_lines:
            f.write("%s\n" % item)
    
    return
