#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read through graph, 
find genes with 'reforund',
extract the corresponding isolate and gene name, 
search for the gene name in the alignments, 
extract the isolate's sequence,
open isolate GFF and organise as dataframe,
add refound to the end of the dataframe,
reorganise by the start and end,
write output'

going to have refound dictionary generation out of loop and then search for each isolate indepepndently to speed up computation time  

SORT OUT THE REVERSE SEAQRCHING. THE POSITON IS ALWAYS RELATIVE TO THE SENSE
"""
import networkx as nx
import time 
import pandas as pd
import glob
import os 
import sys
from tqdm import tqdm

def get_options(): #options for downloading and cleaning
 
    import argparse

    description = 'Use a panaroo-generated graph and mafft alignments to extract updated gffs'
    parser = argparse.ArgumentParser(description=description,
                                  prog='panaroo-update')

    io_opts = parser.add_argument_group('isolate')

    io_opts.add_argument("-d",
                     "--gff_directory",
                     dest="gff_dir",
                     required=True,
                     help='directory for gffs to be updated',
                     type=str)
    
    io_opts.add_argument("-a",
                        "--alignement_dir",
                        dest="aln_dir",
                        required=True,
                        help='directory of panaroo-outputted mafft alignements',
                        type=str)
                        
    io_opts.add_argument("-g",
                     "--graph",
                     dest="graph_dir",
                     required=True,
                     help="the directory where the graph file is located",
                     type=str)
                     
                     
    io_opts.add_argument("--output",
                     dest="output_dir",
                     required=False,
                     help="output directory for updated gff files, default is input directory",
                     type=str,
                     default=" ")
                    
    args = parser.parse_args()

    return (args)

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def position_finder(sequence, fasta_information_forward, fasta_information_reverse):
    sequence = sequence.upper()
    length = len(sequence)
    if fasta_information_forward.count(sequence) == 1:
        start_position = fasta_information_forward.find(sequence) + 1
        end_position = start_position + length - 1
        strand = '+'
    elif fasta_information_reverse.count(sequence) == 1:
        end_position = len(fasta_information_reverse) - int(fasta_information_reverse.find(sequence))
        start_position = end_position - length + 1
        strand = '-'
    else:
        start_position = 0
        end_position = 0
        strand = ''
    return start_position, end_position, strand

def line(start, end, strand, ID, name, sequence_region, file):
    gff_line = sequence_region + '\t' + 'Panaroo' + '\t' + 'CDS' + '\t' + str(start) + '\t' + str(end) + '\t' + '.' + '\t' + str(strand) + '\t' + '.' + '\t' + ID + ';' + 'Name=' + name + ';' + 'gbkey=Gene' + ';' + 'gene_biotype=protein_coding' #+ ';' + 'locus_tag='
    file.append(gff_line)
    return    
    
def generate_library(graph_path, alignement_path):
    
    G = nx.read_gml(graph_path + "final_graph.gml")
    
    refound_genes = []    
    for node in G._node:
        y = G._node[node]
        if 'refound' in y["geneIDs"]:
            refound_genes.append(y["name"])
    
    library = {}
    for refound_gene in tqdm(range(len(refound_genes))):
        with open(alignement_path + refound_genes[refound_gene] + '.aln.fas') as file:
            gene = (file.read().replace("-","")).split(">")[1:]
            isolates = []
            cluster = []
            sequence = []
            for x in range(len(gene)):
                isolate_split = gene[x].split(";")
                isolates.append(isolate_split[0])
                cluster_split = isolate_split[1].split("\n")
                cluster.append(cluster_split[0])
                sequence.append(''.join(cluster_split[1:]))
            dictionary = {"Isolates" : isolates, "clusters": cluster, "Sequences" : sequence}
        library[refound_genes[refound_gene]] = dictionary
    
    return library

def update_gff(isolate, input_gffs, library, output_dir):
        
    isolate_genes = []
    isolate_gene_sequence = []
    
    for key, value in library.items():
        for x in range(len(value['Isolates'])):
            if value['Isolates'][x] == isolate and "refound" in value['clusters'][x]:
                isolate_genes.append(key)
                isolate_gene_sequence.append(value['Sequences'][x])
    
    filename = isolate +'.gff'
    with open(input_gffs + '/' + filename) as gff:
        to_update = gff.read()

    split = to_update.split("##FASTA")
    
    fasta_information_forward = "".join(split[1].splitlines()[2:])
    fasta_information_reverse = reverse_complement(fasta_information_forward)
    gff_information = split[0].splitlines()
        
    index_to_remove = 0
    for split_line in range(len(gff_information)):
        if '##' in gff_information[split_line]:
            index_to_remove = split_line
    title = gff_information[0:index_to_remove + 1]           
    gff_information = gff_information[index_to_remove + 1:]
    
    for title_line in range(len(title)):
       if "sequence-region" in title[title_line]:
           sequence_region = title[title_line].split(" ")[1]
        
    refound_DataFrame = pd.DataFrame(isolate_genes, columns = ['name'])
    refound_DataFrame['sequence'] = isolate_gene_sequence
    refound_DataFrame['ID'] = 'ID='
      
    #refound_DataFrame['start'] , refound_DataFrame['end'], refound_DataFrame['strand'] = refound_DataFrame.apply(lambda row: position_finder(row['sequence'], fasta_information_forward, fasta_information_reverse), axis=1 , result_type="expand")
    
    refound_DataFrame['searched'] = refound_DataFrame['sequence'].apply(position_finder, args=[fasta_information_forward, fasta_information_reverse])
    
    start = []
    end = []
    strand = []
    
    for search in refound_DataFrame['searched']:
        start.append(search[0])
        end.append(search[1])
        strand.append(search[2])
        
    refound_DataFrame['start'] = start
    refound_DataFrame['end'] = end
    refound_DataFrame['strand'] = strand
    
    refound_DataFrame.apply(lambda row: line(row["start"], row["end"], row['strand'], row['ID'], row['name'], sequence_region, gff_information), axis=1)
    
    gff_information = sorted(gff_information, key=lambda x: x.split('\t')[3])
    
    with open(output_dir + '/' + filename, 'w') as f:
        for item in gff_information:
            f.write("%s\n" % item)
                
    return


def main():
    
    start = time.time()
    
    args = get_options()
    
    args.graph_dir = os.path.join(args.graph_dir, "")
    args.aln_dir = os.path.join(args.aln_dir, "")
    
    if args.output_dir == " ":
        args.output_dir = args.gff_dir
        
    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    
    print("Generating library...")
    
    library = generate_library(args.graph_dir, args.aln_dir)
    
    print("Library generated")
    
    paths_gff = args.gff_dir + "/*.gff"

    isolate_files = glob.glob(paths_gff)
    
    print("Updating annotations...")

    for isolate in tqdm(isolate_files):
        isolate = os.path.basename(isolate).split('.gff')[0]
        try:
            update_gff(isolate, args.gff_dir, library, args.output_dir)
        except:
            print("Error with {}.gff".format(isolate))
            pass
    
    end = time.time()
    print(end - start)
    
    sys.exit(0)
