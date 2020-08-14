#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
if multiple names then add synonymn
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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'R': 'R', 'K':'K', 'Y':'Y'}
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

def gff_row(start, end, strand, ID, name, sequence_region):
    gff_line = sequence_region + '\t' + 'Panaroo' + '\t' + 'CDS' + '\t' + str(start) + '\t' + str(end) + '\t' + '.' + '\t' + str(strand) + '\t' + '.' + '\t' + ID + ';' + 'Name=' + name + ';' + 'gbkey=CDS' + ';' + 'gene=' + name + ';' + 'gene_biotype=protein_coding' #+ ';' + 'locus_tag='
    return gff_line
    
def generate_library(graph_path, alignement_path):
    
    G = nx.read_gml(graph_path + "final_graph.gml")

    refound_genes = []    
    for node in G._node:
        y = G._node[node]
        if 'refound' in y["geneIDs"]:
            refound_genes.append(y["name"])
        elif "~" in y["name"]:
            refound_genes.append(y["name"])
        elif "group" in y["name"]:
            refound_genes.append(y["name"])
    
    library = {}
    for refound_gene in tqdm(range(len(refound_genes))):
        try:
            with open(alignement_path + refound_genes[refound_gene] + '.aln.fas') as file:
                gene = (file.read().replace("-","")).split(">")[1:]
        except:
            with open(alignement_path + refound_genes[refound_gene] + '.fasta') as file:
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
    
    return library, G

def update_gff(isolate, input_gffs, library, output_dir):
        
    isolate_genes = []
    isolate_gene_sequence = []
    
    for key, value in library.items():
        for x in range(len(value['Isolates'])):
            if value['Isolates'][x] == isolate: #and "refound" in value['clusters'][x]:
                isolate_genes.append(key)
                isolate_gene_sequence.append(value['Sequences'][x])
    
    filename = isolate +'.gff'
    with open(input_gffs + '/' + filename) as gff:
        to_update = gff.read()

    split = to_update.split("##FASTA")
    split_fasta = (split[1].split(">"))[1:]
    lines = split[0].splitlines()
    
    titles = []
    for line in range(len(lines)):
        if "##sequence-region" in lines[line]:
            titles.append(lines[line])
            end_title_line = line

    annotations = lines[end_title_line + 1:]
        
    annotation_region = []
    annotation_indexes = []
    for annotation in range(len(annotations) - 1):
        if not (annotations[annotation].split("\t")[0]) == (annotations[annotation + 1].split("\t")[0]):
            annotation_region.append(annotations[annotation].split("\t")[0])
            annotation_indexes.append(annotation)
    
    annotation_indexes.append(len(annotations) - 1)
    annotations_split = []
    
    region_to_remove = []
    for region_title in range(len(titles)):
        if not titles[region_title].split(" ")[1] in annotation_region:
            region_to_remove.append(region_title)
    
    for index in sorted(region_to_remove, reverse=True):
        del titles[index]

    start_index = 0
    for index in range(len(annotation_indexes)):
        end = annotation_indexes[index] + 1
        item = slice(start_index, end)
        start_index = end
        annotations_split.append(annotations[item])
    
    contigs = []
    for title in tqdm(range(len(titles))):
        fasta_information_forward = "".join(split_fasta[title].split('\n')[1:])
        fasta_information_reverse = reverse_complement(fasta_information_forward)
        gff_information = annotations_split[title]
        
        sequence_region = titles[title].split(" ")[1]
        refound_DataFrame = pd.DataFrame(isolate_genes, columns = ['name'])
        refound_DataFrame['sequence'] = isolate_gene_sequence
        refound_DataFrame['ID'] = 'ID='
        
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
        
        for annotation_row, row in refound_DataFrame.iterrows():
            if not row['end'] == 0:
                gff_information.append(gff_row(row["start"], row["end"], row['strand'], row['ID'], row['name'], sequence_region))
            else:
                pass
        
        gff_information = sorted(gff_information, key=lambda x: int(x.split('\t')[3]))
        
        positions_out = []
        for line in gff_information:
            positions_out.append((str(line.split('\t')[3]) + "," + str(line.split('\t')[4])))
            
        to_remove = []
        for pos in range(len(positions_out)):
            if positions_out.count(positions_out[pos]) > 1 and not "Panaroo" in gff_information[pos]:
                to_remove.append(pos)
        
        for index in sorted(to_remove, reverse=True):
            del gff_information[index]
            
        gff_information = "\n".join(gff_information)
        contigs.append(titles[title] + "\n" + gff_information + "\n##FASTA" + "".join(split[1:]))
        
    with open(output_dir + '/' + filename, 'w') as f:
        for item in contigs:
            f.write("%s\n" % item)
                
    return

def gene_frequencies(graph, output_dir):
    gene_names = []
    frequencies = []
    
    for value in graph.graph.values():
        num_isolates = len(value)
                           
    for node in tqdm(graph._node):
        y = graph._node[node]
        gene_names.append(y["name"].lower())
        num_sequences = y["seqIDs"]
        unique = set()
        for x in range(len(num_sequences)):
            unique.add(num_sequences[x].split("_")[0])
        frequency = (len(unique)/num_isolates) * 100
        frequencies.append(frequency)
        
    gene_table = pd.DataFrame(gene_names, columns = ["Gene Name"])
    gene_table["Frequency (%)"] = frequencies
    
    gene_table = gene_table.sort_values(by='Gene Name')
    
    gene_table.to_csv(output_dir + "/"+ "Gene_frequencies.csv", index=False)

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
    
    library, G = generate_library(args.graph_dir, args.aln_dir)
    
    print("Library generated")
    
    paths_gff = args.gff_dir + "/*.gff"

    isolate_files = glob.glob(paths_gff)
    
    print("Updating annotations...")

    for isolate in tqdm(isolate_files):
        isolate = os.path.basename(isolate).split('.gff')[0]
        update_gff(isolate, args.gff_dir, library, args.output_dir)
        #except:
           # print("Error with {}.gff".format(isolate))
            #pass
     
    print("Calculating gene frequencies...")
    
    gene_frequencies(G, args.output_dir)
    
    end = time.time()
    print(end - start)
    
    sys.exit(0)

main()
