#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script generates crude annotations from an input fasta file with no pre-existing annotations
"""


import glob 
from tqdm import tqdm 
import re
import os
import time
import numpy as np
import pandas as pd
from multiprocessing import Pool

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])
    
    
def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3] #Defining a codon as 3 bases
            protein+= table[codon] #Translate the codon into an amino acid based on the dictionary and append this to the protein sequence
    return protein

def get_annotations(title, sequence):

    start1 = r"ATG"
    start2 = r"TTG"
    start3 = r"GTG"
    
    start4 = r"CAT"
    start5 = r"CAA"
    start6 = r"CAC"
    
    stop1 = r"TAG"
    stop2 = r"TAA"
    stop3 = r"TGA"
    
    stop4 = r"CTA"
    stop5 = r"TTA"
    stop6 = r"TCA"
    
    startList = [start1, start2, start3, stop4, stop5, stop6]
    start_indexes= []
    start_strand = []
    
    print("Identifying start and stop codons...")
    for start_regex in startList:
        if re.search(start_regex, sequence):
            stop_codon_list = re.finditer(start_regex,sequence)
            for y in stop_codon_list:
                start_indexes.append(y.start())
                if start_regex == "ATG" or start_regex == "TTG" or start_regex == "GTG":
                    start_strand.append("+")
                else:
                    start_strand.append("-")

    stopList = [stop1, stop2, stop3, start4, start5, start6]
    stop_indexes = []
    stop_strand = []

    for stop_regex in stopList:
        if re.search(stop_regex, sequence):
            stop_codon_list = re.finditer(stop_regex,sequence)
            for y in stop_codon_list:
                stop_indexes.append(y.end())
                if stop_regex == "TAG" or stop_regex == "TAA" or stop_regex == "TGA":
                    stop_strand.append("+")
                else:
                    stop_strand.append("-")
    
    start_indexes_forward = []
    stop_indexes_reverse = []
    for x in range(len(start_strand)):
        if start_strand[x] == "+":
            start_indexes_forward.append(start_indexes[x] + 1)
        else:
            stop_indexes_reverse.append(start_indexes[x] + 1)
    
    stop_indexes_forward = []
    start_indexes_reverse = []
    for x in range(len(stop_strand)):
        if stop_strand[x] == "+":
            stop_indexes_forward.append(stop_indexes[x])
        else:
            start_indexes_reverse.append(stop_indexes[x])
    
    start_indexes_forward = sorted(start_indexes_forward)
    stop_indexes_forward = sorted(stop_indexes_forward)
    start_forward_array = np.array(start_indexes_forward)
    stop_forward_array = np.array(stop_indexes_forward)
    
    print("Calculating distances...")
    
    final_forward_list = []
    for start_pos_forward in tqdm(start_forward_array):
        distances_forward = stop_forward_array - start_pos_forward
        distances_forward = distances_forward + 1
        middle_forward = np.where(distances_forward <= 1, 2, distances_forward)
        final_forward = list(np.where((middle_forward > 2) & (middle_forward % 3 == 0), 1, middle_forward))
        if 1 in final_forward:
            final_forward_list.append(final_forward.index(1))
        else:
            final_forward_list.append("")
        
    if len(final_forward_list) > 0:
         itemindex_forward = pd.DataFrame(final_forward_list, columns = ["index"])
    else:
        itemindex_forward = pd.DataFrame()
        
    to_remove_forward = []
    for x in range(len(itemindex_forward["index"])):
        if itemindex_forward["index"][x] == "":
            to_remove_forward.append(x)
    itemindex_forward = itemindex_forward.drop(itemindex_forward.index[to_remove_forward])
        
    locus_number = 0
    annotations = set()
    
    print("Generating annotations...")

    for x in itemindex_forward['index'].index:
        start_codon = start_indexes_forward[x]
        stop_codon = stop_indexes_forward[itemindex_forward['index'][x]]
        if not start_codon > stop_codon and stop_codon - start_codon >30:
            annotation = title
            annotation += "\tCRUDE\tCDS\t"
            annotation += str(start_codon)
            annotation += "\t"
            annotation += str(stop_codon)
            annotation += "\t.\t"
            annotation += "+"
            annotation += "\t0\tID="
            annotation += title
            annotation += "_"
            annotation += str(locus_number)
            annotation += ";product=theoretical_protein"
            annotations.add(annotation)
            locus_number += 1
     
    start_indexes_reverse = sorted(start_indexes_reverse)
    stop_indexes_reverse = sorted(stop_indexes_reverse)
    start_reverse_array = np.array(start_indexes_reverse)
    stop_reverse_array = np.array(stop_indexes_reverse)
    
    final_reverse_list = []
    for start_pos_reverse in tqdm(start_reverse_array):
        distances_reverse = stop_reverse_array - start_pos_reverse
        distances_reverse = distances_reverse - 1
        middle_reverse = np.where(distances_reverse >= -1, -2, distances_reverse)
        final_reverse = list(np.where((middle_reverse < -2) & (middle_reverse % 3 == 0), -1, middle_reverse))
        if -1 in final_reverse:
            final_reverse_list.append((len(final_reverse) - final_reverse[::-1].index(-1) - 1))
        else:
            final_reverse_list.append("")
    
    if len(final_reverse_list) > 0:
         itemindex_reverse = pd.DataFrame(final_reverse_list, columns = ["index"])
    else:
        itemindex_reverse = pd.DataFrame()
        
    to_remove_reverse = []
    for x in range(len(itemindex_reverse["index"])):
        if itemindex_reverse["index"][x] == "":
            to_remove_reverse.append(x)
    itemindex_reverse = itemindex_reverse.drop(itemindex_reverse.index[to_remove_reverse])
        
    for y in itemindex_reverse['index'].index:
        start_codon = start_indexes_reverse[y]
        stop_codon = stop_indexes_reverse[itemindex_reverse["index"][y]]
        if not start_codon < stop_codon and start_codon - stop_codon >30:
            annotation = title
            annotation += "\tCRUDE\tCDS\t"
            annotation += str(stop_codon)
            annotation += "\t"
            annotation += str(start_codon)
            annotation += "\t.\t"
            annotation += "-"
            annotation += "\t0\tID="
            annotation += title
            annotation += "_"
            annotation += str(locus_number)
            annotation += ";product=theoretical_protein"
            annotations.add(annotation)
            locus_number += 1 
            
            
    annotations = sorted(list(annotations), key=lambda x: int(x.split('\t')[3]))
    
    return annotations


def generate_annotations(fasta_files):
    
    total_genes = []
    for header in tqdm(fasta_files):
        print(header)
        with open(header) as f:
            fasta = f.read()
        
        split = fasta.split(">")[1:]
        
        region_titles = []
        annotation_all_regions = []
        sequences_all_regions = ["##FASTA"]
            
        for x in tqdm(range(len(split))):
            title = split[x].split(" ")[0]
            fasta_split = split[x].split("\n")
            fasta_header = ">" + fasta_split[0]
            sequence = "".join(fasta_split[1:])
            sequences_all_regions.append(fasta_header)
            sequences_all_regions.append(sequence)

            annotations = get_annotations(title, sequence)
                  
            region_titles += ["##sequence-region " + title + " 1 " + str(len(sequence))]
            annotation_all_regions += annotations
        
        print(header + " annotations complete")

        gff_file = "\n".join(list(region_titles) + list(annotation_all_regions) + list(sequences_all_regions))
        
        total_genes += annotation_all_regions
        outfile_name = os.path.basename(header).split(".fna")[0]
        
        print("writing file")

        outfile_gff = open("result3_annotated/" + outfile_name + ".gff", "w")
        outfile_gff.write(gff_file)
        outfile_gff.close()
        
        print("file written")
        
    return

start_time = time.time()
fnas = glob.glob("to_crude/*.fna")
generate_annotations(fnas)
end_time = time.time()
print(str(end_time - start_time) + " seconds")

#if __name__ == '__main__':
#    start_time = time.time()
#    fnas = glob.glob("result3_to_annotate/*.fna")
#    chunks = [fnas[i::6] for i in range(6)]
#    pool = Pool(processes=6)
#    total = pool.map(generate_annotations, chunks)
#    end_time = time.time()
#    print(str(end_time - start_time) + " seconds" )
#    print(str(total))
     


