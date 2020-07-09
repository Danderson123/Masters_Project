#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This cript generates crude annotations from an input fasta file if there are no pre-existing annotations. 
While it does not accurately predict genes, the hope is that it will run much faster than prokka and generate more annotations.
This will improve the gene finding potential of panaroo.

It will be tested on the first 10 samples from the Manhattan dataset. 

Information I want to record is:
    - number of annotations generated.
    - time to complete all 10 sequences. 
    - time to complete 1 sequence.
    - number of core genes identified.
    
start codons for now are ATG for + and CAT for - 
stop codons are TAG, TAA, and TGA for +
stop codons are CTA, TTA and TCA for - 
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

def get_forward_annotations(title, sequence, locus_number):
        
    start1 = r"ATG"
    start2 = r"TTG"
    start3 = r"GTG"

    stop1 = r"TAG"
    stop2 = r"TAA"
    stop3 = r"TGA"

    forwardstartList = [start1, start2, start3]
    start_indexes_forward = []

    for start_regex in forwardstartList:
        if re.search(start_regex, sequence):
            stop_codon_list = re.finditer(start_regex,sequence)
            for y in stop_codon_list:
                start_indexes_forward.append(y.start())

    stop_indexes_forward = []

    forwardstopList = [stop1, stop2, stop3]#, stop4, stop5, stop6, stop7, stop8, stop9]

    for stop_regex in forwardstopList:
        if re.search(stop_regex, sequence):
            stop_codon_list = re.finditer(stop_regex,sequence)
            for y in stop_codon_list:
                stop_indexes_forward.append(y.end() -1)
    
    start_indexes_forward = sorted(start_indexes_forward)
    stop_indexes_forward = sorted(stop_indexes_forward)

    start_array = np.array(start_indexes_forward)
    stop_array = np.array(stop_indexes_forward)

    distances = np.subtract(stop_array,start_array[:,None])
    distances1 = distances +1
    middle = np.where(distances1 <= 1, 2, distances1)

    final = np.where((middle > 2) & (middle % 3 == 0), 1, middle)
    
    minimum = np.min(final, axis=1)
        
    itemindex = pd.DataFrame(np.argmin(final,axis=1), columns = ["index"])
    
    to_remove = []

    for x in range(len(minimum)):
        if not minimum[x] == 1:
            to_remove.append(x)
    
    itemindex = itemindex.drop(itemindex.index[to_remove])
        
    annotations = set()
    for x in itemindex['index'].index:
        start_codon = start_indexes_forward[x]
        stop_codon = stop_indexes_forward[itemindex['index'][x]]
        if not start_codon > stop_codon and stop_codon - start_codon >50:
            annotation = title
            annotation += "\tCRUDE\tCDS\t"
            annotation += str(start_codon + 1)
            annotation += "\t"
            annotation += str(stop_codon + 1)
            annotation += "\t.\t"
            annotation += "+"
            annotation += "\t0\tID="
            annotation += title
            annotation += "_"
            annotation += str(locus_number)
            annotation += ";product=putative_protein_region"
            annotations.add(annotation)
            
            locus_number += 1
            
    annotations = sorted(list(annotations), key=lambda x: int(x.split('\t')[3]))
    
    return annotations, locus_number


def get_reverse_annotations(title, sequence, locus_number):

    sequence = reverse_complement(sequence)
    
    start1 = r"ATG"
    start2 = r"TTG"
    start3 = r"GTG"
    
    stop1 = r"TAG"
    stop2 = r"TAA"
    stop3 = r"TGA"
    
    forwardstartList = [start1, start2, start3]
    start_indexes_forward = []

    for start_regex in forwardstartList:
        if re.search(start_regex, sequence):
            stop_codon_list = re.finditer(start_regex,sequence)
            for y in stop_codon_list:
                start_indexes_forward.append(y.start())
    
    stop_indexes_forward = []
    
    forwardstopList = [stop1, stop2, stop3]#, stop4, stop5, stop6, stop7, stop8, stop9]
    
    for stop_regex in forwardstopList:
        if re.search(stop_regex, sequence):
            stop_codon_list = re.finditer(stop_regex,sequence)
            for y in stop_codon_list:
                stop_indexes_forward.append(y.end() -1)
            
    start_indexes_forward = sorted(start_indexes_forward)
    stop_indexes_forward = sorted(stop_indexes_forward)
    
    start_array = np.array(start_indexes_forward)
    stop_array = np.array(stop_indexes_forward)
    
    distances = np.subtract(stop_array,start_array[:,None])
    distances1 = distances +1
    middle = np.where(distances1 <= 1, 2, distances1)
    
    final = np.where((middle > 2) & (middle % 3 == 0), 1, middle)
    
    minimum = np.min(final, axis=1)
    itemindex = pd.DataFrame(np.argmin(final,axis=1), columns = ["index"])
    
    to_remove = []

    for x in range(len(minimum)):
        if not minimum[x] == 1:
            to_remove.append(x)
    
    itemindex = itemindex.drop(itemindex.index[to_remove])
        
    annotations = set()
    
    seq_length = len(sequence)
    for x in itemindex['index'].index:
        start_codon = seq_length - int(start_indexes_forward[x])
        stop_codon = seq_length - int(stop_indexes_forward[itemindex['index'][x]])
        if not start_codon < stop_codon and start_codon - stop_codon >50:
            annotation = title
            annotation += "\tCRUDE\tCDS\t"
            annotation += str(stop_codon)
            annotation += "\t"
            annotation += str(start_codon)
            annotation += "\t.\t"
            annotation += "-"
            annotation += "\t0\tID="
            annotation += title
            annotation += "_rv_"
            annotation += str(locus_number)
            annotation += ";product=putative_protein_region"
            annotations.add(annotation)
            
            locus_number += 1
            
    annotations = sorted(list(annotations), key=lambda x: int(x.split('\t')[3]))
        
    return annotations


def generate_annotations(fnas):
    
    total_genes = []
    for header in fnas:
        print(header)
        with open(header) as f:
            fasta = f.read()
        
        split = fasta.split(">")[1:]
        
        region_titles = []
        annotation_all_regions = []
        sequences_all_regions = ["##FASTA"]
        
        for x in range(len(split)):
            locus_number = 0
            title = split[x].split(" ")[0]
            fasta_split = split[x].split("\n")
            fasta_header = ">" + fasta_split[0]
            sequence = "".join(fasta_split[1:])
            sequences_all_regions.append(fasta_header)
            sequences_all_regions.append(sequence)
            forward_annotations = []
            reverse_annotations = []
            
            try:
                forward_annotations, locus_number = get_forward_annotations(title, sequence, locus_number)
            except:
                pass
            try:
                reverse_annotations = get_reverse_annotations(title, sequence, locus_number)
            except:
                pass
                
            annotations_total = forward_annotations + reverse_annotations
            
            annotations = sorted(annotations_total, key=lambda x: int(x.split('\t')[3]))
        
            region_titles += ["##sequence-region " + title + " 1 " + str(len(sequence))]
            annotation_all_regions += annotations
        
        print(header + " annotations complete")

        gff_file = "\n".join(list(region_titles) + list(annotation_all_regions) + list(sequences_all_regions))
        
        total_genes += annotation_all_regions
        outfile_name = os.path.basename(header).split(".fna")[0]
        
        print("writing file")

        outfile_gff = open("PAnnotate/crudely_annotated_100/" + outfile_name + ".gff", "w")
        outfile_gff.write(gff_file)
        outfile_gff.close()
        
        print("file written")
        
    return "yes"
    

if __name__ == '__main__':

    start_time = time.time()
    fnas = glob.glob("PAnnotate/Results1/fasta_files/*.fna")
    chunks = [fnas[i::6] for i in range(6)]
    pool = Pool(processes=6)
    total = pool.map(generate_annotations, chunks)
    end_time = time.time()
    print(str(end_time - start_time) + " seconds" )
    print(str(total))
    

