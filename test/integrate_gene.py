#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:42:09 2020

@author: danielanderson
"""

import glob 
import pandas as pd
import cobs_index as cobs

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

def get_options(): #options for downloading and cleaning
 
    import argparse

    description = 'Use a panaroo-generated graph and mafft alignments to extract updated gffs'
    parser = argparse.ArgumentParser(description=description,
                                  prog='panaroo-update')

    io_opts = parser.add_argument_group('isolate')

    io_opts.add_argument("-i",
                     "--gene_dir",
                     dest="gene_dir",
                     required=True,
                     help='directory containing genes to integrate',
                     type=str)
    
    io_opts.add_argument("-d",
                        "--panaroo_dir",
                        dest="pan_dir",
                        required=True,
                        help='directory of panaroo output',
                        type=str)
                    
    args = parser.parse_args()

    return (args)

def main():
    args = get_options()
    
    path = glob.glob(args.gene_file)
    
    fastas = []
    for file in path:
        with open(file, "r") as f:
            fastas.append(f.read())
    
    proteins = []
    for fasta in fastas:
        proteins.append(translate(fasta))

gene = pd.read_csv('gene_data.csv')

index = []
for x in range(len(gene["gene_name"])):
    if gene["gene_name"][x] == "rsmI":
        index.append(gene["prot_sequence"][x])
    
indels = index[:20]

kmers = []
for y in indels:
    sub_kmers = []
    for z in range(len(y)):
        sub_kmers.append(y[z:z+31])
    kmers.append(sub_kmers)


