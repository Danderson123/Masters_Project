#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 19:15:05 2020

@author: danielanderson
"""


import subprocess
import glob
from tqdm import tqdm

#title = glob.glob("*.fna")

#headers = []
#for x in tqdm(range(len(title))):
   # out = title[x].split(".fna")[0]
   # headers.append(out)
   # input_fna = out + ".fna"
   # out = out + ".gff"
    #command = "prodigal -f gff -i " + input_fna + " -o " + out + " -q"
   # subprocess.run(command, shell = True)

#subprocess.run("gzip *.gff", shell = True)
 
from retrieve_sequences import *

head = glob.glob("*.gff.gz")
headers = []
for h in head:
    h = h.split(".gff.gz")[0]
    headers.append(h)

#split_contigs(headers,
             # '',
              #'prodigal_joined/')

output_dir = 'prodigal_joined/'
input_dir = ""



import time
import os
from Bio import Entrez
import urllib.request
import re
import sys
import shutil
import tempfile
from tqdm import tqdm
import gzip
import pandas as pd
               
all_gene_ids = []
source = []
type = []
phase = []
attributes = []

all_raw = []
non_cds = []
non_divis = []
duplicates = []
not_start = []
stop_in =[]

for header in tqdm(headers):
    d = {"[": "", "]": "", " ": "", "'": ""}
    header = replace_all(header, d)
    fasta = input_dir + header + '.fna.gz'
    gff = input_dir + header + '.gff.gz'
    
    try:
        fasta = input_dir + header + '.fna.gz'
        with gzip.open(fasta, 'rt') as f:
            stored_fasta = str(f.read())
    except:
        fasta = input_dir + header + '.fna'
        with open(fasta, 'rt') as f:
            stored_fasta = str(f.read())
            
    with gzip.open(gff, 'rt') as g:
        stored_gff = str(g.read())

    stored_fasta = stored_fasta.split('>')
    stored_fasta = stored_fasta[1:]
    
    stored_gff = stored_gff.split("# Sequence Data:")
    stored_gff = stored_gff[1:]
    
    all_region_names = []
    all_region_annotations = []
    all_region_sequences = []
    
    for gff_region in range(len(stored_gff)):

        d = {"Dbxref=": "dbxref=", "Notes=": "note="}
        stored_gff[gff_region] = replace_all(stored_gff[gff_region], d) #Ensure GFF annotation format is consistent with prokka output
        
        cleaned_gene_list = []
        gene_list = stored_gff[gff_region].splitlines()
        title = gene_list[0].split("\n")[0]
        for line in range(len(gene_list)):
          if not '#' in gene_list[line]:
              cleaned_gene_list.append(gene_list[line])

        gene_list = gene_list[2:]
        
        print(cleaned_gene_list)
        genes = []
        for y in range(len(cleaned_gene_list)):
            all_raw.append(cleaned_gene_list[y])
            try:
                split = re.split(r'\t+', cleaned_gene_list[y])
            except:
                print(split)
            end = int(split[4])
            start = int(split[3])
            if split[2] == 'CDS' and (((end - start) + 1) % 3) == 0: #process_pokka input only looks for CDSs and returns duplicate error when the exon is split.
                genes.append(cleaned_gene_list[y])
            elif not split[2] == 'CDS':
                non_cds.append(gene_list[y])
            elif split[2] == 'CDS' and not (((end - start) + 1) % 3) == 0:
                non_divis.append(cleaned_gene_list[y])
        
        Ids= []
        for gene in range(len(genes)):
            Id = re.search('ID=(.*?);', genes[gene]).group(1)
            Ids.append(Id)

        output_cds = []
        start_index = []
        end_index = []
        sense = []
        for Id in range(len(Ids)):
            if Ids.count(Ids[Id]) == 1:
              output_cds.append(genes[Id])
              no_duplicates_split = re.split(r'\t+', genes[Id])
              start_index.append(int(no_duplicates_split[3]) - 1)
              end_index.append(int(no_duplicates_split[4]) - 1)
              sense.append(no_duplicates_split[6])
            else:
                duplicates.append(Ids[Id])
        
        nucleotides = ''.join(stored_fasta[gff_region].split('\n')[1:])

        to_remove = []
        for sign in range(len(sense)):
            if sense[sign] == '+':
                start_codon = nucleotides[start_index[sign]: start_index[sign] + 3]
                protein = translate(nucleotides[start_index[sign]: end_index[sign] - 2])
                if not (start_codon == "ATG" or start_codon == "TTG" or start_codon == "GTG") or "*" in protein or protein == "":
                    to_remove.append(sign)
                    if not (start_codon == "ATG" or start_codon == "TTG" or start_codon == "GTG"):
                        not_start.append(start_codon)
                    if "*" in protein:
                        stop_in.append(protein)
            else:
                reverse_start_codon = nucleotides[end_index[sign] - 2: end_index[sign] + 1]
                protein_reverse = translate(reverse_complement(nucleotides[start_index[sign] + 3: end_index[sign] + 1]))
                if not (reverse_start_codon == "CAT" or reverse_start_codon == "CAA" or reverse_start_codon == "CAC") or "*" in protein_reverse or protein_reverse == "":
                    to_remove.append(sign)
                    if not (reverse_start_codon == "CAT" or reverse_start_codon == "CAA" or reverse_start_codon == "CAC"):
                        not_start.append(start_codon)
                    if "*" in protein_reverse:
                        stop_in.append(protein_reverse)
        
        for index in sorted(to_remove, reverse=True):
            del output_cds[index]
                
        if len(output_cds) == 0:
          continue

        all_region_names.append("##sequence-region " + title)
                    
        cleaned_gffs = "\n".join(str(z) for z in output_cds)
        
        for annotation_line in range(len(output_cds)):
            tab_splitted = output_cds[annotation_line].split('\t')
            source.append(tab_splitted[1])
            type.append(tab_splitted[2])
            phase.append(tab_splitted[7])
            attributes.append(tab_splitted[8])
            all_gene_ids.append(re.search('ID=(.*?);', tab_splitted[8]).group(1))
            
        all_region_annotations.append(cleaned_gffs)
        #Concatenate downloaded fasta and reformatted GFF
        fasta_file = ">" + stored_fasta[gff_region]
        all_region_sequences.append(fasta_file)
        
    annotated_file = "\n".join(all_region_names + all_region_annotations) + "\n##FASTA\n" + "".join(all_region_sequences)
    filename_cleaned = output_dir + header + '.gff'
    outfile = open(filename_cleaned,'w')
    outfile.write(annotated_file)
    outfile.close()

all_files_annotations = pd.DataFrame({'ID': all_gene_ids, 'source':source,'type':type, 'phase': phase, 'attributes': attributes})

all_files_annotations.to_csv(output_dir + "all_annotations.csv", index=False)

print("raw annotations:" + str(len(all_raw)))
print("Non-CDS:" + str(len(non_cds)))
print("Non-divisible:" + str(len(non_divis)))
print("duplicates:" + str(len(duplicates)))
print("not_start:" + str(len(not_start)))
print("stop_in:" + str(len(stop_in)))
print("all_valid:" + str(len(all_gene_ids)))
