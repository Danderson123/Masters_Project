#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to remove non-gene annotations from NCBI GFF files  
"""

import re 

gff = "sequence.gff3"
o = open(gff, 'r')
D39V = str(o.read()) 

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

d = {"Dbxref=": "dbxref=", "Notes=": "note="}
D39V = replace_all(D39V, d) #Ensure GFF annotation format is consistent with prokka output

gene_list = D39V.splitlines() 
title = gene_list[0]
gene_list = gene_list[2:-1] #remove non-essential content

genes = []

for x in range(len(gene_list)):
    if re.split(r'\t+', gene_list[x])[2] == 'gene':
        genes.append(gene_list[x])

       
cleaned_gffs = "\n".join(str(x) for x in genes)
cleaned_gffs = title + '\n' + cleaned_gffs

#Concatenate downloaded fasta and reformatted GFF
cleaned_gffs = cleaned_gffs + '\n' + '##FASTA' + '\n' + str((open('15902044.fasta','r').read())) 
outfile = open('genes_cleaned.gff','w')          
outfile.write(cleaned_gffs)
outfile.close()

