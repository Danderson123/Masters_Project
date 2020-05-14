#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 11:37:24 2020

@author: danielanderson
"""

import re 

gff = "D39V.gff3"
o = open(gff, 'r')
D39V = str(o.read()) 

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

d = {"Dbxref=": "dbxref=", "Notes=": "note=", "old_locus_tag":"remove"}
D39V = replace_all(D39V, d) #Ensure GFF annotation format is consistent with prokka output

gene_list = D39V.splitlines() 
title = gene_list[0]
gene_list = gene_list[2:-1] #remove non-essential content

clean = [] 
unique_ids = []
IDs_seen = set()

for x in gene_list: #Remove features that are not included in the output of prokka
    pattern = "ID=(.+?);"
    ID = re.search(pattern, x).group(1)
    if ID not in unique_ids:
       unique_ids.append(ID)
       IDs_seen.add(ID)
    else: 
        x = replace_all(x, {"ID=":""})
    line = list(x.split(';'))    
    newlist = [] 
    for y in line:
        if 'ID=' in y or 'Name=' in y or 'dbxref=' in y or 'gene=' in y or 'inference=' in y or 'locus_tag=' in y or 'product=' in y or 'eC_number=' in y or 'note=' in y:
            newlist.append(y)
        filtered = ";".join(str(z) for z in newlist)
    if 'ID=' in filtered:
        clean.append(filtered)
    else:
        pass
       
cleaned_gffs = "\n".join(str(x) for x in clean)
cleaned_gffs = title + '\n' + cleaned_gffs

#Concatenate downloaded fasta and reformatted GFF
cleaned_gffs = cleaned_gffs + '\n' + '##FASTA' + '\n' + str((open('D39V.fasta','r').read())) 
outfile = open('D39V_cleaned.gff','w')          
outfile.write(cleaned_gffs)
outfile.close()

f = str((open('D39V.fasta','r').read())) 
len(f) / 3
list(set(line))
line = f.strip('\n')

len(str((open('D39V.fasta','r').read())))
