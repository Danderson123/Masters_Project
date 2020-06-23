#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:40:01 2020

@author: danielanderson
"""
import subprocess

with open("Spn_artifical_query.fasta", 'rt') as f:
    stored_fasta = str(f.read())

split_fasta = stored_fasta.split(">")[1:]

title = []
sequences = []
for x in range(len(split_fasta)):
    split_fasta[x] = split_fasta[x].upper()  
    split = split_fasta[x].split("\n")
    title.append(split[0])
    sequences.append("".join(split[1:]))
    
for x in range(len(title)):
     outfile = open("split_contigs/" + title[x] + '.fasta','w')
     outfile.write(">" + title[x] + "\n" + sequences[x])
     outfile.close()
     command = "prokka split_contigs/" + title[x] + ".fasta --outdir reannotated/" + title[x] + ' --prefix ' + title[x]
     subprocess.run(command, shell = True)
     
     out_gff = open("test_gffs/" + title[x] + ".gff", 'w')
     
     in_gff = open("reannotated/" + title[x] + "/" + title[x] + ".gff", 'r').read()
     out_gff.write(in_gff)
     out_gff.close()
