#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script was use to retrieve the 616 SPARC project genomic sequences and functional annotations. 
Using retrieve sequences tends to crash after ~400 isolates so this makes it easier to resume downloading if crashes occur.
The same exclusion criteria as in "retrieve_sequences.py" are applied. 
"""
from Bio import Entrez
import os
import urllib.request
from tqdm import tqdm

with open("SPARC_2.txt", "r") as s:
    sample_accessions = s.read()

sample_list = sample_accessions.split("\n")

cleaned = []
for x in range(len(sample_list)):
    if not sample_list[x] == "":
        cleaned.append(sample_list[x])
sample_input = ",".join(cleaned)

Entrez.email = ""

success = []
for item in tqdm(range(len(cleaned))):
    
    handle = Entrez.read(Entrez.esearch(db="assembly", term=cleaned[item], retmax = 1))
    
    assembly_ids = handle['IdList']
    
    #try:
    for assembly_id in assembly_ids:
        esummary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
    
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        
        if url == '':
            url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
            
        label = os.path.basename(url)
        fasta_link = os.path.join(url,label+'_genomic.fna.gz')
        gff_link = os.path.join(url,label+'_genomic.gff.gz')
        
        gff_file = '{}.gff.gz'.format(label)
        fasta_file = '{}.fna.gz'.format(label)
        
        urllib.request.urlretrieve(gff_link, gff_file)
        urllib.request.urlretrieve(fasta_link, fasta_file)
        success.append(cleaned[item])
    #except:
       # continue

success_out = open("successes.txt", "w")
success_out.write("\n".join(success))
success_out.close()

from retrieve_sequences import split_contigs
import glob
import os

headers = glob.glob("*.gff.gz")

clen = []

for x in headers:
    x = os.path.basename(x).split(".gff")[0]
    clen.append(x)
    
split_contigs(clen,
              '',
              '')
