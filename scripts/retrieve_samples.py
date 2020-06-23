#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 15:27:42 2020

@author: danielanderson
"""
from Bio import Entrez
import os
import urllib.request
from tqdm import tqdm

with open("data/biosample_result.txt", "r") as s:
    sample_accessions = s.read()

sample_list = sample_accessions.split("\n")

cleaned = []
for x in range(len(sample_list)):
    if not sample_list[x] == "":
        cleaned.append(sample_list[x])
sample_input = ",".join(cleaned)

Entrez.email = "danielanderson1@hotmail.com"

for item in tqdm(range(len(cleaned))):
    
    handle = Entrez.read(Entrez.esearch(db="assembly", term=cleaned[item], retmax = 1))
    
    assembly_ids = handle['IdList']
    
    try:                        
        for assembly_id in assembly_ids:#tqdm(assembly_ids):
            esummary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
            esummary_record = Entrez.read(esummary_handle, validate = False)
        
            url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
            
            if url == '':
                url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                
            label = os.path.basename(url)
            #get the fasta link - change this to get other formats
            fasta_link = os.path.join(url,label+'_genomic.fna.gz')
            fasta_file = '{}.fna.gz'.format(label)
            
            urllib.request.urlretrieve(fasta_link, fasta_file)
    except:
        print(cleaned[item])
        continue
