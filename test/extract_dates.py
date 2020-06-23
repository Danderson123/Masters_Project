#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 14:10:07 2020

@author: danielanderson
"""
import re
import pandas as pd
from Bio import Entrez
import time
from tqdm import tqdm

with open("NC_003098.1.gff", 'r') as f:
    seq = f.read()

seq_split = seq.split("##FASTA")[0]
seq_split = seq_split.split("\n")[1:-1]

gene = []
db = []
Id = []

for x in range(len(seq_split)):
    try:
        gene.append(re.search('gene=(.*?);', seq_split[x]).group(1))
    except:
        try:
            gene.append(re.search('Name=(.*?);', seq_split[x]).group(1))
        except:
            gene.append(re.search('ID=(.*?);', seq_split[x]).group(1))
    try:
        refseq = re.search('RefSeq:(.*?);', seq_split[x]).group(1)
        db.append('RefSeq')
        Id.append(refseq)
    except:
        dbxref = re.search('dbxref=Genbank:(.*?);', seq_split[x]).group(1)
        db.append('Genbank')
        Id.append(dbxref)


information = pd.DataFrame(gene, columns = ["gene"])
information['database'] = db
information['ID'] = Id

create_list = []
update_list = []

start = time.time()
Entrez.email = 'danielanderson1@hotmail.com'
for x in tqdm(range(len(Id))):
    handle = Entrez.read(Entrez.esearch(db="protein", term= Id[x]))
    assembly_id = handle['IdList'][0]
    esummary_handle = Entrez.esummary(db="protein", id=assembly_id, report="full")
    esummary_record = Entrez.read(esummary_handle, validate = False)[0]
    create_list.append(esummary_record['CreateDate'])
    update_list.append(esummary_record['UpdateDate'])       
end = time.time()
print(end - start)

c = open("create.txt", 'w')
u = open("update.txt", 'w')

c.write("\n".join(create_list))
u.write("\n".join(update_list))

