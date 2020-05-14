#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:49:12 2020

@author: danielanderson
"""

from Bio import Entrez
from Bio import SeqIO

#Download the genomes in fasta format  

Entrez.email = "danielanderson1@hotmail.com"
Accessions = ['NZ_CP027540.1', 'NC_017592.1', 'NC_012469.1', 'NC_011900.1', 'NZ_AKVY01000001.1', 'NC_003098.1']
search = " ".join(Accessions)

handle = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
genome_ids = handle['IdList']
for genome_id in genome_ids:
    record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
    filename = 'genBankRecord_{}.fasta'.format(genome_id)
    print('Writing:{}'.format(filename))
    with open(filename, 'w') as f:
        f.write(record.read())
    genome = SeqIO.read(('genBankRecord_{}.fasta'.format(genome_id)), "fasta")
    sequence = str(genome.seq)

