#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 21:00:10 2020

@author: danielanderson
"""


import glob 
from tqdm import tqdm

alns = glob.glob("aligned_gene_sequences/*.aln.fas")
fastas = glob.glob("aligned_gene_sequences/*.fasta")

core = []
accessory = []
for aln in tqdm(alns):
    with open(aln, "r") as a:
        file = a.read().split(">")[1:]
    cleaned_seq = []
    for seq in file:
        if not "x" in seq:
            cleaned_seq.append(seq)
    for clean in cleaned_seq:
        sequence = "".join(clean.splitlines()[1:])
        if len(cleaned_seq) > 609:
            core.append(len(sequence))
        elif len(cleaned_seq) <= 609:
            accessory.append(len(sequence))
                
for fasta in tqdm(fastas):
    with open(fasta, "r") as f:
        file = f.read().split(">")[1:]
    for seq in file:
        if not "x" in seq:
            sequence = "".join(clean.splitlines()[1:])
            accessory.append(len(sequence))
            
print("Core: " + str(sum(core)))
print("Accessory: " + str(sum(accessory)))
