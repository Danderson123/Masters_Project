#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Using combined_DNA_CDS.fasta but there are random letters interspersed. Removing these for now but should check. 
"""
import numpy as np
import re
from tqdm import tqdm

with open("blpN~~~blpK.aln.fas", "r") as f:
    aligned = f.read().split(">")[1:]

alphabet = {"A":1, "C":2, "T":3, "G":4, "-":"-", "a":1, "c":2, "t":3, "g":4, "N": "*", "K": "*", "R":"*"}

def encode(string):
    string = re.sub(r'[^ACTG]', '', "".join(string.splitlines()[1:]))
    coded_string = ''
    for x in range(len(string)):
        coded_string += str(alphabet[string[x]])
    return coded_string
     

def encode_2(string):
    string = re.sub(r'[^ACTG]', '',string)
    coded_string = ''
    for x in range(len(string)):
        coded_string += str(alphabet[string[x]])
    return coded_string
   
cleaned = []
lengths = []
coded = ''
coded_1 = ''

meh = "".join(aligned[2].splitlines()[1:])
for i in range(len(meh)):
    cleaned = meh[i]
    coded += str(alphabet[meh[i]])

meh = "".join(aligned[4].splitlines()[1:])
for i in range(len(meh)):
    cleaned = meh[i]
    coded_1 += str(alphabet[cleaned])

coded = list(coded)        
coded_1 = list(coded_1)

to_remove = []
for x in range(len(coded)):
    if coded[x] == '-':
        to_remove.append(x)

for x in range(len(coded_1)):
    if coded_1[x] == '-':
        to_remove.append(x)
        
for index in sorted(to_remove, reverse=True):
    del coded_1[index]
    del coded[index]

arr_1 = np.array([int(x) for x in coded])
arr_2 = np.array([int(x) for x in coded_1])

distances = np.subtract(arr_1,arr_2[:,None])

diagonals = np.diagonal(distances)

with open("combined_DNA_CDS.fasta", "r") as f:
    totals = f.read().split(">")[1:]

cleaned_total = []
for x in tqdm(totals):
    cleaned_total.append(encode(x))

query = "atggatacaaaaattatggaacaatttcatgaaatggatataactatgttatctagtattgaaggaggcaagaataattggcaaactaatgtcttagaaggtggtggtgctgcttttggttgtgcagcgggtggagtgaaatatgggagacttctaggaccatggggcgctgcaataggaggaattggaggagcagtggtttgtggatatttagcctataccgctacatcataa".upper()
encoded = np.array([int(x) for x in list(encode_2(query))])

sums = []
for linear in tqdm(cleaned_total):
    arr = np.array([int(x) for x in list(linear)])
    dsds = np.subtract(arr,encoded[:,None])
    sums.append((sum(np.diagonal(dsds))))

found = []
for x in range(len(sums)):
    if sums[x] == 0:
        found.append(x)
    #sums.append(sum(np.diagonal(np.subtract(x,query[:,None]))))
    