#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Search an index for all non-reference sequences in an "all_annotations.csv" file at different query thresholds. 
"""

import cobs_index as cobs
from tqdm import tqdm
import pandas as pd

def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary
    """"Translate nucleotide sequences"""
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    protein =""
    try:
        if len(seq)%3 == 0:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3] #Defining a codon as 3 bases
                protein+= table[codon] #Translate the codon into an amino acid based on the dictionary and append this to the protein sequence
    except:
        protein = ""
    return protein

def identify(strand, sequence, index, thresh, isolate):
    #if cora == "accessory":
    if len(sequence) >= 30 and not "_ASM" in isolate:
        try:
            if strand == "+":
                sequence = translate(sequence)
                if not sequence == "":
                    result = index.search(sequence, threshold = thresh)
                else:
                    result = "CHARACTERS"
                if not len(result) == 0:
                    result = (str(result[0]).split("'")[1]).split(";")[0]
                else:
                    result = "[]"
            elif strand == "-":
                sequence = translate(sequence)
                if not sequence == "":
                    result = index.search(sequence, threshold = thresh)
                else:
                    result = "CHARACTERS"
                if not len(result) == 0:
                    result = (str(result[0]).split("'")[1]).split(";")[0]
                else:
                    result = "[]"
            #else:
               # result = ""
            if result == "[]":
                result = "NULL"
        except:
            result = "ERROR"
    else:
        result = "LENGTH/ISOLATE"
    return result
    
index = cobs.Search("INDEX/10_index_index.cobs_compact")

all_annotations = pd.read_csv("reference_sparc_merged/all_annotations.csv")

tqdm.pandas()

#all_annotations["COBS 75%"] = all_annotations.progress_apply(lambda row: identify(row["strand"], row["sequence"], index, 0.75, row["Isolate"]), axis = 1)
#all_annotations["COBS 80%"] = all_annotations.progress_apply(lambda row: identify(row["strand"], row["sequence"], index, 0.8, row["Isolate"]), axis = 1)
all_annotations["COBS 85%"] = all_annotations.progress_apply(lambda row: identify(row["strand"], row["sequence"], index, 0.85, row["Isolate"]), axis = 1)
#all_annotations["COBS 90%"] = all_annotations.progress_apply(lambda row: identify(row["strand"], row["sequence"], index, 0.9, row["Isolate"]), axis = 1)
#all_annotations["COBS 95%"] = all_annotations.progress_apply(lambda row: identify(row["strand"], row["sequence"], index, 0.95, row["Isolate"]), axis = 1)

all_annotations.to_csv("reference_sparc_merged/85_SEARCHED_all_annotations.csv", index = False)

def count(name, list):
    length = 616
    list = "\n".join(list)
    sum = list.count(name)
    freq = sum / 616 * 100
    return freq
    
#all_annotations["freq 75%"] = all_annotations(lambda row: count(row["COBS 75%"], list(all_annotations["COBS 75%"])), axis = 1)
#all_annotations["freq 80%"] = all_annotations(lambda row: count(row["COBS 80%"], list(all_annotations["COBS 80%"])), axis = 1)
all_annotations["freq 85%"] = all_annotations(lambda row: count(row["COBS 85%"], list(all_annotations["COBS 85%"])), axis = 1)
#all_annotations["freq 90%"] = all_annotations(lambda row: count(row["COBS 90%"], list(all_annotations["COBS 90%"])), axis = 1)
#all_annotations["freq 95%"] = all_annotations(lambda row: count(row["COBS 95%"], list(all_annotations["COBS 95%"])), axis = 1)

all_annotations.to_csv("reference_sparc_merged/85_SEARCHED_all_annotations.csv", index = False)
