#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Replace the isolate names in alignment files so they are consistent with the ML phylogeny of the overall population structure
"""

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

import glob

alns = glob.glob("SPARC_AMR_ALIGNMENTS/*.aln.fas")
alns +=  glob.glob("SPARC_AMR_ALIGNMENTS/*.fasta")

name = []
freq = []

for aln in alns:
    with open(aln, "r") as f:
        alnnnt = f.read()
    file = alnnnt.split(">")[1:]
    gene = aln.split("/")[1]
    gene = gene.split(".")[0]
    name.append(gene)
    freq.append((len(file)/616)*100)
    
    for thing in file:
        long_title = thing.splitlines()[0]
        title = long_title.split(";")[0]
        try:
            with open("SPARC_AMR_ALIGNMENTS/data/" + title + ".fna", "r") as da:
                fna = da.read()
            try:
                fna = fna.split("isolate ")[1].split(",")[0]
            except:
                fna = fna.split("strain ")[1].split(",")[0] 
            print(fna)
            d = {long_title: fna}
            alnnnt = replace_all(alnnnt, d)
            out_file = open(aln, "w")
            out_file.write(alnnnt)
            out_file.close()
        except:
            pass
