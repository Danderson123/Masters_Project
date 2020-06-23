#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
this script was used to run and time prokka to annotate 10 sequences and subsequently count the number of annotations.
"""
import subprocess
import glob
from tqdm import tqdm
import time

start = time.time()

title = glob.glob("Results1/*.fna")

for x in tqdm(range(len(title))):
     command = "prokka --outdir Results1/PROKKA_annotated" + title[x] + " --genus Streptococcus --species pneumoniae --cpus 0 --quiet --norrna --notrna " + title[x]
     subprocess.run(command, shell = True)

end = time.time()
print(str(end - start) + " seconds")

gffs = glob.glob("prokka_gffs/*.gff")

all_annotations = []

for file in gffs:
    with open(file, 'r') as f:
        infile = f.read()
    annotations = infile.splitlines()
    index = 0
    for line in range(len(annotations)):
        if "##sequence-region" in annotations[line] or "##gff-version 3" in annotations[line]:
            index += 1
    cleaned = "\n".join(annotations[index:])
    fasta_split = cleaned.split("##FASTA")[0].splitlines()
    all_annotations += fasta_split
    
print(len(all_annotations))
