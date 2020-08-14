#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This was used to run prodigal on the 616 S. pneumoniae genomes to generate non-functional anntoations
"""

import subprocess
import glob
from tqdm import tqdm

title = glob.glob("*.fna")

headers = []
for x in tqdm(range(len(title))):
    out = title[x].split(".fna")[0]
    headers.append(out)
    input_fna = out + ".fna"
    out = out + ".gff"
    command = "prodigal -f gff -i " + input_fna + " -o " + out + " -q"
    subprocess.run(command, shell = True)

subprocess.run("gzip *.gff", shell = True)
