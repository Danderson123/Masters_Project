# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from retrieve_sequences import split_contigs
import glob
import os

headers = glob.glob("sample_gffs/raw_files/*.gff.gz")

clen = []

for x in headers:
    x = os.path.basename(x).split(".gff")[0]
    clen.append(x)
    
split_contigs(clen,
              'sample_gffs/raw_files/',
              'sample_gffs/')
