#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script counts the number of annotation lines in a prokka-outputted gff  
"""


import glob 

files = glob.glob("*.gff")
total_length = []

for file in files:
    with open(file, 'r') as f:
        gff = f.read()
    
    gff = (((gff.split("##FASTA")[0]).split("##sequence-region")[-1]).splitlines()[1:])
    length = len(gff)
    total_length.append(length)
    
sum(total_length)
