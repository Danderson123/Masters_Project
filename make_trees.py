#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses subprocess to generate fast trees for each gene in an alignment file.
"""


import subprocess
import glob

paths_aln = "*.aln.fas"
aln_files = glob.glob(paths_aln)

for files in range(len(aln_files)):
    #filename = aln_files[files].split('/')[1]
    #filename = filename.split('.')[0]
    filename = aln_files[files].split('.')[0]
    subprocess_command = ("FastTree -gtr -nt < {}.aln.fas >").format(filename, )
    subprocess.run(subprocess_command, shell = True)
    
