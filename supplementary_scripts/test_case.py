#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
this script was used to run and time prokka to annotate 10 sequences and subsequently count the number of annotations.
"""
import subprocess
import glob
from tqdm import tqdm
import time

title = glob.glob("result3_to_annotate/*.fna")

start = time.time()
for x in tqdm(range(len(title))):
     command = "prokka --outdir Results3_PROKKA_annotated_" + title[x] + " --cpus 1  --norrna --notrna " + title[x]
     subprocess.run(command, shell = True)
end = time.time()
print("Prokka: " + str(end - start) + " seconds")

#start = time.time()

#for x in tqdm(range(len(title))):
 #   out = title[x].split(".fna")[0] + ".gff"
  #  command = "prodigal -f gff -i " + title[x] + " -o " + out
   # subprocess.run(command, shell = True)

#end = time.time()
#print("Prodigal: " + str(end - start) + " seconds")

