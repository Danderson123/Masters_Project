import sys, os
import argparse
from fasta import FastaReader
import numpy as np
import re
import pandas as pd
from tqdm import tqdm

def reverse_complement(dna):
    try:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        reversed = ''.join([complement[base] for base in dna[::-1]])
    except:
        reversed = ""
    return reversed

def convert(gfffile, fastafile, outputfile, count, header):

    all_gene_ids = []
    source = []
    type = []
    strand = []
    phase = []
    attributes = []
    all_sequences = []
    isolates = []
    
    with open(outputfile, 'w') as outfile:
        outfile.write("##gff-version 3\n")
        for h, s in FastaReader(fastafile):
            h = h.split()[0]
            outfile.write(" ".join(["##sequence-region", h, "1",
                                    str(len(s))]) + "\n")

        cleaned_gff = []
        with open(gfffile, 'rU') as infile:
            for line in infile:
                if line[0] != "#":
                    line = line.split('\t')
                    old_ID = re.search('ID=(.*?);', line[8]).group(1)
                    new_ID = str(count) + "_" + old_ID
                    line[8] = line[8].replace(old_ID, new_ID)
                    line = "\t".join(line)
                    cleaned_gff.append(line)
                    outfile.write(line)
    
        outfile.write("##FASTA\n")

        for h, s in FastaReader(fastafile):
            h = h.split()[0]
            for x in cleaned_gff:
                if h in x:
                    tab_splitted = x.split('\t')
                    source.append(tab_splitted[1])
                    type.append(tab_splitted[2])
                    phase.append(tab_splitted[7])
                    attributes.append(tab_splitted[8])
                    all_gene_ids.append(re.search('ID=(.*?);', tab_splitted[8]).group(1))
                    annotation_seq = s[int(tab_splitted[3]) - 1: int(tab_splitted[4]) - 3]
                    strand.append(tab_splitted[6])
                    isolates.append(header)
                    if tab_splitted[6] == "+":
                        all_sequences.append(annotation_seq)
                    elif tab_splitted[6] == "-":
                        all_sequences.append(reverse_complement(annotation_seq))

            outfile.write(">" + h + "\n" + s + "\n")
    
    count += 1
    return count, all_gene_ids, source, type, strand, phase, attributes, all_sequences, isolates


def main():
    all_gene_ids = []
    source = []
    type = []
    strand = []
    phase = []
    attributes = []
    all_sequences = []
    headers = []
    
    import glob
    
    fastas = glob.glob("*.fna")
    gffs = glob.glob("*.gff")
    
    count = 0
    for fasta in tqdm(fastas):
        gff = fasta.split(".fna")[0] + ".gff"
        header = fasta.split(".fna")[0]
        out = "prodigal_joined/" + gff
        count,agi,so,ty,st,ph,att,all,h = convert(gff, fasta, out, count, header)
        all_gene_ids +=agi
        source +=so
        type +=ty
        strand +=st
        phase +=ph
        attributes +=att
        all_sequences +=all
        headers += h
    all_files_annotations = pd.DataFrame({"Isolate": headers, 'ID': all_gene_ids, 'source':source,'type':type, 'strand' : strand, 'phase': phase, 'attributes': attributes, 'sequence': all_sequences})
    all_files_annotations.to_csv("prodigal_joined/" + "all_annotations.csv", index=False)
    return


if __name__ == '__main__':
    main()
