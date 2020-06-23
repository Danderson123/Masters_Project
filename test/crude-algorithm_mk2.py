#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests for the crude_annotate script
"""

import re
import numpy as np

with open("string.txt") as s:
    string = s.read()

sequence = "".join((string.split(">")[5].splitlines())[1:])

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])
    
    
def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary  
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] #Defining a codon as 3 bases 
            protein+= table[codon] #Translate the codon into an amino acid based on the dictionary and append this to the protein sequence 
    return protein     

def get_forward_annotations(sequence):
        
    start1 = r"ATG"
    start2 = r"TTG"
    start3 = r"GTG"
    
    stop1 = r"TAG"
    stop2 = r"TAA"
    stop3 = r"TGA"
    
    forwardstartList = [start1, start2, start3]
    start_indexes_forward = []

    for start_regex in forwardstartList:
        if re.search(start_regex, sequence):
            stop_codon_list = re.finditer(start_regex,sequence)
            for y in stop_codon_list:
                start_indexes_forward.append(y.start())
    #start_codons_forward = re.finditer(r"ATG", sequence)
    
    stop_indexes_forward = []

    #for x in start_codons_forward:
        #start_indexes_forward.append(x.start() + 1)
    
    forwardstopList = [stop1, stop2, stop3]#, stop4, stop5, stop6, stop7, stop8, stop9]
    
    for stop_regex in forwardstopList:
        if re.search(stop_regex, sequence):
            stop_codon_list = re.finditer(stop_regex,sequence)
            for y in stop_codon_list:
                stop_indexes_forward.append(y.end() -1)
            
    start_indexes_forward = sorted(start_indexes_forward)
    stop_indexes_forward = sorted(stop_indexes_forward)
    
    start_array = np.array(start_indexes_forward)
    stop_array = np.array(stop_indexes_forward)
    
    distances = np.subtract(stop_array,start_array[:,None])
    #medial = np.where(((distances < 0) | (distances == 0)), 2, distances)
    distances1 = distances +1
    middle = np.where(distances1 <= 1, 2, distances1)
    
    final = np.where((middle > 2) & (middle % 3 == 0), 1, middle)
    
    import pandas as pd
    minimum = np.min(final, axis=1)
    itemindex = pd.DataFrame(np.argmin(final,axis=1), columns = ["index"])
    
    to_remove = []

    for x in range(len(minimum)):
        if not minimum[x] == 1:
            to_remove.append(x)
    
    itemindex = itemindex.drop(itemindex.index[to_remove])
        
    annotations = []
    sbtr = []
    sequence_list = []
    
    locus_number = 0
            
    for x in range(len(itemindex['index'])):
        start_codon = start_indexes_forward[x]
        stop_codon = stop_indexes_forward[itemindex['index'][x]]
        if not start_codon > stop_codon:
            annotation = title
            annotation += "\tCRUDE\tCDS\t"
            annotation += str(start_codon + 1)
            annotation += "\t"
            annotation += str(stop_codon + 1)
            annotation += "\t.\t"
            annotation += "+"
            annotation += "\t0\tID="
            annotation += title.split(".")[0]
            annotation += "_"
            annotation += str(locus_number)
            annotation += ";product=putative_protein_region"
            annotations.append(annotation)
            
            locus_number += 1
            
            substring = sequence[start_codon:stop_codon + 1]
            sbtr.append(substring)
            substring = translate(substring)
            sequence_list.append(substring)
            
    annotations = sorted(annotations, key=lambda x: int(x.split('\t')[3]))
    
    return annotations, locus_number


def get_reverse_annotations(sequence, locus_number):

    sequence = reverse_complement(sequence)
    
    start1 = r"ATG"
    start2 = r"TTG"
    start3 = r"GTG"
    
    stop1 = r"TAG"
    stop2 = r"TAA"
    stop3 = r"TGA"
    
    forwardstartList = [start1, start2, start3]
    start_indexes_forward = []

    for start_regex in forwardstartList:
        if re.search(start_regex, sequence):
            stop_codon_list = re.finditer(start_regex,sequence)
            for y in stop_codon_list:
                start_indexes_forward.append(y.start())
    #start_codons_forward = re.finditer(r"ATG", sequence)
    
    stop_indexes_forward = []

    #for x in start_codons_forward:
        #start_indexes_forward.append(x.start() + 1)
    
    forwardstopList = [stop1, stop2, stop3]#, stop4, stop5, stop6, stop7, stop8, stop9]
    
    for stop_regex in forwardstopList:
        if re.search(stop_regex, sequence):
            stop_codon_list = re.finditer(stop_regex,sequence)
            for y in stop_codon_list:
                stop_indexes_forward.append(y.end() -1)
            
    start_indexes_forward = sorted(start_indexes_forward)
    stop_indexes_forward = sorted(stop_indexes_forward)
    
    start_array = np.array(start_indexes_forward)
    stop_array = np.array(stop_indexes_forward)
    
    distances = np.subtract(stop_array,start_array[:,None])
    #medial = np.where(((distances < 0) | (distances == 0)), 2, distances)
    distances1 = distances +1
    middle = np.where(distances1 <= 1, 2, distances1)
    
    final = np.where((middle > 2) & (middle % 3 == 0), 1, middle)
    
    import pandas as pd
    minimum = np.min(final, axis=1)
    itemindex = pd.DataFrame(np.argmin(final,axis=1), columns = ["index"])
    
    to_remove = []

    for x in range(len(minimum)):
        if not minimum[x] == 1:
            to_remove.append(x)
    
    itemindex = itemindex.drop(itemindex.index[to_remove])
        
    annotations = []
    sbtr = []
    sequence_list = []
    
    seq_length = len(sequence)
    for x in range(len(itemindex['index'])):
        start_codon = seq_length - int(start_indexes_forward[x])
        stop_codon = seq_length - int(stop_indexes_forward[itemindex['index'][x]])       
        if not start_codon < stop_codon:
            annotation = title
            annotation += "\tCRUDE\tCDS\t"
            annotation += str(stop_codon)
            annotation += "\t"
            annotation += str(start_codon)
            annotation += "\t.\t"
            annotation += "-"
            annotation += "\t0\tID="
            annotation += title.split(".")[0]
            annotation += "_"
            annotation += str(locus_number)
            annotation += ";product=putative_protein_region"
            annotations.append(annotation)
            
            locus_number += 1
            
            substring = reverse_complement(sequence)
            substring = substring[stop_codon-1:start_codon]
            sbtr.append(substring)
            aa = translate(reverse_complement(substring))
            sequence_list.append(aa)
            
    annotations = sorted(annotations, key=lambda x: int(x.split('\t')[3]))
    
    return annotations


































def forward_stop_starts(title, sequence):
    start_codons_forward = re.finditer(r"ATG", sequence)
    
    start_indexes_forward = []
    stop_indexes_forward = []

    for x in start_codons_forward:
        start_indexes_forward.append(x.start() + 1)
    
    stop1 = r"TAG"
    stop2 = r"TAA"
    stop3 = r"TGA"
    
    forwardstopList = [stop1, stop2, stop3]
    
    
    for stop_regex in forwardstopList:
        if re.search(stop_regex, sequence):
            stop_codon_list = re.finditer(stop_regex,sequence)
            for y in stop_codon_list:
                stop_indexes_forward.append(y.end())
                
    start_array = np.array(start_indexes_forward)
    stop_array = np.array(stop_indexes_forward)
    
    distances = np.subtract(stop_array,start_array[:,None])
    medial = np.where(((distances < 0) | (distances == 0)), 2, distances)
    final = np.where((medial > 2) & (medial % 3 == 0), True, medial)
    itemindex = np.argmax(final==1, axis=1)
    
    to_remove = []
    for x in range(len(itemindex)):
        if itemindex[x] < itemindex[x-1] and itemindex[x] == 0:
            to_remove.append(x)
    
    itemindex = list(itemindex)
    for index in sorted(to_remove, reverse=True):
        del itemindex[index]
        
    annotations = []
    
    locus_number = 0
    for x in range(len(itemindex)):
        annotation = title 
        annotation += "\tCRUDE\tCDS\t" 
        annotation += str(start_indexes_forward[x]) 
        annotation += "\t" 
        annotation += str(stop_indexes_forward[itemindex[x]])
        annotation += "\t.\t" 
        annotation += "+" 
        annotation += "\t0\tID=" 
        annotation += title.split(".")[0] 
        annotation += "_" 
        annotation += str(locus_number) 
        annotation += ";product=putative_protein_region"
        annotations.append(annotation)
        
        locus_number += 1
        
    return locus_number, annotations
    


def reverse_stop_starts(title, sequence, locus_number):
    
    start1 = r"CTA"
    start2 = r"TTA"
    start3 = r"TCA"
    
    reversestoplist = [start1, start2, start3]
    
    stop_indexes_reverse = []

    for stop_regex in reversestoplist:
        if re.search(stop_regex, sequence):
            stop_codon_list = re.finditer(stop_regex,sequence)
            for y in stop_codon_list:
                stop_indexes_reverse.append(y.start() + 1)
                
    start_codons_reverse= re.finditer(r"CAT", sequence)
    
    start_indexes_reverse = []
    
    for x in start_codons_reverse:
        start_indexes_reverse.append(x.end())
    
   
    start_array = np.array(stop_indexes_reverse)
    stop_array = np.array(start_indexes_reverse)
    
    distances = np.subtract(stop_array,start_array[:,None])
    medial = np.where(((distances < 0) | (distances == 0)), 2, distances)
    final = np.where((medial > 2) & (medial % 3 == 0), True, medial)
    itemindex = np.argmax(final==1, axis=1)
    
    to_remove = []
    for x in range(len(itemindex)):
        if itemindex[x] < itemindex[x-1] and itemindex[x] == 0:
            to_remove.append(x)
    
    itemindex = list(itemindex)
    for index in sorted(to_remove, reverse=True):
        del itemindex[index]
        
    annotations = []

    for x in range(len(itemindex)):
        annotation = title 
        annotation += "\tCRUDE\tCDS\t" 
        annotation += str(stop_indexes_reverse[x]) 
        annotation += "\t" 
        annotation += str(start_indexes_reverse[itemindex[x]])
        annotation += "\t.\t" 
        annotation += "-" 
        annotation += "\t0\tID=" 
        annotation += title.split(".")[0] 
        annotation += "_" 
        annotation += str(locus_number) 
        annotation += ";product=putative_protein_region"
        annotations.append(annotation)
        
        locus_number += 1
    return annotations


for x in tqdm(range(len(split))):
    locus_number = 0
    title = split[x].split(" ")[0]
    fasta_split = split[x].split("\n")
    fasta_header = ">" + fasta_split[0]
    sequence = "".join(fasta_split[1:])
    sequences_all_regions.append(fasta_header)
    sequences_all_regions.append(sequence)
    #region_titles.append("##sequence-region " + title + " 1 " + str(len(sequence)))
    annotations = []
    
    title = "booobaler"
    locus_number, forward_annotations = forward_stop_starts(title, string)
    reverse_annotations = reverse_stop_starts(title, string, locus_number)
    
    annotations = forward_annotations + reverse_annotations
    
    annotations = sorted(annotations, key=lambda x: int(x.split('\t')[3]))
    
    
    
    for regex in regexList:
        if re.findall(regex, sequence):
            start_list = re.finditer(regex,sequence)
            start_codon = regex.split("(")[0]
            
            if start_codon == "ATG":
                gene_strand = "+"
                for stop_regex in forwardstopList:
                    if re.search(stop_regex, sequence):
                        stop_codon_list = re.finditer(stop_regex,sequence)
                        for result in start_list:
                            stop_result_end_list = []
                            for stop_result in stop_codon_list:
                                if (stop_result.end() - result.start()) % 3 == 0 and stop_result.end() - result.start() >=1:
                                    stop_result_end_list.append(stop_result.end())
                            if not stop_result_end_list == []:
                                stop_end = min(stop_result_end_list)
                                annotation = title + "\tCRUDE\tCDS\t" + str(result.start() - 2) + "\t" + str(stop_result.end() + 3) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                                annotations.append(annotation)
                                locus_number += 1
                            else:
                                pass
                
                
            elif start_codon == "CTA" or start_codon == "TTA" or start_codon == "TCA":
                gene_strand = "-"
                for stop_regex in reversestoplist:
                    if re.search(stop_regex, sequence):
                        stop_codon_list = re.finditer(stop_regex,sequence)
                        for result in start_list:
                            stop_result_end_list = []
                            for stop_result in stop_codon_list:
                                if (stop_result.end() - result.start()) % 3 == 0 and stop_result.end() - result.start() >=1:
                                    stop_result_end_list.append(stop_result.end())
                            if not stop_result_end_list == []:
                                stop_end = min(stop_result_end_list)
                                annotation = title + "\tCRUDE\tCDS\t" + str(result.start() - 2) + "\t" + str(stop_result.end() + 3) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                                annotations.append(annotation)
                                locus_number += 1
                            else:
                                pass
            else:
                print("oh no")
        else:
            print("no regex")
                                      
    annotation_all_regions += ["##sequence-region " + title + " 1 " + str(len(sequence))]
    annotation_all_regions = sorted(annotations, key=lambda x: x.split('\t')[3])
            
        
        
    #gff_file = "\n".join(region_titles + annotation_all_regions + sequences_all_regions)
    gff_file = "\n".join(list(annotation_all_regions))
    fasta_out = "\n".join(list(sequences_all_regions))
    
    total_genes += annotation_all_regions
    #outfile_name = os.path.basename(header).split(".fna")[0] + ".gff"
    outfile_name = os.path.basename(header).split(".fna")[0]
    
    outfile_fasta = open("crudely_annotated/" + outfile_name + ".fna", "w")
    outfile_fasta.write(fasta_out)
    outfile_fasta.close()
    
    outfile_gff = open("crudely_annotated/" + outfile_name + ".gff", "w")
    outfile_gff.write(gff_file)
    outfile_gff.close()
