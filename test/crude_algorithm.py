#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Initial testing script for crude_annotate.py

"""


string = "TGATGCTGAGACTGACGCACTTGTTGATGCTGAAGCTGACGCACTGGTTGATGCGGATTCAGAAGCTGATGTGCTTGCTGAGGCTGACGCACTCGTTGATGCTGAGACTGACGCACTTGTTGACGCTGAGGCTGAGGCACTTGTCGACGCTGAAGCTGACGCACTCGTTGATGCTGAAGCCGAGGCACTGGTTGACGCTGAAGCTGAAGCACTGGTTGACGCTGAGGCTGACACACTCGTTGATGCTGAGGCTGACGCACTGGTACTTGCTGAAGCCGAGGCGCTAGTGCTTGCTG"

genes = []
start_indexes = []
end_indexes = []

for x in range(len(string)):
    gene = ''
    index = 0
    if string[x:x+3] == "ATG":
        start = x #start is the index of ATG
        gene = string[x:x+3] #gene is ATG 
        cut = string[x:] #ATG to end of string 
        
        end = 0
        for y in range(len(cut)): # for character in the cut
            if not cut[y:y+3] == "CTA": #if the character:character + 3 not CTA
                gene += cut[y:y+3] #append the codon to the gene
                end = y+2 #new y becomes y + 3
                print(end)
            else:
                genes.append(gene + "CTA")
                start_indexes.append(start+1)
                end_of_cut = start + y
                end_indexes.append(y + 7)
                continue
            


genes = []
start_indexes = []
end_indexes = []
cuts = []
strands = []

for x in range(len(string)):
    index = x
    if string[index:index+3] == "ATG":
        gene = ''
        end_codon = ["TAG", "TAA", "TGA"]
        gene_strand = "+"
        start = index #start is the index of ATG
        start_indexes.append(start + 1)
        cut = string[index:] #ATG to end of string 
        cuts.append(cut)
        for y in range(len(cut)): # for character in the cut
            if (y + 3) % 3 ==0:
                gene+=cut[y:y+3]
                if any(cut[y:y+3] == motif for motif in end_codon):
                    genes.append(gene)
                    strands.append(gene_strand)
                    end_indexes.append(start + (y+3))
                    
                    
              
    elif any(string[index:index+3] == rev_stop for rev_stop in ["CTA", "TTA", "TCA"]):
        gene = ''
        gene_strand = "-"
        start = index #start is the index of ATG
        start_indexes.append(start + 1)
        cut = string[index:] #ATG to end of string 
        cuts.append(cut)
        for y in range(len(cut)): # for character in the cut
            if (y + 3) % 3 ==0:
                gene+=cut[y:y+3]
                if cut[y:y+3] == "CAT":
                    genes.append(gene)
                    strands.append(gene_strand)
                    end_indexes.append(start + (y+3))
                    
            

for x in tqdm(range(len(split))):
        locus_number = 0
        title = split[x].split(" ")[0] 
        fasta_split = split[x].split("\n")
        fasta_header = ">" + fasta_split[0] 
        sequence = "".join(fasta_split[1:])
        sequences_all_regions.append(fasta_header)
        sequences_all_regions.append(sequence)
        #region_titles.append("##sequence-region " + title + " 1 " + str(len(sequence)))
        annotations = set()
        for regex in regexList:
                if re.search(regex, sequence):
                    some_list= re.finditer(regex,sequence)
                    for result in some_list:
                        if (result.end() - result.start()) % 3 == 0:
                            #codon = regex.split("(")[0]
                            start_codon = regex[regex.find("=")+1:regex.find(")")]
                            print(start_codon)
                            if start_codon == "ATG":
                                gene_strand = "+"
                            elif start_codon == "CTA" or start_codon == "TTA" or start_codon == "TCA":
                                gene_strand = "-"
                            else:
                                gene_strand = ""
                            annotation = title + "\tCRUDE\tCDS\t" + str(result.start() - 2) + "\t" + str(result.end() + 3) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                            annotations.add(annotation)
                            locus_number += 1
        annotation_all_regions += ["##sequence-region " + title + " 1 " + str(len(sequence))]
        annotation_all_regions += sorted(list(annotations), key=lambda x: x.split('\t')[3])
        
        
        
for x in tqdm(range(len(split))):
    locus_number = 0
    title = split[x].split(" ")[0] 
    fasta_split = split[x].split("\n")
    fasta_header = ">" + fasta_split[0] 
    sequence = "".join(fasta_split[1:])
    sequences_all_regions.append(fasta_header)
    sequences_all_regions.append(sequence)
    #region_titles.append("##sequence-region " + title + " 1 " + str(len(sequence)))
    annotations = set()
    
    for base in range(len(sequence)):
        index = base
        if string[index:index+3] == "ATG":
            gene = ''
            end_codon = ["TAG", "TAA", "TGA"]
            gene_strand = "+"
            start = index #start is the index of ATG
            cut = string[index:] #ATG to end of string 
            for y in range(len(cut)): # for character in the cut
                if (y + 3) % 3 ==0:
                    gene+=cut[y:y+3]
                    if any(cut[y:y+3] == motif for motif in end_codon):
                        annotation = title + "\tCRUDE\tCDS\t" + str(start + 1) + "\t" + str(start + (y+3)) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                        annotations.add(annotation)
                        locus_number += 1
                        
        elif any(string[index:index+3] == rev_stop for rev_stop in ["CTA", "TTA", "TCA"]):
            gene = ''
            gene_strand = "-"
            start = index #start is the index of ATG
            cut = string[index:] #ATG to end of string 
            for y in range(len(cut)): # for character in the cut
                if (y + 3) % 3 ==0:
                    gene+=cut[y:y+3]
                    if cut[y:y+3] == "CAT":
                        annotation = title + "\tCRUDE\tCDS\t" + str(start + 1) + "\t" + str(start + (y+3)) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                        annotations.add(annotation)
                        locus_number += 1
                    
        annotation_all_regions += ["##sequence-region " + title + " 1 " + str(len(sequence))]
        annotation_all_regions += sorted(list(annotations), key=lambda x: x.split('\t')[3])
        
        
        
        
region_titles = []
annotation_all_regions = []
sequences_all_regions = ["##FASTA"]

import re
from tqdm import tqdm 

regex1 = r"ATG(.*)(?=TAG)"
regex2 = r"ATG)(.*)(?=TAA)"
regex3 = r"ATG)(.*)(?=TGA)"
regex4 = r"CTA)(.*)(?=CAT)"
regex5 = r"TTA)(.*)(?=CAT)"
regex6 = r"TCA)(.*)(?=CAT)"

stop1 = r"(?<=ATG)(.*)TAG"
stop2 = r"(?<=ATG)(.*)TAA"
stop3 = r"(?<=ATG)(.*)TGA"
stop4 = r"(?<=CTA)(.*)CAT"
stop5 = r"(?<=TTA)(.*)CAT"
stop6 = r"(?<=TCA)(.*)CAT"

regexList = [regex1, regex2, regex3, regex4, regex5, regex6]
stopList = [stop1, stop2, stop3, stop4, stop5, stop6]

for x in tqdm(range(len(split))):
    locus_number = 0
    title = split[x].split(" ")[0]
    fasta_split = split[x].split("\n")
    fasta_header = ">" + fasta_split[0]
    sequence = "".join(fasta_split[1:])
    sequences_all_regions.append(fasta_header)
    sequences_all_regions.append(sequence)
    #region_titles.append("##sequence-region " + title + " 1 " + str(len(sequence)))
    annotations = set()
    for regex in regexList:
            if re.search(regex, sequence):
                start_list = re.finditer(regex,sequence)
            for stop_regex in stopList:
                if re.search(stop_regex, sequence):
                    stop_codon_list = re.finditer(stop_regex,sequence)
                    
                    for result in start_list:
                        for stop_result in stop_codon_list:
                            if (stop_result.end() - result.start()) % 3 == 0:
                                codon = regex.split("(")[0]
                                start_codon = regex[regex.find("=")+1:regex.find(")")]
                                print(start_codon)
                                if start_codon == "ATG":
                                    gene_strand = "+"
                                elif start_codon == "CTA" or start_codon == "TTA" or start_codon == "TCA":
                                    gene_strand = "-"
                                else:
                                    gene_strand = ""
                                annotation = title + "\tCRUDE\tCDS\t" + str(result.start() - 2) + "\t" + str(stop_result.end() + 3) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                                annotations.add(annotation)
                                locus_number += 1
                        
                        
    annotation_all_regions += ["##sequence-region " + title + " 1 " + str(len(sequence))]
    annotation_all_regions += sorted(list(annotations), key=lambda x: x.split('\t')[3])


for x in tqdm(range(len(split))):
    locus_number = 0
    title = split[x].split(" ")[0]
    fasta_split = split[x].split("\n")
    fasta_header = ">" + fasta_split[0]
    sequence = "".join(fasta_split[1:])
    sequences_all_regions.append(fasta_header)
    sequences_all_regions.append(sequence)
    #region_titles.append("##sequence-region " + title + " 1 " + str(len(sequence)))
    annotations = set()
    
    for base in range(len(sequence)):
        index = base
        if sequence[index:index+3] == "ATG":
            gene = ''
            end_codon = ["TAG", "TAA", "TGA"]
            gene_strand = "+"
            start = index #start is the index of ATG
            cut = sequence[index:] #ATG to end of string
            for y in range(len(cut)): # for character in the cut
                if (y + 3) % 3 ==0:
                    gene+=cut[y:y+3]
                    if any(cut[y:y+3] == motif for motif in end_codon):
                        annotation = title + "\tCRUDE\tCDS\t" + str(start + 1) + "\t" + str(start + (y+3)) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                        annotations.add(annotation)
                        locus_number += 1
                    else:
                        pass
                        
        elif any(sequence[index:index+3] == rev_stop for rev_stop in ["CTA", "TTA", "TCA"]):
            gene = ''
            gene_strand = "-"
            start = index #start is the index of ATG
            cut = sequence[index:] #ATG to end of string
            for y in range(len(cut)): # for character in the cut
                if (y + 3) % 3 ==0:
                    gene+=cut[y:y+3]
                    if cut[y:y+3] == "CAT":
                        annotation = title + "\tCRUDE\tCDS\t" + str(start + 1) + "\t" + str(start + (y+3)) + "\t.\t" + gene_strand + "\t0\tID=" + title.split(".")[0] + "_" + str(locus_number) + ";product=putative_protein_region"
                        annotations.add(annotation)
                        locus_number += 1
                    else:
                        pass
                    
    annotation_all_regions += ["##sequence-region " + title + " 1 " + str(len(sequence))]
    annotation_all_regions += sorted(list(annotations), key=lambda x: x.split('\t')[3])
        
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
    
end_time = time.time()

print(str(end_time - start_time) + " seconds" )
print(len(total_genes))
        
        
        
