#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script was used to produce the results for ROC curves assessing the perfomance of COBS index under different k-mer lengths and query thresholds. 
Indexes can be built from the translated "centroids" associated with each Panaroo node or from a random subset of 100,000 alignment files of all clustered genes.  
A Panaroo-updated GFF file for a reference strain is then provided and searches are conducted for all genes, then core genes and accessory genes separately (based on whether the genes are core or accessory in the Panaroo output)
"""
from random import randint
import cobs_index as cobs
import glob
import pandas as pd
from tqdm import tqdm
import os 
import tempfile
import re
import networkx as nx

def get_characteristic(gene, found, result, partial_result, query, file_dir):
    result = str(result)
    partial_result = str(partial_result)
    if not partial_result == "[]":
        filename = str(partial_result.split("'")[1]) + ".txt"
        with open(file_dir + filename, "r") as s:
            length_seq = len(s.read())
        if len(query) < (length_seq/2):
            partial_result = "[]"
        else:
            partial_result = partial_result
    if gene in partial_result and found == 1:
        value = "TP"
    elif not gene in result and partial_result == "[]":
        value = "TN"
    elif not partial_result == "[]" and found == 0:
        value = "FP"
    elif partial_result == "[]" and gene in result:
        value = "FN"
    else:
        value = ""
    return value
    
def create_indexes_from_kmers(directory):
    test_kmers = [1,2,4,6,8,10]
    for kmer in tqdm(test_kmers):
        params = cobs.CompactIndexParameters()
        params.term_size = kmer
        params.clobber = True               # overwrite output and temporary files
        params.false_positive_rate = 0.01    # higher false positive rate -> smaller index
        cobs.compact_construct(directory,"ini_cluster_indexes/" + str(kmer) + "_index.cobs_compact", index_params=params)
    return test_kmers

def search_index(kmer, sequence, name, file_dir):
    found = ''
    index = cobs.Search("ini_cluster_indexes/" + str(kmer) + "_index.cobs_compact")
    result_zero= index.search(sequence, threshold = 0.0)
    result_zero_five= index.search(sequence, threshold = 0.05)
    result_one= index.search(sequence, threshold = 0.1)
    result_one_five= index.search(sequence, threshold = 0.15)
    result_two= index.search(sequence, threshold = 0.2)
    result_two_five= index.search(sequence, threshold = 0.25)
    result_three= index.search(sequence, threshold = 0.3)
    result_three_five= index.search(sequence, threshold = 0.35)
    result_four= index.search(sequence, threshold = 0.4)
    result_four_five= index.search(sequence, threshold = 0.45)
    result_five= index.search(sequence, threshold = 0.5)
    result_five_five= index.search(sequence, threshold = 0.55)
    result_six = index.search(sequence, threshold = 0.6)
    result_six_five= index.search(sequence, threshold = 0.65)
    result_seven= index.search(sequence, threshold = 0.7)
    result_seven_five= index.search(sequence, threshold = 0.75)
    result_eight=index.search(sequence, threshold = 0.8)
    result_eight_five= index.search(sequence, threshold = 0.85)
    result_nine= index.search(sequence, threshold = 0.9)
    result_nine_five= index.search(sequence, threshold = 0.95)
    result_ten= index.search(sequence, threshold = 1.0)
    
    if name == (result_zero[0][1]).split(";")[0]:
        found = 1
    else:
        found = 0
    
    result_zero_five = get_characteristic(name, found, result_zero, result_zero_five, sequence, file_dir)
    result_one = get_characteristic(name, found, result_zero, result_one, sequence, file_dir)
    result_one_five = get_characteristic(name, found, result_zero, result_one_five, sequence, file_dir)
    result_two = get_characteristic(name, found, result_zero, result_two, sequence, file_dir)
    result_two_five = get_characteristic(name, found, result_zero, result_two_five, sequence, file_dir)
    result_three = get_characteristic(name, found, result_zero, result_three, sequence, file_dir)
    result_three_five = get_characteristic(name, found, result_zero, result_three_five, sequence, file_dir)
    result_four = get_characteristic(name, found, result_zero, result_four, sequence, file_dir)
    result_four_five = get_characteristic(name, found, result_zero, result_four_five, sequence, file_dir)
    result_five = get_characteristic(name, found, result_zero, result_five, sequence, file_dir)
    result_five_five = get_characteristic(name, found, result_zero, result_five_five, sequence, file_dir)
    result_six = get_characteristic(name, found, result_zero, result_six, sequence, file_dir)
    result_six_five = get_characteristic(name, found, result_zero, result_six_five, sequence, file_dir)
    result_seven = get_characteristic(name, found, result_zero, result_seven, sequence, file_dir)
    result_seven_five = get_characteristic(name, found, result_zero, result_seven_five, sequence, file_dir)
    result_eight = get_characteristic(name, found, result_zero, result_eight, sequence, file_dir)
    result_eight_five = get_characteristic(name, found, result_zero, result_eight_five, sequence, file_dir)
    result_nine = get_characteristic(name, found, result_zero, result_nine, sequence, file_dir)
    result_nine_five = get_characteristic(name, found, result_zero, result_nine_five, sequence, file_dir)
    result_ten = get_characteristic(name, found, result_zero, result_ten, sequence, file_dir)
    result_zero = get_characteristic(name, found, result_zero, result_zero, sequence, file_dir)

    return (found, [result_zero, result_zero_five, result_one, result_one_five, result_two, result_two_five, result_three, result_three_five, result_four, result_four_five, result_five, result_five_five, result_six, result_six_five, result_seven, result_seven_five, result_eight, result_eight_five, result_nine, result_nine_five, result_ten])

def get_found(in_tuple):
    return in_tuple[0]

def result_zero(in_tuple):
    return in_tuple[1][0]
def result_zero_five(in_tuple):
    return in_tuple[1][1]
def result_one(in_tuple):
    return in_tuple[1][2]
def result_one_five(in_tuple):
    return in_tuple[1][3]
def result_two(in_tuple):
    return in_tuple[1][4]
def result_two_five(in_tuple):
    return in_tuple[1][5]
def result_three(in_tuple):
    return in_tuple[1][6]
def result_three_five(in_tuple):
    return in_tuple[1][7]
def result_four(in_tuple):
    return in_tuple[1][8]
def result_four_five(in_tuple):
    return in_tuple[1][9]
def result_five(in_tuple):
    return in_tuple[1][10]
def result_five_five(in_tuple):
    return in_tuple[1][11]
def result_six(in_tuple):
    return in_tuple[1][12]
def result_six_five(in_tuple):
    return in_tuple[1][13]
def result_seven(in_tuple):
    return in_tuple[1][14]
def result_seven_five(in_tuple):
    return in_tuple[1][15]
def result_eight(in_tuple):
    return in_tuple[1][16]
def result_eight_five(in_tuple):
    return in_tuple[1][17]
def result_nine(in_tuple):
    return in_tuple[1][18]
def result_nine_five(in_tuple):
    return in_tuple[1][19]
def result_ten(in_tuple):
    return in_tuple[1][20]
    
def search_for_query(kmer_lengths, protein_sequence, protein_name, file_dir):
    search_df = pd.DataFrame(kmer_lengths, columns = ["Kmer Length"])
    search_df["Found and result"] = search_df.apply(lambda row: search_index(row["Kmer Length"], protein_sequence, protein_name, file_dir), axis =1)
    search_df["Found"] = search_df["Found and result"].apply(get_found)
    search_df["Result_0"] = search_df["Found and result"].apply(result_zero)
    search_df["Result_05"] = search_df["Found and result"].apply(result_zero_five)
    search_df["Result_1"] = search_df["Found and result"].apply(result_one)
    search_df["Result_15"] = search_df["Found and result"].apply(result_one_five)
    search_df["Result_2"] = search_df["Found and result"].apply(result_two)
    search_df["Result_25"] = search_df["Found and result"].apply(result_two_five)
    search_df["Result_3"] = search_df["Found and result"].apply(result_three)
    search_df["Result_35"] = search_df["Found and result"].apply(result_three_five)
    search_df["Result_4"] = search_df["Found and result"].apply(result_four)
    search_df["Result_45"] = search_df["Found and result"].apply(result_four_five)
    search_df["Result_5"] = search_df["Found and result"].apply(result_five)
    search_df["Result_55"] = search_df["Found and result"].apply(result_five_five)
    search_df["Result_6"] = search_df["Found and result"].apply(result_six)
    search_df["Result_65"] = search_df["Found and result"].apply(result_six_five)
    search_df["Result_7"] = search_df["Found and result"].apply(result_seven)
    search_df["Result_75"] = search_df["Found and result"].apply(result_seven_five)
    search_df["Result_8"] = search_df["Found and result"].apply(result_eight)
    search_df["Result_85"] = search_df["Found and result"].apply(result_eight_five)
    search_df["Result_9"] = search_df["Found and result"].apply(result_nine)
    search_df["Result_95"] = search_df["Found and result"].apply(result_nine_five)
    search_df["Result_10"] = search_df["Found and result"].apply(result_ten)
    search_df["Gene Name"] = protein_name
    search_df = search_df[["Gene Name", "Kmer Length", "Found", "Result_0", "Result_05","Result_1","Result_15", "Result_2", "Result_25","Result_3", "Result_35","Result_4", "Result_45","Result_5", "Result_55", "Result_6", "Result_65","Result_7", "Result_75","Result_8", "Result_85","Result_9", "Result_95","Result_10"]]
    return search_df

def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary
    """"Translate nucleotide sequences"""
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    protein =""
    try:
        if len(seq)%3 == 0:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3] #Defining a codon as 3 bases
                protein+= table[codon] #Translate the codon into an amino acid based on the dictionary and append this to the protein sequence
    except:
        protein = ""
    return protein

def translate_alignements(alignement_path):
    """"Translate all alignement files in path and add to single list"""
    core_names = []
    accessory_names = []
    all_proteins = []
    for file in tqdm(alignement_path):
        with open(file, "r") as f:
            aln = f.read()
        gene_name = (file.split("/")[2]).split(".")[0]
        aln = aln.split(">")[1:]
        cleaned = []
        for x in aln:
           sequence = "".join(x.splitlines()[1:]).replace("-", "")
           sequence = ">"+gene_name + "\n" + translate(sequence.upper())
           if not sequence == (">"+gene_name + "\n"):
               cleaned.append(sequence)
               
        if len(cleaned) > 609:
            core_names.append(gene_name)
        else:
            accessory_names.append(gene_name)
        
        cleaned = "\n".join(cleaned)
        all_proteins.append(cleaned)
    
    all_proteins = "\n".join(all_proteins)
    return core_names, accessory_names, all_proteins

def long_centroid(graph):
    """"Translate all alignement files in path and add to single list"""
    core_names = []
    accessory_names = []
    all_proteins = []
    
    G = nx.read_gml(graph)
    for node in G._node:
        y = G._node[node]
        seq = y["protein"].split(";")
        if not isinstance(y["members"], int):
            if len(y["members"]) > 609:
                core_names.append(y["name"])
            else:
                accessory_names.append(y["name"])
        else:
            accessory_names.append(y["name"])
        for x in seq:
            all_proteins.append(">" + y["name"] + "\n" + x)
            
    all_proteins = "\n".join(all_proteins)
    return core_names, accessory_names, all_proteins

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'R': 'R', 'K':'K'}
    return ''.join([complement[base] for base in dna[::-1]])

def subset_attributes(subset_gffs_path, all_core_names, all_accessory_names):
    
    subset_core_names = []
    subset_core_sequences = []
    subset_accessory_names = []
    subset_accessory_sequences = []
    
    unannotated_name_list = []
    unannotated_sequence_list = []
    core_unannotated_name_list = []
    core_unannotated_sequence_list = []
    accessory_unannotated_sequence_list = []
    accessory_unannotated_name_list = []

    for file in tqdm(subset_gffs_path):
        with open(file, "r") as s:
            gff_file = s.read().split("##FASTA")
        gff_info = gff_file[0].splitlines()
        fasta_info = "".join((gff_file[1].splitlines())[2:])
    
        for x in gff_info:
            if not "##sequence-region" in x:
                split = x.split("\t")
                if re.search("gene=(.*);", split[8]):
                    gene_name = re.search("gene=(.*);", split[8]).group(1)
                    gene_name = str(gene_name.split(";")[0])
                    if gene_name in all_core_names:
                        subset_core_names.append(gene_name)
                        if not split[6] == "-":
                            sequence = translate(fasta_info[int(split[3])-1:int(split[4])])
                            subset_core_sequences.append(sequence)
                        elif split[6] == "-":
                            sequence = translate(reverse_complement(fasta_info[int(split[3])-1:int(split[4])]))
                            subset_core_sequences.append(sequence)
                    elif gene_name in all_accessory_names:
                        subset_accessory_names.append(gene_name)
                        if not split[6] == "-":
                            sequence = translate(fasta_info[int(split[3])-1:int(split[4])])
                            subset_accessory_sequences.append(sequence)
                        elif split[6] == "-":
                            sequence = translate(reverse_complement(fasta_info[int(split[3])-1:int(split[4])]))
                            subset_accessory_sequences.append(sequence)
                else:
                    unannotated_name_list.append("Unannotated")
                    if not split[6] == "-":
                        unannotated_sequence = translate(fasta_info[int(split[3])-1:int(split[4])])
                        unannotated_sequence_list.append(unannotated_sequence)
                    elif split[6] == "-":
                        unannotated_sequence = translate(reverse_complement(fasta_info[int(split[3])-1:int(split[4])]))
                        unannotated_sequence_list.append(unannotated_sequence)
    
    core_unannotated_name_list = unannotated_name_list[:len(subset_core_names)]
    core_unannotated_sequence_list = unannotated_sequence_list[:len(subset_core_names)]
    accessory_unannotated_name_list = unannotated_name_list[:len(subset_accessory_names)]
    accessory_unannotated_sequence_list = unannotated_sequence_list[:len(subset_accessory_names)]
    
    subset_core_names += core_unannotated_name_list
    subset_core_sequences += core_unannotated_sequence_list
    subset_accessory_names += accessory_unannotated_name_list
    subset_accessory_sequences += accessory_unannotated_sequence_list
    return subset_core_names, subset_core_sequences, subset_accessory_names, subset_accessory_sequences
        
    
alignement_path = glob.glob("reference_pangenome/aligned_gene_sequences/*.aln.fas")
alignement_path += glob.glob("reference_pangenome/aligned_gene_sequences/*.fasta")

#all_core_names, all_accessory_names, all_proteins = translate_alignements(alignement_path)
all_core_names, all_accessory_names, all_proteins = long_centroid("reference_pangenome/final_graph.gml")

os.mkdir("cluster_files")
os.mkdir("cluster_indexes")
os.mkdir("cluster_results")

protein_file = all_proteins.split(">")[1:]

titles = []
cleaned = []
for x in protein_file:
    cleaned.append("".join(x.splitlines()[1:]))
    titles.append(x.splitlines()[0])

#hundred_thousand_indexes = [randint(0,(len(cleaned) -1)) for i in range(100000)]
hundred_thousand_indexes = []
for x in range(len(protein_file)):
    hundred_thousand_indexes.append(x)
    
proteins = []
gene_names = []
for protein in hundred_thousand_indexes:
    info = cleaned[protein]
    proteins.append(info)
    gene_names.append(titles[protein])

temp_dir = os.path.join(tempfile.mkdtemp(dir="cluster_files"), "")

duplicate_set = set()
duplicate_count = 1
for prot_index in tqdm(range(len(hundred_thousand_indexes))):
    info = proteins[prot_index]
    if not gene_names[prot_index] in duplicate_set:
        thou = open(temp_dir + gene_names[prot_index] + ".txt", "w")
        thou.write(info)
        thou.close()
        duplicate_set.add(gene_names[prot_index])
    else:
        thou = open(temp_dir + gene_names[prot_index] + ";" + str(duplicate_count) + ".txt", "w")
        thou.write(info)
        thou.close()
        duplicate_count += 1

kmer_lengths = create_indexes_from_kmers(temp_dir)

subset_gffs_path = ["crude_annotated/GCF_000006885.1_ASM688v1_mk2.gff"]
subset_core_names, subset_core_sequences, subset_accessory_names, subset_accessory_sequences = subset_attributes(subset_gffs_path, all_core_names, all_accessory_names)

all_names = subset_core_names + subset_accessory_names
all_sequences = subset_core_sequences + subset_accessory_sequences
 
search_names = [subset_accessory_names, subset_core_names, all_names]
search_sequences = [subset_accessory_sequences, subset_core_sequences, all_sequences]

for name_list in range(len(search_names)):
    temp_results = os.path.join(tempfile.mkdtemp(dir="cluster_results"), "")
    for query in tqdm(range(len(search_names[name_list]))):
        query_prot_sequence = search_sequences[name_list][query]
        query_prot_name = search_names[name_list][query]
        kmer_df = search_for_query(kmer_lengths, query_prot_sequence, query_prot_name, temp_dir)
        kmer_df.to_csv(temp_results + str(query) + ".csv")
        
print("Queries complete")

