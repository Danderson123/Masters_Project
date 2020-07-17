#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script was used to generate plots for figure ..... It is used to determine the sensitivity and specificity for different: kmer lengths and number of files to determine kmer length.
"""
from random import randint
import cobs_index as cobs
import networkx as nx
import pandas as pd
from tqdm import tqdm
import os 
import tempfile
import shutil

def get_cluster_name(cluster):
	name = None
	for node in graph._node:
		y = graph._node[node]
		if cluster in y["geneIDs"]:
			name = y["name"].lower()	
	return name

def create_indexes_from_kmers(directory):
    test_kmers = [1,2,5,10,15,20]
    for kmer in test_kmers:
        params = cobs.CompactIndexParameters()
        params.term_size = kmer
        params.clobber = True               # overwrite output and temporary files
        params.false_positive_rate = 0.4    # higher false positive rate -> smaller index
        cobs.compact_construct(directory,"indexes/" + str(kmer) + "_index.cobs_compact", index_params=params)
    return test_kmers

def search_index(kmer, sequence, name):
    found = ''
    index = cobs.Search("indexes/" + str(kmer) + "_index.cobs_compact")
    result_zero= index.search(sequence, threshold = 0.0)
    result_one= index.search(sequence, threshold = 0.1)
    result_two= index.search(sequence, threshold = 0.2)
    result_three= index.search(sequence, threshold = 0.3)
    result_four= index.search(sequence, threshold = 0.4)
    result_five= index.search(sequence, threshold = 0.5)
    result_six = index.search(sequence, threshold = 0.6)
    result_seven= index.search(sequence, threshold = 0.7)
    result_eight=index.search(sequence, threshold = 0.8)
    result_nine= index.search(sequence, threshold = 0.9)
    result_ten= index.search(sequence, threshold = 1.0)
    if name == (result_zero[0][1]).split(";")[0]:
        found = 1
    else:
        found = 0
    return (found, [result_zero, result_one, result_two, result_three, result_four, result_five, result_six, result_seven, result_eight, result_nine, result_ten])

def get_found(in_tuple):
    return in_tuple[0]

def result_zero(in_tuple):
    return in_tuple[1][0]
def result_one(in_tuple):
    return in_tuple[1][1]
def result_two(in_tuple):
    return in_tuple[1][2]
def result_three(in_tuple):
    return in_tuple[1][3]
def result_four(in_tuple):
    return in_tuple[1][4]
def result_five(in_tuple):
    return in_tuple[1][5]
def result_six(in_tuple):
    return in_tuple[1][6]
def result_seven(in_tuple):
    return in_tuple[1][7]
def result_eight(in_tuple):
    return in_tuple[1][8]
def result_nine(in_tuple):
    return in_tuple[1][9]
def result_ten(in_tuple):
    return in_tuple[1][10]
    
def search_for_query(kmer_lengths, protein_sequence, protein_name):
    search_df = pd.DataFrame(kmer_lengths, columns = ["Kmer Length"])
    search_df["Found and result"] = search_df.apply(lambda row: search_index(row["Kmer Length"], protein_sequence, protein_name), axis =1)     
    search_df["Found"] = search_df["Found and result"].apply(get_found)
    search_df["Result_0"] = search_df["Found and result"].apply(result_zero)
    search_df["Result_1"] = search_df["Found and result"].apply(result_one)
    search_df["Result_2"] = search_df["Found and result"].apply(result_two)
    search_df["Result_3"] = search_df["Found and result"].apply(result_three)
    search_df["Result_4"] = search_df["Found and result"].apply(result_four)
    search_df["Result_5"] = search_df["Found and result"].apply(result_five)
    search_df["Result_6"] = search_df["Found and result"].apply(result_six)
    search_df["Result_7"] = search_df["Found and result"].apply(result_seven)
    search_df["Result_8"] = search_df["Found and result"].apply(result_eight)
    search_df["Result_9"] = search_df["Found and result"].apply(result_nine)
    search_df["Result_10"] = search_df["Found and result"].apply(result_ten)
    search_df["Gene Name"] = protein_name
    search_df = search_df[["Gene Name", "Kmer Length", "Found", "Result_0", "Result_1", "Result_2", "Result_3", "Result_4", "Result_5", "Result_6", "Result_7", "Result_8", "Result_9", "Result_10"]]
    return search_df

with open("reference_pangenome/combined_protein_cdhit_out.txt", "r") as f:
	protein_file = f.read().split(">")[1:]

graph = nx.read_gml("reference_pangenome/final_graph.gml")

titles = []
cleaned = []
for x in protein_file:
	cleaned.append("".join(x.splitlines()[1:]))
	titles.append(x.splitlines()[0])

tqdm.pandas()

cluster_name = pd.DataFrame(titles, columns = ["Cluster"])
cluster_name["gene"] = cluster_name["Cluster"].progress_apply(get_cluster_name)
cluster_name = cluster_name.dropna()

cleaned = [cleaned[i] for i in list(cluster_name.index)]
cluster_name = cluster_name.reset_index(drop=True)

#The following randomly selects subsets of documents. This will later be changed so that multiple subsets are created to total 1000 documents for each.
onethousand_indexes = [randint(0,(len(cleaned) -1)) for i in range(1000)]

os.mkdir("files")
os.mkdir("indexes")
#os.mkdir("results")

proteins = []
gene_names = []

filename_set = set()
print("Separating proteins...")
for protein in onethousand_indexes:
    info = cleaned[protein]
    proteins.append(info)
    gene_names.append(cluster_name["gene"][protein])
    
ten_indexes_filenames = [onethousand_indexes[i::100] for i in range(100)]
onehundred_indexes_filenames = [onethousand_indexes[i::10] for i in range(10)]
onethousand_indexes_filenames = [onethousand_indexes]

#filename_list = [ten_indexes_filenames, onehundred_indexes_filenames, onethousand_indexes_filenames]
filename_list = [onethousand_indexes_filenames]

print("Generating and searching indexes...")
for file in filename_list:
    count = 0
    for a in range(len(file)):
        for c in tqdm(range(len(file[a]))):
            temp_dir = os.path.join(tempfile.mkdtemp(dir="files"), "") 
            without_y = file[a][:]
            del without_y[c]

            duplicate_set = set()
            duplicate_count = 1
            for prot_index in without_y:
                info = cleaned[prot_index]
                if not cluster_name["gene"][protein] in duplicate_set:
                    thou = open(temp_dir + cluster_name["gene"][prot_index] + ".txt", "w")
                    thou.write(info)
                    thou.close()
                    duplicate_set.add(cluster_name["gene"][prot_index])
                else:
                    thou = open(temp_dir + cluster_name["gene"][prot_index] + ";" + str(count) + ".txt", "w")
                    thou.write(info)
                    thou.close()
                    duplicate_count += 1
            kmer_lengths = create_indexes_from_kmers(temp_dir)
    
            query_prot_sequence = proteins[c]
            query_prot_name = gene_names[c]

            kmer_df = search_for_query(kmer_lengths, query_prot_sequence, query_prot_name)
            kmer_df.to_csv("results/" + str(count) + "_" + str(len(file[0])) + ".csv")
            count += 1
            shutil.rmtree(temp_dir)


