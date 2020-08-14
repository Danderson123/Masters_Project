#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This uses the gene_presence_absence output of panaroo to associate all input annotations with their assigned name in the pangenomic graph. 
This is important because node organisation changes whenever you add a genome to the graph. 
It takes two CSV files containing all of the features inputted into two panaroo graphs that have been subsequently merged
""" 
import pandas as pd
from io import StringIO
from tqdm import tqdm

def CorA(freq):
    if freq >= 99:
        value = "core"
    else:
        value = "accessory"
    return value
    
all_annotations = pd.read_csv("REAL_prodigal_unannotated_SENSITIVE_pangenome/all_annotations.csv")
reference_annotations = pd.read_csv("REFERENCE_GFFS/all_annotations.csv")

with open("reference_sparc_merged/gene_presence_absence.csv", "r") as p:
    correct_headers = p.read()

split_file = correct_headers.splitlines()
ID_length = list(range(1, len(split_file[1].split(",")) - 2))
New_headers = split_file[0] + "," + str(ID_length).replace("[", "").replace("]","").replace(" ", "")
information = "\n".join([New_headers] + split_file[1:])

data = StringIO(information)
presence_absence = pd.read_csv(data, sep=",")

IDs = presence_absence.iloc[:, 3:]

IDs = list(pd.Series(IDs.fillna('').values.tolist()).str.join(''))
    
tqdm.pandas()

def get_gene_name(ID):
    index = ([i for i, x in enumerate(IDs) if ID in x])
    if not len(index) == 0:
        index = presence_absence['Gene'][index[0]]
    else:
        index = float("NaN")
    return index

all_annotations = all_annotations.append(reference_annotations, ignore_index=True)
all_annotations["Name in graph"] = all_annotations["ID"].progress_apply(get_gene_name)

all_annotations.dropna(subset = ["Name in graph"], inplace=True)
all_annotations = all_annotations.sort_values(by='Name in graph')
all_annotations = all_annotations.reset_index(drop=True)


gene_names = []
frequencies = []

import networkx as nx
graph = nx.read_gml("reference_sparc_merged/final_graph.gml")
for value in graph.graph.values():
    num_isolates = len(value)

for node in tqdm(graph._node):
    y = graph._node[node]
    gene_names.append(y["name"].lower())
    num_sequences = y["seqIDs"]
    unique = set()
    for x in range(len(num_sequences)):
        unique.add(num_sequences[x].split("_")[0])
    frequency = (len(unique)/num_isolates) * 100
    frequencies.append(frequency)

gene_table = pd.DataFrame(gene_names, columns = ["Gene Name"])
gene_table["Frequency (%)"] = frequencies
gene_table = gene_table.drop_duplicates(subset=['Gene Name'])

gene_table = gene_table.sort_values(by='Gene Name')

def get_freq(name):
    name = name.lower()
    freq = gene_table["Frequency (%)"][list(gene_table["Gene Name"]).index(name)]
    return freq

ID_list = []
number = 0
for x in all_annotations["Name in graph"].index:
    ID_list.append("PN_" + str(number))
    number += 1
all_annotations["New ID"] = ID_list
all_annotations.to_csv("reference_sparc_merged/all_annotations.csv", index = False)

all_annotations["Frequency (%)"] = all_annotations["Name in graph"].apply(get_freq)
all_annotations["CorA"] = all_annotations["Frequency (%)"].apply(CorA)
all_annotations.to_csv("reference_sparc_merged/all_annotations.csv", index = False)
