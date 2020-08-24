#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 12:02:40 2020

@author: danielanderson
"""


import pandas as pd 
file = pd.read_csv("reference_all_annotations.csv")

import re 
genes = []
for index, row in tqdm(file.iterrows()):
    if ";product=" in row["attributes"]:
        genes.append(re.search('product=(.*?);', row["attributes"]).group(1))

genes.count("hypothetical protein")


all_gen = []
yes = []
for node in g._node:
    y = g._node[node]
    if "hypothetical protein" in y['description']:
        if not isinstance(y["members"],int):
            all_gen.append(len(y["members"]))
        else:
            all_gen.append(1)
    if y['description'] == "hypothetical protein":
        if not isinstance(y["members"],int):
            yes.append(len(y["members"]))
        else:
            yes.append(1)
            
        
from tqdm import tqdm
seq = []
for index, row in tqdm(file.iterrows()):
    if row["Name in graph"] == "pspA":
        seq.append(row["sequence"])
        
lengths = []
for x in seq:
    lengths.append(len(x))
    

crouch = pd.read_csv("croucher_all_annotations.csv")
from tqdm import tqdm
name = []
l2 = []
for index, row in tqdm(crouch.iterrows()):
    #if row["sequence"] in seq:
    if "pspA" in row["attributes"]:
        name.append(row["attributes"])
        l2.append(len(row["sequence"]))

ps = []
len_ps = []
isol = set()            
for index, row in tqdm(crouch.iterrows()):
    if "pspA" in row["attributes"]:
        ps.append(row["attributes"])
        len_ps.append(len(row["sequence"]))
        isol.add(row["Isolate"])
        
        
import networkx as nx
g = nx.read_gml("final_graph.gml")

pspA = []

all_gen = []
psps = []
for node in g._node:
    y = g._node[node]
    if "psp" in y['name']:
        psps.append((y["name"], y["members"]))
    if not isinstance(y["members"],int):
        all_gen.append(len(y["members"]))
        if not "group_" in y["name"]:
            pspA.append(len(y["members"]))
    else:
        all_gen.append(1)
        if not "group_" in y["name"]:
            pspA.append(1)
        
        
isol = list(isol)

isol_cleaned = set()

for x in isol:
    isol_cleaned.add("_".join(x.split("_")[:2]))
    

file = pd.read_csv("gene_data.csv")

pspA = []
memA = set()
pspC = []
memC = set()
for index, row in tqdm(file.iterrows()):
    if "pspA" in str(row["gene_name"]):
        pspA.append(row["gene_name"])
        memA.add(row["gff_file"])
    if "pspC" in str(row["gene_name"]):
        pspC.append(row["gene_name"])
        memC.add(row["gff_file"])
        

memA = list(memA)
memC = list(memC)

mem = memA + memC
mem = set(mem)



g = nx.read_gml("ref_graph.gml")
g_int = nx.read_gml("int_graph.gml")

all_g = []
for node in g._node:
    y = g._node[node]
    all_g.append(y["name"])

all_int = []
for node in g_int._node:
    y = g_int._node[node]
    all_int.append(y["name"])

not_in = []
for x in all_int:
    if not x in all_g and not "~" in x:
        not_in.append(x)
   
    if "hypothetical protein" in y['description']:
        if not isinstance(y["members"],int):
            all_gen.append(len(y["members"]))
        else:
            all_gen.append(1)
    if y['description'] == "hypothetical protein":
        if not isinstance(y["members"],int):
            yes.append(len(y["members"]))
        else:
            yes.append(1)
            
import glob 
gff = glob.glob("data/comp/*.gff.gz")
import gzip

num = []
for x in gff:
    with gzip.open(x, 'rt') as f:
        stored = str(f.read())
    stored = stored.splitlines()
    for line in stored:
        if not "#" in line and "CDS" in line:
            num.append(line)