#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 11:11:18 2020

@author: danielanderson
"""


import networkx as nx

G = nx.read_gml("final_graph.gml")
Gc = nx.read_gml("crouch_graph.gml")

g_names = []
c_names = []
g_core = []
c_core = []

g_all = []
c_all = []
g_named = []
c_named = []
g_name_core = []
c_name_core = []
for node in G._node:
    y = G._node[node]
    g_all.append(y["name"])
    if not "group_" in y["name"]:
        g_named.append(y["name"])
        if not isinstance(y["members"], int):
            if len(y["members"]) > 637:
               g_name_core.append(y["name"])
    if not y["description"] == "" and not y["description"] == "hypothetical protein":
        g_names.append(y["name"])
        if not isinstance(y["members"], int):
            if len(y["members"]) > 637:
               g_core.append(y["name"])
               
for node in Gc._node:
    y = Gc._node[node]
    c_all.append(y["name"])
    if not "group_" in y["name"]:
        c_named.append(y["name"])
        if not isinstance(y["members"], int):
            if len(y["members"]) > 637:
               c_name_core.append(y["name"])
    if not y["description"] == "" and not y["description"] == "hypothetical protein":
        c_names.append(y["name"])
        if not isinstance(y["members"], int):
            if len(y["members"]) > 637:
               c_core.append(y["name"])
num = []
length = []

for x in g_names:
    try:
        with open("aligned_gene_sequences/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
       # with open("aligned_gene_sequences/" + x + ".fasta", "r") as f:
           # alns = f.read()
        continue
    alns = alns.split(">")[1:]
    num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    length.append(sum(sub_len))

average_len = sum(length)/len(length)
avergae_num = sum(num)/len(num)

c_num = []
c_length = []

for x in c_names:
    try:
        with open("croucher_alignments/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
        #with open("croucher_alignments/" + x + ".fasta", "r") as f:
            #alns = f.read()
        continue
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    c_length.append(sum(sub_len))

c_average_len = sum(c_length)/len(c_length)
c_avergae_num = sum(c_num)/len(c_num)


g_core_length = []
c_core_length = []

for x in g_core:
    with open("aligned_gene_sequences/" + x + ".aln.fas", "r") as f:
        alns = f.read()
    alns = alns.split(">")[1:]
    num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    g_core_length.append(sum(sub_len))

for x in c_core:
    try:
        with open("croucher_alignments/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
        with open("croucher_alignments/" + x + ".fasta", "r") as f:
            alns = f.read()
        continue
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    c_core_length.append(sum(sub_len))


g_all_seq = []
for x in g_all:
    try:
        with open("aligned_gene_sequences/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
        with open("aligned_gene_sequences/" + x + ".fasta", "r") as f:
            alns = f.read()
        continue
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    g_all_seq.append(sum(sub_len))
    

c_all_seq = []
for x in c_all:
    try:
        with open("croucher_alignments/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
        with open("croucher_alignments/" + x + ".fasta", "r") as f:
            alns = f.read()
        continue
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    c_all_seq.append(sum(sub_len))
    
annot_prop_g = 1122359362/1184477453
annot_prop_c = 1066334521/1116592720

core_prop_g = 892818303/1184477453
core_prop_c = 815269080/1116592720

accessory_prop_g = 229541059/1184477453
accessory_prop_c = 251065441/1116592720

g_named = []
for x in g_all:
    if not "group_" in x:
        g_named.append(x)
c_named = []
for x in c_all:
    if not "group_" in x:
        c_named.append(x)
        
g_name_seq = []
for x in g_named:
    try:
        with open("aligned_gene_sequences/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
        with open("aligned_gene_sequences/" + x + ".fasta", "r") as f:
            alns = f.read()
        continue
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    g_name_seq.append(sum(sub_len))
    

c_name_seq = []
for x in c_named:
    try:
        with open("croucher_alignments/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
        with open("croucher_alignments/" + x + ".fasta", "r") as f:
            alns = f.read()
        continue
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    c_name_seq.append(sum(sub_len))

g_name_core_seq = []
for x in g_name_core:
    with open("aligned_gene_sequences/" + x + ".aln.fas", "r") as f:
        alns = f.read()
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    g_name_core_seq.append(sum(sub_len))
    

c_name_core_seq = []
for x in c_name_core:
    try:
        with open("croucher_alignments/" + x + ".aln.fas", "r") as f:
            alns = f.read()
    except:
        with open("croucher_alignments/" + x + ".fasta", "r") as f:
            alns = f.read()
        continue
    alns = alns.split(">")[1:]
    c_num.append(len(alns))
    
    sub_len = []
    for line in alns:
        line = "".join(line.splitlines()[1:])
        line = line.replace("-","")
        sub_len.append(len(line))
    c_name_core_seq.append(sum(sub_len))