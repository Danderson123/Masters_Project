#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 23:26:18 2020

@author: danielanderson
"""

import pandas as pd
from tqdm import tqdm

def con(pan, cob):
    if not pan == "group_1890" and not pan == "group_1938":
        if pan == cob:
            value=pan
        else:
            value = 0
    elif pan == "group_1890" and cob == "group_2880":
        value = pan
    elif pan == "group_1938" and cob == "group_2939":
        value = pan
    else:
        value = 0
    return value

eight = pd.read_csv("85_SEARCHED_all_annotations.csv")

tqdm.pandas()
eight["con"] = eight.progress_apply(lambda row: con(row["Name in graph"], row["COBS 85%"]), axis =1)

isolates = []
core = []
cobs = []
pan_freqs = []
consensus = []
onea_names = []
onea_isols = set()
twoa_names = []
twoa_isols = set()
secy2 = set()
con_secy2 = set()
con_onea = set()
con_twoa = set()

for index, row in tqdm(eight.iterrows()):
    #if row["CorA"] == "core":
    pan_freqs.append(row["Frequency (%)"])
    core.append(row["Name in graph"])
    #cobs.append(row["Isolate"])
    cobs.append(row["COBS 85%"])
    consensus.append(row["con"])
    if row["Name in graph"] == "group_1890":
        onea_names.append(row["COBS 85%"])
    if row["COBS 85%"] == "group_2880":
        onea_isols.add(row["Isolate"])
    if row["Name in graph"] == "group_1938":
        twoa_names.append(row["COBS 85%"]) 
    if row["COBS 85%"] == "group_2939":
        twoa_isols.add(row["Isolate"])
    if row["COBS 85%"] == "secY2":
        secy2.add(row["Isolate"])
    if row["con"] == "secY2":
        con_secy2.add(row["Isolate"])
    if row["con"] == "group_1890":
       con_onea.add(row["Isolate"])
    if row["con"] == "group_1938":
       con_twoa.add(row["Isolate"])

       # cobs_freqs.append(ls.count(row["COBS 85%"])/616*100)
        
        
pbponea = list(eight["Name in graph"]).count("secY2")
pbptwoa = list(eight["COBS 85%"]).count("group_2939")

import networkx as nx

g = nx.read_gml("final_graph.gml")

unique = []
desc = []
mem = []
nam = []
for node in g._node:
    y = g._node[node]
    if "enicillin" in y["description"]:
        desc.append(y["description"])
        mem.append(y["members"])
        nam.append(y["name"])

c = nx.read_gml("crouch_graph.gml")

c_unique = []
c_desc = []
c_mem = []
c_nam = []
for node in c._node:
    y = c._node[node]
    if "enicillin" in y["description"]:
        c_desc.append(y["description"])
        c_mem.append(y["members"])
        c_nam.append(y["name"])


    
core_set = set(core)
core_set = list(core_set)    

pan_freqs = []
for x in core_set:
    for y, row in eight.iterrows():
        pan_freqs
        