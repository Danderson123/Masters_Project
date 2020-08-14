#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Search for AMR/virulent genes in two graphs and compare their frequencies
"""


import networkx as nx

G = nx.read_gml("final_graph.gml")

erm = []
G_names = []
zp = []

for node in G._node:
    y = G._node[node]
    print(len([y["members"]]))
    if len([y["members"]]) > 26:
        G_names.append(y["name"])
    if y["name"] == "tetM":
        tetm = y["members"]
    if y["name"] == "ermB":
        ermb = y["members"]
    if "erm" in y["name"]:
       erm.append([y["name"], y["members"]])
    if y["name"] == "mefA":
        mefa = y["members"]  
    if y["name"] == "apha3IIIa":
       apha3IIIa = y["members"] 
    if "aph" in y["name"]:
       apha3IIIa = y["members"] 
    if "zmp" in y["name"]:
       zmp = y["members"]
       zp.append([y["name"], y["members"]])
       
       
       
G_crouch = nx.read_gml("croucher_merged.gml")

G_names = []
erm = []
zp = []
pbp = []

for node in G_crouch._node:
    y = G_crouch._node[node]
    G_names.append(y["name"])
    if y["name"] == "tetM":
        tetm = y["members"][10:]
    if y["name"] == "ermB":
        ermb = y["members"]
    if "erm" in y["name"]:
       erm.append([y["name"], y["members"]])
    if y["name"] == "mefA":
        mefa = y["members"]
    if y["name"] == "apha3IIIa":
       apha3IIIa = y["members"] 
    if "aph" in y["name"]:
        apha3IIIa = y["members"] 
    if "zmp" in y["name"]:
       zmp = y["members"]
       zp.append([y["name"], y["members"]])
    if "pbp" in y["name"]:
       pbp = y["members"]
       pbp.append([y["name"], y["members"]])
       
count = 0
for x in G._node:
    if "group" in G._node[x]["name"]:
        count += 1
        
        45/616*100
