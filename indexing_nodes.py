#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script was used to determine whether the panaroo node indexes remain consistent when merging graphs. They do not.
"""



import networkx as nx

G_original = nx.read_gml('not_split/final_graph.gml')

Gene_indexes= []
Gene_name = []
for node in G_original._node:
    Gene_indexes.append(node)
    Gene_name.append(G_original._node[node]["name"])
    
original = dict(zip(Gene_indexes, Gene_name))

G_merged = nx.read_gml('merged/final_graph.gml')

Merged_indexes= []
Merged_name = []
for node in G_merged._node:
    Merged_indexes.append(node)
    Merged_name.append(G_merged._node[node]["name"])

merged = dict(zip(Merged_indexes, Merged_name))

