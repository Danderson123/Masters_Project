#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 14:13:12 2020

@author: danielanderson
"""


with open("aph3IIIa.aln.fas", "r") as f:
    aln = f.read()

aln = aln.split(">")[1:]

cleaned = []

for x in aln:
    x = x.splitlines()[1:]
    x = "".join(x)
    cleaned.append(x)

mins = []

for z in range(len(cleaned)):
    cleaneddddddd = cleaned[:]
    count = 0
    ref = cleaneddddddd[z]
    ref_len = len(ref)
    del cleaneddddddd[z]
    same_same = []
    for x in cleaneddddddd:
        same = []
        for y in range(len(x)):
            if x[y] == ref[y]:
                same.append(x[y])
        same_same.append(len(same))
    mins.append(min(same_same))

min(mins) / len(cleaned[0]) *100

for x in cleaned:
    print(len(x))