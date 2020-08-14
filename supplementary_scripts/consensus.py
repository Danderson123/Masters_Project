# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
from tqdm import tqdm 

all_annotations = pd.read_csv("INDEX_SEARCHED_all_annotations.csv")

list(all_annotations["COBS 85%"]).count("pspC~~~~~~cbpA")

tqdm.pandas()

def consensus(panaroo, cobs):
    if not "group_" in panaroo:
        if panaroo == cobs:
            value = cobs
        else:
            value = "disagree"
    else:
        value = "grp"
    return value 

all_annotations["consensus 75%"] = all_annotations.progress_apply(lambda row: consensus(row["Name in graph"], row["COBS 75%"]), axis = 1)
all_annotations["consensus 80%"] = all_annotations.progress_apply(lambda row: consensus(row["Name in graph"], row["COBS 80%"]), axis = 1)
all_annotations["consensus 85%"] = all_annotations.progress_apply(lambda row: consensus(row["Name in graph"], row["COBS 85%"]), axis = 1)
all_annotations["consensus 90%"] = all_annotations.progress_apply(lambda row: consensus(row["Name in graph"], row["COBS 90%"]), axis = 1)
all_annotations["consensus 95%"] = all_annotations.progress_apply(lambda row: consensus(row["Name in graph"], row["COBS 95%"]), axis = 1)

list(all_annotations["consensus 80%"]).count("ermB")

264/616*100

def COBS_freq(consensus, total):
    if not consensus == "grp" or not consensus == "disagree":
        num = total.count(consensus)
        freq = num / 616 * 100
    else:
        freq = ""
    return freq 

all_annotations["COBS 75% Frequency"] = all_annotations.progress_apply(lambda row: COBS_freq(row["consensus 75%"], list(all_annotations["consensus 75%"])), axis = 1)
all_annotations["COBS 80% Frequency"] = all_annotations.progress_apply(lambda row: COBS_freq(row["consensus 80%"], list(all_annotations["consensus 80%"])), axis = 1)
all_annotations["COBS 85% Frequency"] = all_annotations.progress_apply(lambda row: COBS_freq(row["consensus 85%"], list(all_annotations["consensus 85%"])), axis = 1)
all_annotations["COBS 90% Frequency"] = all_annotations.progress_apply(lambda row: COBS_freq(row["consensus 90%"], list(all_annotations["consensus 90%"])), axis = 1)
all_annotations["COBS 75% Frequency"] = all_annotations.progress_apply(lambda row: COBS_freq(row["consensus 95%"], list(all_annotations["consensus 95%"])), axis = 1)

list(all_annotations["consensus"]).count("disagree")


con = []
for x in tqdm(range(len(list(all_annotations["CorA"])))):
    if all_annotations["CorA"][x] == "core" and all_annotations["consensus"][x] == "disagree":
        con.append(x)

list(all_annotations["COBS 85%"]).count("zmpC")/ 616 * 100
list(all_annotations["COBS 85%"]).count("pspC~~~~~~cbpA") / 616 * 100

for x in range(len(list(all_annotations["Name in graph"]))):
        if "erm" in all_annotations["Name in graph"][x]:
            print(x)
            
print(all_annotations["COBS 85%"][52439])

list(all_annotations["consensus 95%"]).count("zmpC")/ 616 * 100

aph3IIIa=0.16,0.16,0.16,0.16,0.16
mefA=13.15,13.15,13.15,13.15,13.15  
tetM=9.74,9.74,9.74,9.74,9.09
pspA=22.73,22.56,19.97,16.40,14.12
pspC=49.19,48.86,42.86,34.09,30.03
zmpC=14.61,14.61,14.61,14.61,14.61