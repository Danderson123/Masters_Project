#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 13:34:41 2020

@author: danielanderson
"""

import glob
import pandas as pd
import matplotlib.pyplot as plt 

one_kmer_length = glob.glob("results/*100_1.csv")
two_kmer_length = glob.glob("results/*100_2.csv")
four_kmer_length = glob.glob("results/*100_4.csv")
eight_kmer_length = glob.glob("results/*100_8.csv")
sixteen_kmer_length = glob.glob("results/*100_16.csv")

def calculate_sens(TP,FN):
    return TP/(TP+FN)

def calculate_fpr(FP,TN):
    FPR = 0
    if TN == -1:
        TN = 0
    FPR = FP/(TN+FP)
    return FPR

total_thresholds = []
total_TP = []
total_FP = []
total_FN = []
total_TN = []

for file in one_kmer_length:
    result = pd.read_csv(file)
    total_thresholds += list(result["threshold"])
    total_TP += list(result["TP"])
    total_FP += list(result["FP"])
    total_FN += list(result["FN"])
    total_TN += list(result["TN"])
    
total_one = pd.DataFrame(total_thresholds, columns = ["thresholds"])
total_one["TP"] = total_TP
total_one["FP"] = total_FP
total_one["FN"] = total_FN
total_one["TN"] = total_TN
total_one["TPR"] = total_one.apply(lambda row: calculate_sens(row["TP"], row["FN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_one["FPR"] = total_one.apply(lambda row: calculate_fpr(row["FP"], row["TN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_one["FPR"] = total_one["FPR"].fillna(0)
total_one["kmer_length"] = 1


total_thresholds = []
total_TP = []
total_FP = []
total_FN = []
total_TN = []

for file in two_kmer_length:
    result = pd.read_csv(file)
    total_thresholds += list(result["threshold"])
    total_TP += list(result["TP"])
    total_FP += list(result["FP"])
    total_FN += list(result["FN"])
    total_TN += list(result["TN"])

total_two = pd.DataFrame(total_thresholds, columns = ["thresholds"])
total_two["TP"] = total_TP
total_two["FP"] = total_FP
total_two["FN"] = total_FN
total_two["TN"] = total_TN
total_two["TPR"] = total_two.apply(lambda row: calculate_sens(row["TP"], row["FN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_two["FPR"] = total_two.apply(lambda row: calculate_fpr(row["FP"], row["TN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_two["FPR"] = total_two["FPR"].fillna(0)

total_two["kmer_length"] = 2



total_thresholds = []
total_TP = []
total_FP = []
total_FN = []
total_TN = []

for file in four_kmer_length:
    result = pd.read_csv(file)
    total_thresholds += list(result["threshold"])
    total_TP += list(result["TP"])
    total_FP += list(result["FP"])
    total_FN += list(result["FN"])
    total_TN += list(result["TN"])

total_four = pd.DataFrame(total_thresholds, columns = ["thresholds"])
total_four["TP"] = total_TP
total_four["FP"] = total_FP
total_four["FN"] = total_FN
total_four["TN"] = total_TN
total_four["TPR"] = total_four.apply(lambda row: calculate_sens(row["TP"], row["FN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_four["FPR"] = total_four.apply(lambda row: calculate_fpr(row["FP"], row["TN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_four["FPR"] = total_four["FPR"].fillna(0)


total_four["kmer_length"] = 4



total_thresholds = []
total_TP = []
total_FP = []
total_FN = []
total_TN = []

for file in eight_kmer_length:
    result = pd.read_csv(file)
    total_thresholds += list(result["threshold"])
    total_TP += list(result["TP"])
    total_FP += list(result["FP"])
    total_FN += list(result["FN"])
    total_TN += list(result["TN"])

total_eight = pd.DataFrame(total_thresholds, columns = ["thresholds"])
total_eight["TP"] = total_TP
total_eight["FP"] = total_FP
total_eight["FN"] = total_FN
total_eight["TN"] = total_TN
total_eight["TPR"] = total_eight.apply(lambda row: calculate_sens(row["TP"], row["FN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_eight["FPR"] = total_eight.apply(lambda row: calculate_fpr(row["FP"], row["TN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_eight["FPR"] = total_eight["FPR"].fillna(0)

total_eight["kmer_length"] = 8



total_thresholds = []
total_TP = []
total_FP = []
total_FN = []
total_TN = []

for file in sixteen_kmer_length:
    result = pd.read_csv(file)
    total_thresholds += list(result["threshold"])
    total_TP += list(result["TP"])
    total_FP += list(result["FP"])
    total_FN += list(result["FN"])
    total_TN += list(result["TN"])

total_sixteen = pd.DataFrame(total_thresholds, columns = ["thresholds"])
total_sixteen["TP"] = total_TP
total_sixteen["FP"] = total_FP
total_sixteen["FN"] = total_FN
total_sixteen["TN"] = total_TN
total_sixteen["TPR"] = total_sixteen.apply(lambda row: calculate_sens(row["TP"], row["FN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_sixteen["FPR"] = total_sixteen.apply(lambda row: calculate_fpr(row["FP"], row["TN"]), axis=1) #Output notes of key 5' guide RNA characterstics to a new column
total_sixteen["FPR"] = total_sixteen["FPR"].fillna(0)

total_sixteen["kmer_length"] = 16


# =============================================================================
# total_one.to_csv("total_one.csv")
# total_two.to_csv("total_two.csv")
# total_four.to_csv("total_four.csv")
# total_eight.to_csv("total_eight.csv")
# total_sixteen.to_csv("total_sixteen.csv")
# =============================================================================


# =============================================================================
# import seaborn as sns
# 
# sns.lineplot(x="FPR", y="TPR", hue="thresholds", markers=True, dashes=False, data=total_one)
# sns.lineplot(x="FPR", y="TPR", hue="thresholds", markers=True, dashes=False, data=total_two)
# sns.lineplot(x="FPR", y="TPR", hue="thresholds", markers=True, dashes=False, data=total_four)
# sns.lineplot(x="FPR", y="TPR", hue="thresholds", markers=True, dashes=False, data=total_eight)
# sns.lineplot(x="FPR", y="TPR", hue="thresholds", markers=True, dashes=False, data=total_sixteen)
# 
# =============================================================================

thou_point_two = total_sixteen[total_sixteen["thresholds"] == 0.2]
thou_point_four = total_sixteen[total_sixteen["thresholds"] == 0.4]
thou_point_six = total_sixteen[total_sixteen["thresholds"] == 0.6]
thou_point_eight = total_sixteen[total_sixteen["thresholds"] == 0.8]
thou_one = total_sixteen[total_sixteen["thresholds"] == 1.0]

plt.plot(thou_point_two["FPR"], thou_point_two["TPR"])
plt.plot(thou_point_four["FPR"], thou_point_four["TPR"])
plt.plot(thou_point_six["FPR"], thou_point_six["TPR"])
plt.plot(thou_point_eight["FPR"], thou_point_eight["TPR"])
plt.plot(thou_one["FPR"], thou_one["TPR"])


groups = total_one.groupby('thresholds')
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.FPR, group.TPR, marker='o', linestyle='', ms=3, label=name)
ax.legend()
plt.show()

groups = total_two.groupby('thresholds')
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.FPR, group.TPR, marker='o', linestyle='', ms=3, label=name)
ax.legend()
plt.show()

groups = total_four.groupby('thresholds')
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.FPR, group.TPR, marker='o', linestyle='', ms=3, label=name)
ax.legend()
plt.show()

groups = total_eight.groupby('thresholds')
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.FPR, group.TPR, marker='o', linestyle='', ms=3, label=name)
ax.legend()
plt.show()

groups = total_sixteen.groupby('thresholds')
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.FPR, group.TPR, marker='o', linestyle='', ms=3, label=name)
ax.legend()
plt.show()
