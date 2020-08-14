# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 18:50:47 2020

@author: danielanderson
"""
import glob
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

tenthousand = glob.glob("cluster_results/core/*.csv")

files_dict = []

kmer_dict = {}
for kmer in tqdm([1,2,4,6,8,10]):
    sens_fpr_dict = {}
    for threshold in tqdm(["0","05","1","15","2","25", "3","35","4","45", "5","55","6","65", "7", "75", "8","85", "9","95", "10"]):
        Result = []
        print("importing files...")
        for x in tenthousand:
            #try:
            df = pd.read_csv(x, engine= "python")
            if "group_" in df["Gene Name"][0]:
                continue
            #except:
                #continue
            if kmer == 1:
                Result.append(df["Result_"+threshold][0])
            elif kmer == 2:
                Result.append(df["Result_"+threshold][1])
            elif kmer == 4:
                Result.append(df["Result_"+threshold][2])
            elif kmer == 6:
                Result.append(df["Result_"+threshold][3])
            elif kmer == 8:
                Result.append(df["Result_"+threshold][4])
            elif kmer == 10:
                Result.append(df["Result_"+threshold][5])
        
        print("calculating sensitivity and specificity...")
        
        try:
            kmer_sens = (Result.count("TP")) / (Result.count("TP") + Result.count("FN"))
        except:
            kmer_sens = 0
            
            
        kmer_fpr = (Result.count("FP")) / (Result.count("FP") + Result.count("TN"))
   
        sens_fpr_local = {"sensitivity" : kmer_sens, "false positive rate" : kmer_fpr}
        sens_fpr_dict.update({threshold: sens_fpr_local})           

    kmer_dict.update({kmer : sens_fpr_dict})
files_dict.append(kmer_dict)
    
print("plotting...")
to_plot_files = files_dict[0]
one_kmers = to_plot_files[1]
two_kmers =to_plot_files[2]
four_kmers =to_plot_files[4]
six_kmers =to_plot_files[6]
eight_kmers =to_plot_files[8]
ten_kmers =to_plot_files[10]
fig, ax = plt.subplots()
fig_2, ax_2 = plt.subplots()
kmer_list_dict = {1:one_kmers,2:two_kmers,4:four_kmers, 6:six_kmers,8:eight_kmers,10:ten_kmers}#,20:twenty_kmers}
#kmer_list_dict = {10:ten_kmers}#,20:twenty_kmers}
for kmer_length, kmer in kmer_list_dict.items():
    sense = []
    fpr = []
    thresh_list = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]
    for key, value in kmer.items():
        for k, v in value.items():
            if k == "sensitivity":
                sense.append(v)
            else:
                fpr.append(v)
                
# =============================================================================
#     sense_str = [str(i) for i in sense]
#     fpr_str = [str(i) for i in fpr]
#+  
#     out_sense = open(str(kmer_length)+"_sense.txt", "w")
#     out_sense.write("\n".join(sense_str))
#     out_sense.close()
#      
#     out_fpr = open(str(kmer_length)+"_fpr.txt", "w")
#     out_fpr.write("\n".join(fpr_str))
#     out_fpr.close()
# =============================================================================
      

    ax.plot(fpr, sense, label = str(kmer_length))
    ax.set_xlim([0.2,0.5])
    ax.set_ylim([0.4,1])
    ax.scatter(fpr, sense, s=15)
# =============================================================================
#     thresh_labels = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
#     if kmer_length == 10:
#         for i, txt in enumerate(thresh_labels):
#             ax.annotate(txt, (fpr[i], sense[i]))
#     else:
#         pass
     
    #ax_2.plot(thresh_list, sense, label = str(kmer_length))
   # ax_2.set_xlim([0,1])
   # ax_2.set_ylim([0,1])
   # ax_2.scatter(thresh_list, sense, s=15)
    
#ax.legend(title="k-mer length", loc = 'lower right')
ax.set_xlabel("false positive rate")
ax.set_ylabel("sensitivity")
import numpy as np
x = np.linspace(-5,5,100)
ax.plot(x, x, '-k', linestyle='--', linewidth = 0.7)
# =============================================================================
# ax_2.legend(title="k-mer length")
# ax_2.set_xlabel("threshold")
# ax_2.set_ylabel("sensitivity")
# =============================================================================

#fig.savefig('centroid_core_zoomed.pdf')

#fig_2.savefig('accessory_sensitivity-threshold.pdf')

