#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 12:15:10 2020

@author: danielanderson
"""

from random import randint
import cobs_index as cobs


#def optimal_kmer_length(protein_list, kmers, num_samples):

    #one = []
    #two = []
    #four = []
    #eight = []
    #sixteen = []
    
    #for y in subset:
        #sub_one = []
        #sub_two = []
        #sub_four = []
        #sub_eight = []
        #sub_sixteen = []        
        #for z in range(len(y)):
            #sub_one.append(y[z])
            #if len(y[z:z+2]) == 2:
                #sub_two.append(y[z:z+2])
            #if len(y[z:z+4]) == 4:
                #sub_four.append(y[z:z+4])
            #if len(y[z:z+8]) == 8:
                #sub_eight.append(y[z:z+8])
            #if len(y[z:z+16]) == 16:
                #sub_sixteen.append(y[z:z+16])
        #one.append(sub_one)
        #two.append(sub_two)
        #four.append(sub_four)
        #eight.append(sub_eight)
        #sixteen.append(sub_sixteen)
    
    #dictionary = {"one":one,"two":two,"four":four,"eight":eight,"sixteen":sixteen}
    #return dictionary

def index_creation(test_kmers):
    
   # 1_param = cobs.CompactIndexParameters()
    #1_param.term_size = 1               # k-mer size
    #1_param.clobber = True               # overwrite output and temporary files
    #1_param.false_positive_rate = 0.4    # higher false positive rate -> smaller index

    for num in test_kmers:
        ten_param = cobs.CompactIndexParameters()
        ten_param.term_size = num               # k-mer size
        ten_param.clobber = True               # overwrite output and temporary files
        ten_param.false_positive_rate = 0.4    # higher false positive rate -> smaller index
        cobs.compact_construct("/ten", str(num) +"_" + "ten_index.cobs_compact", index_params=ten_param)

    #10_param = cobs.CompactIndexParameters()
    #10_param.term_size = 1               # k-mer size
    #10_param.clobber = True               # overwrite output and temporary files
   # 10_param.false_positive_rate = 0.4    # higher false positive rate -> smaller index

    #1_param = cobs.CompactIndexParameters()
    #1_param.term_size = 1               # k-mer size
    #1_param.clobber = True               # overwrite output and temporary files
    #1_param.false_positive_rate = 0.4    # higher false positive rate -> smaller index

 #   cobs.compact_construct("/one", "one_index.cobs_compact")
  #  cobs.compact_construct("/hundred", "hundred_index.cobs_compact")
   # cobs.compact_construct("/thousand", "thousand_index.cobs_compact")
    return

def main():
    
    with open("combined_protein_cdhit_out.txt", "r") as f:
        protein_file = f.read().split(">")[1:]
    
    titles = []
    cleaned = []
    for x in protein_file:
        cleaned.append("".join(x.splitlines()[1:]))
        titles.append(x.splitlines()[0])
    
# =============================================================================
#     #Generate a random list of indices of range num_samples
#     one_indexes = [randint(0,(len(cleaned) -1)) for i in range(1)]
#     ten_indexes = [randint(0,(len(cleaned) -1)) for i in range(10)]
#     onehundred_indexes = [randint(0,(len(cleaned) -1)) for i in range(100)]
#     onethousand_indexes = [randint(0,(len(cleaned) -1)) for i in range(1000)]
#     #subset = [cleaned[i] for i in indexes]
# 
#     for a in one_indexes:
#         info = titles[a] + "\n" + cleaned[a]
#         one = open("one/" + titles[a] + ".txt", "w")
#         one.write(info)
#         one.close()
#         
#     for b in ten_indexes:
#         info = titles[b] + "\n" + cleaned[b]
#         ten = open("ten/" + titles[b] + ".txt", "w")
#         ten.write(info)
#         ten.close()
#         
#     for c in onehundred_indexes:
#         info = titles[c] + "\n" + cleaned[c]
#         hun = open("onehundred/" + titles[c] + ".txt", "w")
#         hun.write(info)
#         hun.close()
# 
#     for d in onethousand_indexes:
#         info = titles[d] + "\n" + cleaned[d]
#         thou = open("onethousand/" + titles[d] + ".txt", "w")
#         thou.write(info)
#         thou.close()
# 
#     for y in range(len(titles)):
#         information = titles[y] + "\n" + cleaned[y]
#         file = open("DocumentList/" + titles[y] + ".txt", "w")
#         file.write(information)
#         file.close()
# =============================================================================
        
    test_kmers = [1,2,4,8,16]
    index_creation(test_kmers)
    #ls = optimal_kmer_length(cleaned, test_kmers, 10)
    
    