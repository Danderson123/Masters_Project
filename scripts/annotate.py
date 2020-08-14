#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 15:31:59 2020

@author: danielanderson
"""

import glob 
import subprocess 

def get_options(): #options for downloading and cleaning
    
    import argparse
    description = 'Annotate a genomic sequence in fasta format'
    parser = argparse.ArgumentParser(description=description,
                                     prog='annotate')
    io_opts = parser.add_argument_group('Inputs')
    
    io_opts.add_argument("-i",
                        "--input_index",
                        dest="input_index",
                        required=True,
                        help='the index to annotate from',
                        type=str)
                        
    io_opts.add_argument("-d",
                        "--input_dir",
                        dest="input_dir",
                        required=True,
                        help='the directory of FASTA files to annotate',
                        type=str)
                        
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,  
                        help="output directory for annotated sequences",
                        type=str)
                        
    io_opts.add_argument("--threshold",
                        dest="threshold",
                        help="proportion of matching query kmers (default = 0.8)",
                        default=0.8,
                        type=int)
                        
    io_opts.add_argument("--predict_tool",
                      dest="predict_tool",
                      help="Software to predict open reading frames (prodigal or crude, default = prodigal)",
                      choices=['prodigal', 'crude'],
                      default="prodigal",
                      type=str)
        
    args = parser.parse_args()
    return (args)

def predict(tool, sequence):
    if tool == "crude":
    elif tool == "prokka":
    elif tool == "prodigal":
    elif tool == "glimmer":
    else:
        raise ValueError("Invalid tool specified")     
        
        
        
        
def main():
    args = get_options()
    
    genome_path = glob.glob(args.input_dir + "/*.fna")
    
    for genome in genome path:
        with open(genome, "r") as g:
            seq = g.read()
        
        cds = predict(args.predict_tool, seq)
        
        annotate()