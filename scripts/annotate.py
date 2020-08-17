#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functionally annotate a genomic sequence using a Panaroo-output and COBS index. 
"all_annotations.csv" file is created by retrieve_sequences.py prior to Panaroo input. 
"""

import glob 
import subprocess 
import pandas as pd
from io import StringIO
from tqdm import tqdm
import networkx as nx
import cobs_index as cobs
import tempfile
import os 
import re
import time

def predict(temp_dir, title):
    header = os.path.basename(title).split(".fna")[0]
    out =  temp_dir + header + ".gff"
    command = "prodigal -f gff -i " + title + " -o " + out + " -q"
    subprocess.run(command, shell = True)
    return out

def reverse_complement(dna):
    try:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        reversed = ''.join([complement[base] for base in dna[::-1]])
    except:
        reversed = ""
    return reversed

def convert(prod, fasta, count):

    all_gene_ids = []
    source = []
    type = []
    strand = []
    phase = []
    attributes = []
    all_sequences = []
    isolates = []
    regions = []
    starts = []
    ends = []
    
    header = os.path.basename(prod).split(".gff")[0]
    
    with open(prod, "r") as f:
        gff = f.read().splitlines()
    with open(fasta, "r") as fas:
        fasta_file = fas.read()
        
    titles = []
    features = []
    for y in gff:
        if "##sequence-region" in y:
            titles.append(y)
        if not "#" in y:
            features.append(y)

    fas_split = fasta_file.split(">")[1:]
    cleaned_fas_split = []
    for x in fas_split:
        seq = x.splitlines()[1:]
        cleaned_fas_split.append("".join(seq))
        
    for x in features:
        tab_splitted = x.split('\t')
        for s in range(len(fas_split)):
            if tab_splitted[0] in fas_split[s]:
                annotation_seq = cleaned_fas_split[s][int(tab_splitted[3]) - 1: int(tab_splitted[4]) - 3]
                strand.append(tab_splitted[6])
                starts.append(tab_splitted[3])
                ends.append(tab_splitted[4])
                isolates.append(header)
                if tab_splitted[6] == "+":
                  all_sequences.append(annotation_seq)
                elif tab_splitted[6] == "-":
                  all_sequences.append(reverse_complement(annotation_seq))
                regions.append(tab_splitted[0])
                source.append(tab_splitted[1])
                type.append(tab_splitted[2])
                phase.append(tab_splitted[7])
                attributes.append(tab_splitted[8])
                all_gene_ids.append(re.search('ID=(.*?);', tab_splitted[8]).group(1))
        
    titles += features
    titles += ["##FASTA"]
    titles += [fasta_file]
    
    outfile = open(prod, "w")
    outfile.write("\n".join(titles))
    outfile.close()
    
    count += 1
    return count, all_gene_ids, source, type, strand, phase, attributes, all_sequences, isolates, regions, starts, ends
    
def CorA(freq):
    if freq >= 99:
        value = "core"
    else:
        value = "accessory"
    return value

def get_gene_name(ID, IDs, presence_absence):
    index = ([i for i, x in enumerate(IDs) if ID in x])
    if not len(index) == 0:
        index = presence_absence['Gene'][index[0]]
    else:
        index = float("NaN")
    return index

def get_freq(name, gene_table):
    name = name.lower()
    freq = gene_table["Frequency (%)"][list(gene_table["Gene Name"]).index(name)]
    return freq

def get_desc(name, gene_table):
    name = name.lower()
    desc = gene_table["description"][list(gene_table["Gene Name"]).index(name)]
    return desc

def update_file(ref_dir, reference, all_annotations):
    reference_annotations = pd.read_csv(reference + "/all_annotations.csv")

    with open(ref_dir + "/gene_presence_absence.csv", "r") as p:
        correct_headers = p.read()

    split_file = correct_headers.splitlines()
    ID_length = list(range(1, len(split_file[1].split(",")) - 2))
    New_headers = split_file[0] + "," + str(ID_length).replace("[", "").replace("]","").replace(" ", "")
    information = "\n".join([New_headers] + split_file[1:])

    data = StringIO(information)
    presence_absence = pd.read_csv(data, sep=",")

    IDs = presence_absence.iloc[:, 3:]

    IDs = list(pd.Series(IDs.fillna('').values.tolist()).str.join(''))
        
    tqdm.pandas()

    all_annotations = all_annotations.append(reference_annotations, ignore_index=True)
    all_annotations["Name in graph"] = all_annotations.progress_apply(lambda row: get_gene_name(row["ID"], IDs, presence_absence), axis = 1)

    all_annotations.dropna(subset = ["Name in graph"], inplace=True)
    all_annotations = all_annotations.sort_values(by='Name in graph')
    all_annotations = all_annotations.reset_index(drop=True)

    gene_names = []
    frequencies = []
    descriptions = []

    graph = nx.read_gml(ref_dir + "/final_graph.gml")
    for value in graph.graph.values():
        num_isolates = len(value)

    for node in tqdm(graph._node):
        y = graph._node[node]
        gene_names.append(y["name"].lower())
        if not y["description"] == "":
            descriptions.append(y["description"])
        else:
            descriptions.append("hypothetical protein")
        num_sequences = y["seqIDs"]
        unique = set()
        for x in range(len(num_sequences)):
            unique.add(num_sequences[x].split("_")[0])
        frequency = (len(unique)/num_isolates) * 100
        frequencies.append(frequency)

    gene_table = pd.DataFrame(gene_names, columns = ["Gene Name"])
    gene_table["Frequency (%)"] = frequencies
    gene_table["description"] = descriptions
    gene_table = gene_table.drop_duplicates(subset=['Gene Name'])

    gene_table = gene_table.sort_values(by='Gene Name')

    ID_list = []
    number = 0
    for x in all_annotations["Name in graph"].index:
        ID_list.append("PN_" + str(number))
        number += 1
    all_annotations["New ID"] = ID_list
    all_annotations.to_csv(ref_dir + "/fixed_all_annotations.csv", index = False)

    all_annotations["Frequency (%)"] = all_annotations.apply(lambda row: get_freq(row["Name in graph"], gene_table), axis = 1)
    all_annotations["description"] = all_annotations.apply(lambda row: get_desc(row["Name in graph"], gene_table), axis = 1)

    all_annotations["CorA"] = all_annotations["Frequency (%)"].apply(CorA)
    all_annotations.to_csv(ref_dir + "/fixed_all_annotations.csv", index = False)

    return all_annotations

def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary
    """"Translate nucleotide sequences"""
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    protein =""
    try:
        if len(seq)%3 == 0:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3] #Defining a codon as 3 bases
                protein+= table[codon] #Translate the codon into an amino acid based on the dictionary and append this to the protein sequence
    except:
        protein = ""
    return protein

def identify(strand, sequence, index, thresh, isolate):
    #if cora == "accessory":
    if len(sequence) >= 30 and not "_ASM" in isolate:
        try:
            if strand == "+":
                sequence = translate(sequence)
                if not sequence == "":
                    result = index.search(sequence, threshold = thresh)
                else:
                    result = "CHARACTERS"
                if not len(result) == 0:
                    result = (str(result[0]).split("'")[1]).split(";")[0]
                else:
                    result = "[]"
            elif strand == "-":
                sequence = translate(sequence)
                if not sequence == "":
                    result = index.search(sequence, threshold = thresh)
                else:
                    result = "CHARACTERS"
                if not len(result) == 0:
                    result = (str(result[0]).split("'")[1]).split(";")[0]
                else:
                    result = "[]"
            #else:
               # result = ""
            if result == "[]":
                result = "NULL"
        except:
            result = "ERROR"
    else:
        result = "LENGTH/ISOLATE"
    return result

def consensus(panaroo, cobs, CorA):
    if CorA == "acessory":
        if panaroo == cobs:
            value = cobs
        else:
            value = "disagree"
    else:
        value = "core"
    return value

def gff_row(region, source, what_type, start, end, score, strand, phase, ID, description, name, cobs, CorA):
    if CorA == "accessory":
        if name == cobs:
            gff_line = region + '\t' + source + '\t' + what_type + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + str(strand) + '\t' + str(phase) + '\tID=' + ID + ";gbkey=CDS;" + ";gene=" + name + ';gene_biotype=protein_coding' + ";product=" + description + ";locus_tag=" + ID 
        else:
            gff_line = region + '\t' + source + '\t' + what_type + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + str(strand) + '\t' + str(phase) + '\tID=' + ID + ";gbkey=CDS;" + ";gene(Panaroo)=" + name + ";gene(COBS)=" + cobs + ';gene_biotype=protein_coding' + ";product(Panaroo)=" + description + ";locus_tag=" + ID 
    else:
        gff_line = region + '\t' + source + '\t' + what_type + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + str(strand) + '\t' + str(phase) + '\tID=' + ID + ";gbkey=CDS;" + ";gene=" + name + ';gene_biotype=protein_coding' + ";product=" + description + ";locus_tag=" + ID
    return str(gff_line)

def write_gff(gff, source, fasta, isolate):
    
    with open(fasta, "r") as f:
        split = f.read()
        
    lines = split.split(">")[1:]
        
    titles = []
    region_names = []
    for l in lines:
        reg_name = (l.splitlines()[0]).split(" ")[0]
        region_names.append(reg_name)
        titles.append("##sequence-region " + reg_name)
    
    contigs = []
    contigs.append("##isolate " + isolate)
    for title in tqdm(range(len(titles))):        
        sequence_region = region_names[title]
        gff_information = source[source["region"] == sequence_region]
        gff_information = gff_information[["region", "type", "start", "end", "strand", "phase", "attributes", "New ID", "description", "Name in graph", "COBS", "CorA"]]
        gff_information = gff_information.sort_values(by=['start'])
        gff_information["source"] = "Panaroo"
        gff_information["score"] = "."
        print(gff_information)
        try:
            gff_information["line"] = gff_information.apply(lambda row: gff_row(row["region"], row["source"], row["type"], row["start"], row["end"], row["score"], row["strand"], row["phase"], row["New ID"], row["description"], row["Name in graph"], row["COBS"], row["CorA"]), axis = 1)
            gff_information = list(gff_information["line"])
            gff_information = [titles[title]] + gff_information
            gff_information = "\n".join(gff_information)
            contigs.append(gff_information)
        except:
            pass
   
    contigs += ["##FASTA", split]

    with open(gff, 'w') as f:
        f.write("\n".join(contigs))
            
    return 

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
    
    io_opts.add_argument("-g",
                           "--graph_dir",
                           dest="graph_dir",
                           required=True,
                           help='directory of reference Panaroo output',
                           type=str)
                           
    io_opts.add_argument("-f",
                        "--input_FASTA",
                        dest="input_FASTA",
                        required=True,
                        help='the directory of FASTA files to annotate',
                        type=str)
    
    io_opts.add_argument("-a",
                        dest="all_annot",
                        help="file containing all panaroo-updated functional annotations",
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
                        
    args = parser.parse_args()
    return (args)
        
def main():
    args = get_options()
    start = time.time()
    
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    genome_path = glob.glob(args.input_FASTA + "/*.fna")
    
    all_gene_ids = []
    source = []
    type = []
    strand = []
    phase = []
    attributes = []
    all_sequences = []
    isolates = []
    regions = []
    starts = []
    ends = []
    
    print("Predicting CDSs...")
    count = 0
    for genome in tqdm(genome_path):
        out = predict(temp_dir, genome) #use prodigal to predict from fasta files
        gff_file = out
        count,agi,so,ty,st,ph,att,all,h,reg,star,en = convert(gff_file, genome, count) #convert prodigal to prokka
        
        all_gene_ids += agi
        source += so
        type += ty
        strand += st
        phase += ph
        attributes += att
        all_sequences += all
        isolates += h
        regions += reg
        starts += star
        ends += en
        
    print("Writing all annotations file...")
    
    all_files_annotations = pd.DataFrame({"Isolate": isolates, 'ID': all_gene_ids, 'region' : regions, 'source':source,'type':type, 'start': starts, 'end': ends, 'strand' : strand, 'phase': phase, 'attributes': attributes, 'sequence': all_sequences})
    
    all_files_annotations.to_csv(temp_dir + "all_annotations.csv", index = False)

#construct panaroo graph
    print("Generating initial graph...")
    panaroo_command = "panaroo -i " + temp_dir + "*.gff -o " + temp_dir + "unannotated_pangenome --clean-mode moderate"
    subprocess.run(panaroo_command, shell = True)
    
#merge with reference
    print("Merging with reference...")
    merge_command = "panaroo-merge -d " + temp_dir + "unannotated_pangenome " + args.graph_dir +  " -o " + args.output_dir
    subprocess.run(merge_command, shell = True)
    
#update all annotations from graph
    print("Updating all annotations file...")
    all_annotations = update_file(args.output_dir, args.graph_dir, all_files_annotations)

#search index
    print("Searching index...")
    index = cobs.Search(args.input_index)
    tqdm.pandas()
    all_annotations["COBS"] = all_annotations.progress_apply(lambda row: identify(row["strand"], row["sequence"], index, 0.85, row["Isolate"]), axis = 1)
#get consensus
    all_annotations["consensus"] = all_annotations.progress_apply(lambda row: consensus(row["Name in graph"], row["COBS"], row["CorA"]), axis = 1)
    all_annotations.to_csv(args.output_dir + "/annotated_all_annotations.csv", index = False)

#write out gff
    print("Writing annotations...")

    strains = list(set(isolates))
    for isolate in tqdm(range(len(strains))):
        fasta = args.input_FASTA + "/" + strains[isolate] + ".fna"
        gff = args.output_dir + "/" + strains[isolate] + ".gff"
        source = all_annotations[all_annotations["Isolate"] == strains[isolate]]
        write_gff(gff, source, fasta, strains[isolate])
    end = time.time()
    print(str(end-start) + " seconds")
    
    return 

main()
