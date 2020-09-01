#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script downloads fasta and gff files for accessions of interest. 
These are then joined and cleaned for direct input into panaroo.
"""
import time
import os
from Bio import Entrez
import urllib.request
import re
import sys
import shutil
import tempfile
from tqdm import tqdm
import gzip
import pandas as pd

def translate(seq): #Translate the exon sequence of the gene into its respective amino acid codes using a dictionary
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
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

def reverse_complement(dna):
    try:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        reversed = ''.join([complement[base] for base in dna[::-1]])
    except:
        reversed = ""
    return reversed


def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

def get_options(): #options for downloading and cleaning
    
    import argparse

    description = 'Download and process assemblies for panaroo input'
    parser = argparse.ArgumentParser(description=description,
                                     prog='retrieve')

    io_opts = parser.add_argument_group('Entrez')
    
    io_opts.add_argument("-s",
                        "--search_term",
                        dest="accessions",
                        required=True,
                        help='sequences to download, specify "genues species" for all accessions or specific accessions of interest. Separate accessions by "," with no spaces',
                        type=str) #Specify the search term for entrez

    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="output directory for cleaned gffs",
                        type=str)
                        
    io_opts.add_argument("-e",
                        "--email",
                        dest="email",
                        required=True,
                        help="specify email for entrez access",
                        type=str)
    
    io_opts.add_argument("--number",
                         dest="number",
                         required=False,
                         help="maximum number of sequences to retrieve",
                         default=100)
     
    io_opts.add_argument("--dirty",
                      dest="dirty",
                      help="keep temporary directory containing raw genome files and split contigs",
                      action='store_true',
                      default=False)
        
    args = parser.parse_args()
    
    return (args)

def genome_downloader(email, accession, number, output_dir):
    """Searches for available GFF files and retrieves fastas only for these"""

    Entrez.email = email

    handle = Entrez.read(Entrez.esearch(db="assembly", term=accession, retmax = int(number)))

    assembly_ids = handle['IdList']
                    
    os.chdir(output_dir + 'raw_files')

    files = []
    
    outfile = open("Accessions.txt",'w')
    outfile.write('\n'.join(assembly_ids))
    outfile.close()
    
    for assembly_id in assembly_ids:
        esummary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)

        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        
        if url == '':
            url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
            
        label = os.path.basename(url)
        files.append(label)
        #get the fasta link - change this to get other formats
        gff_link = os.path.join(url,label+'_genomic.gff.gz') ##https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
        fasta_link = os.path.join(url,label+'_genomic.fna.gz')
        gff_file = '{}.gff.gz'.format(label)
        fasta_file = '{}.fna.gz'.format(label)
        
        urllib.request.urlretrieve(gff_link, gff_file)
        urllib.request.urlretrieve(fasta_link, fasta_file)
        
    os.chdir('../../..')
    return files
    
def split_contigs(headers, input_dir, output_dir):
    
    """Most annotations are incomplete and consist of multiple contigs. Panaroo only accepts prokka-formatted genomic regions. Similarly, CDSs not starting in ATG are not accepted."""
    all_gene_ids = []
    source = []
    type = []
    phase = []
    attributes = []
    all_sequences = []
    strand = []
    isolates = []
    regions = []
    starts = []
    ends = []
    
    all_raw = []
    non_cds = []
    non_divis = []
    duplicates = []
    not_start = []
    stop_in =[]
    
    for header in tqdm(headers):
        d = {"[": "", "]": "", " ": "", "'": ""}
        header = replace_all(header, d)
        fasta = input_dir + header + '.fna.gz'
        gff = input_dir + header + '.gff.gz'
        
        fasta = input_dir + header + '.fna.gz'
        with gzip.open(fasta, 'rt') as f:
            stored_fasta = str(f.read())
       
            #fasta = input_dir + header + '.fna'
            #with open(fasta, 'rt') as f:
               # stored_fasta = str(f.read())
                
        with gzip.open(gff, 'rt') as g:
            stored_gff = str(g.read())
        
        stored_fasta = stored_fasta.split('>')
        stored_fasta = stored_fasta[1:]
        
        stored_gff = stored_gff.split("##sequence-region ")
        stored_gff = stored_gff[1:]
        
        all_region_names = []
        all_region_annotations = []
        all_region_sequences = []
        
        for gff_region in range(len(stored_gff)):

            d = {"Dbxref=": "dbxref=", "Notes=": "note="}
            stored_gff[gff_region] = replace_all(stored_gff[gff_region], d) #Ensure GFF annotation format is consistent with prokka output
              
            gene_list = stored_gff[gff_region].splitlines()
            title = gene_list[0].split("\n")[0]
            for line in range(len(gene_list)):
              if '###' in gene_list[line]:
                  gene_list = gene_list[:line]

            gene_list = gene_list[2:]

            genes = []
            for y in range(len(gene_list)):
                all_raw.append(gene_list[y])
                split = re.split(r'\t+', gene_list[y])
                end = int(split[4])
                start = int(split[3])
                if split[2] == 'CDS' and (((end - start) + 1) % 3) == 0: #process_pokka input only looks for CDSs and returns duplicate error when the exon is split.
                    genes.append(gene_list[y])
                elif not split[2] == 'CDS':
                    non_cds.append(gene_list[y])
                elif split[2] == 'CDS' and not (((end - start) + 1) % 3) == 0:
                    non_divis.append(gene_list[y])
            
            Ids= []
            for gene in range(len(genes)):
                Id = re.search('ID=(.*?);', genes[gene]).group(1)
                Ids.append(Id)

            output_cds = []
            start_index = []
            end_index = []
            sense = []
            for Id in range(len(Ids)):
                if Ids.count(Ids[Id]) == 1:
                  output_cds.append(genes[Id])
                  no_duplicates_split = re.split(r'\t+', genes[Id])
                  start_index.append(int(no_duplicates_split[3]) - 1)
                  end_index.append(int(no_duplicates_split[4]) - 1)
                  sense.append(no_duplicates_split[6])
                else:
                    duplicates.append(Ids[Id])
            
            nucleotides = ''.join(stored_fasta[gff_region].split('\n')[1:])

            to_remove = []
            for sign in range(len(sense)):
                if sense[sign] == '+':
                    start_codon = nucleotides[start_index[sign]: start_index[sign] + 3]
                    protein = translate(nucleotides[start_index[sign]: end_index[sign] - 2])
                    if not (start_codon == "ATG" or start_codon == "TTG" or start_codon == "GTG") or "*" in protein or protein == "":
                        to_remove.append(sign)
                        if not (start_codon == "ATG" or start_codon == "TTG" or start_codon == "GTG"):
                            not_start.append(start_codon)
                        if "*" in protein:
                            stop_in.append(protein)
                else:
                    reverse_start_codon = nucleotides[end_index[sign] - 2: end_index[sign] + 1]
                    protein_reverse = translate(reverse_complement(nucleotides[start_index[sign] + 3: end_index[sign] + 1]))
                    if not (reverse_start_codon == "CAT" or reverse_start_codon == "CAA" or reverse_start_codon == "CAC") or "*" in protein_reverse or protein_reverse == "":
                        to_remove.append(sign)
                        if not (reverse_start_codon == "CAT" or reverse_start_codon == "CAA" or reverse_start_codon == "CAC"):
                            not_start.append(start_codon)
                        if "*" in protein_reverse:
                            stop_in.append(protein_reverse)
            
            for index in sorted(to_remove, reverse=True):
                del output_cds[index]
                    
            if len(output_cds) == 0:
              continue
    
            all_region_names.append("##sequence-region " + title)
                        
            cleaned_gffs = "\n".join(str(z) for z in output_cds)
            
            for annotation_line in range(len(output_cds)):
                tab_splitted = output_cds[annotation_line].split('\t')
                regions.append(tab_splitted[0])
                source.append(tab_splitted[1])
                type.append(tab_splitted[2])
                phase.append(tab_splitted[7])
                attributes.append(tab_splitted[8])
                all_gene_ids.append(re.search('ID=(.*?);', tab_splitted[8]).group(1))
                strand.append(tab_splitted[6])
                isolates.append(header)
                annotation_seq = nucleotides[int(tab_splitted[3]) - 1: int(tab_splitted[4]) - 3]
                starts.append(tab_splitted[3])
                ends.append(tab_splitted[4])
                if tab_splitted[6] == "+":
                    all_sequences.append(annotation_seq)
                elif tab_splitted[6] == "-":
                    all_sequences.append(reverse_complement(annotation_seq))
                
            all_region_annotations.append(cleaned_gffs)
            #Concatenate downloaded fasta and reformatted GFF
            fasta_file = ">" + stored_fasta[gff_region]
            all_region_sequences.append(fasta_file)
            
        annotated_file = "\n".join(all_region_names + all_region_annotations) + "\n##FASTA\n" + "".join(all_region_sequences)
        filename_cleaned = output_dir + header + '.gff'
        outfile = open(filename_cleaned,'w')
        outfile.write(annotated_file)
        outfile.close()

    all_files_annotations = pd.DataFrame({"Isolate": isolates, 'ID': all_gene_ids, 'region' : regions, 'source':source,'type':type, 'start': starts, 'end': ends, 'strand' : strand, 'phase': phase, 'attributes': attributes, 'sequence': all_sequences})
    
    all_files_annotations.to_csv(output_dir + "all_annotations.csv", index=False)
    
    print("raw annotations:" + str(len(all_raw)))
    print("Non-CDS:" + str(len(non_cds)))
    print("Non-divisible:" + str(len(non_divis)))
    print("duplicates:" + str(len(duplicates)))
    print("not_start:" + str(len(not_start)))
    print("stop_in:" + str(len(stop_in)))
    print("all_valid:" + str(len(all_gene_ids)))
    
    return

def main():
    start = time.time()
    args = get_options()
    
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
        
    args.output_dir = os.path.join(args.output_dir, "")
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    
    os.makedirs(str(temp_dir + 'raw_files'))
    
    accessions = args.accessions.replace('"', "")
    if "," in accessions:
        accessions = args.accessions.split(',')
    
    failed = []
    
    print("Retrieving Sequences")

    if isinstance(accessions, str):
        headers = genome_downloader(args.email,
                                    [accessions],
                                    int(args.number),
                                    temp_dir)
    elif isinstance(accessions, list):
        headers = []
        for accession in tqdm(accessions):
            try:
                headers.append(str(genome_downloader(args.email,
                                                     accession,
                                                     int(args.number),
                                                     temp_dir)))
            except:
                failed.append(accession)
    else:
        raise RuntimeError("Accessions are invalid")
        
    
    failed_out = open(temp_dir + "failed_accessions.txt", "w")
    failed_out.write(",".join(failed))
    failed_out.close()
    
    print(('found {} ids').format(len(headers)))
    print("Processing GFF files")
    
    split_contigs(headers,
                  temp_dir + 'raw_files/',
                  args.output_dir)
   
    end = time.time()
    print(end - start)
    print("Done")
    
    #remove temporary directory if dirty = False
    if not args.dirty:
        shutil.rmtree(temp_dir)
        
    sys.exit(0)

if __name__ == '__main__':
    main()
