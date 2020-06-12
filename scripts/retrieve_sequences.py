#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script downloads fasta and gff files in fasta format for accessions of interest. These are then joined and cleaned for direct input into panaroo.
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

from .integrate import replace_all

def get_options(): #options for downloading and cleaning
    
    import argparse

    description = 'Download and process assemblies for panaroo input'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo-retrieve')

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
    
    for assembly_id in tqdm(assembly_ids):
            
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
    
    """Most annotations are incomplete and consist of multiple contigs. Panaroo only accepts single genomic regions. Similarly, CDSs not starting in ATG are not accepted."""
    for header in tqdm(headers):
        d = {"[": "", "]": "", " ": "", "'": ""}
        header = replace_all(header, d)
        fasta = input_dir + header + '.fna.gz'
        gff = input_dir + header + '.gff.gz'
        
        with gzip.open(fasta, 'rt') as f:
            stored_fasta = str(f.read())
        
        with gzip.open(gff, 'rt') as g:
            stored_gff = str(g.read())

        stored_fasta = stored_fasta.split('>')
        stored_fasta = stored_fasta[1:]
        
        stored_gff = stored_gff.split("##sequence-region ")
        stored_gff = stored_gff[1:]
        
        for gff_region in range(len(stored_gff)):

            d = {"Dbxref=": "dbxref=", "Notes=": "note="}
            stored_gff[gff_region] = replace_all(stored_gff[gff_region], d) #Ensure GFF annotation format is consistent with prokka output
              
            gene_list = stored_gff[gff_region].splitlines()

            title = gene_list[0].split(" ")[0]
            for line in range(len(gene_list)):
              if '###' in gene_list[line]:
                  gene_list = gene_list[:line]

            gene_list = gene_list[2:]

            genes = []
            for y in range(len(gene_list)):
              split = re.split(r'\t+', gene_list[y])
              end = int(split[4])
              start = int(split[3])
              if split[2] == 'CDS' and (((end - start) + 1) % 3) == 0: #process_pokka input only looks for CDSs and returns duplicate error when the exon is split.
                  genes.append(gene_list[y])
              else:
                  pass

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
            
            nucleotides = ''.join(stored_fasta[gff_region].split('\n')[1:])
            
            to_remove = []
            for sign in range(len(sense)):
                if sense[sign] == '+':
                    if not nucleotides[start_index[sign]: start_index[sign] + 3] == "ATG":
                        to_remove.append(sign)
                else:                     
                    if not nucleotides[end_index[sign] - 2: end_index[sign] + 1] == "CAT":
                        to_remove.append(sign)
            
            for index in sorted(to_remove, reverse=True):
                del output_cds[index]
                    
            if len(output_cds) == 0:
              continue

            cleaned_gffs = "\n".join(str(z) for z in output_cds)
            cleaned_gffs = title + '\n' + cleaned_gffs
            #Concatenate downloaded fasta and reformatted GFF
            fasta_file = stored_fasta[gff_region]
            cleaned_gffs = cleaned_gffs + '\n' + '##FASTA' + '\n' + ">" + fasta_file
            filename_cleaned = output_dir + title + '.gff'
            outfile = open(filename_cleaned,'w')
            outfile.write("##sequence-region " + cleaned_gffs)
            outfile.close()
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
    print("Retrieving Sequences")
    
    if isinstance(accessions, str):
        headers = genome_downloader(args.email,
                                    [accessions],
                                    int(args.number),
                                    temp_dir)
    elif isinstance(accessions, list):
        headers = []
        for accession in accessions:
            headers.append(str(genome_downloader(args.email,
                                                 accession,
                                                 int(args.number),
                                                 temp_dir)))
    else:
        raise RuntimeError("Accessions are invalid")
        
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
