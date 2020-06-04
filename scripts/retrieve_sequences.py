#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script downloads fasta and gff files in fasta format for accessions of interest. These are then joined and cleaned for direct input into panaroo.
"""
import time
import os

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
     
   # io_opts.add_argument("--complete",
                        #  dest="complete",
                         # help="only retrieve complete genomes",
                         # action='complete',
                         # default=None) #Option to only retrieve complete genomes  
    
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,
                        help="output directory for cleaned gffs",
                        type=str) #Specify the search term for entrez

    io_opts.add_argument("-e",
                        "--email",
                        dest="email",
                        required=True,
                        help="specify email for entrez access",
                        type=str) #Specify the search term for entrez
    
    io_opts.add_argument("--number",
                         dest="number",
                         required=False,
                         help="maximum number of sequences to retrieve (default = 100000)",
                         default=100000) #Option to only retrieve complete genomes  
    
    args = parser.parse_args()
    
    return (args)

def genome_downloader(email, accession, number):
    from Bio import Entrez
    import os
    import urllib.request
    import subprocess

    Entrez.email = email
    
    handle = Entrez.read(Entrez.esearch(db="assembly", term=accession, retmax = number))
    
    assembly_ids = handle['IdList']
                    
    os.chdir('raw_files')

    files = []
    for assembly_id in assembly_ids:
            
        esummary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        esummary_record = Entrez.read(esummary_handle, validate = False)
        
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        files.append(label)
        #get the fasta link - change this to get other formats
        gff_link = os.path.join(url,label+'_genomic.gff.gz') ##https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
        fasta_link = os.path.join(url,label+'_genomic.fna.gz')
        gff_file = '{}.gff.gz'.format(label)
        fasta_file = '{}.fna.gz'.format(label)

        print('Writing:{}'.format(gff_file))
        print('Writing:{}'.format(fasta_file))
        
        urllib.request.urlretrieve(gff_link, gff_file)
        urllib.request.urlretrieve(fasta_link, fasta_file)
        gunzip_gff_command = 'gunzip ' + gff_file
        gunzip_fasta_command = 'gunzip ' + fasta_file
        subprocess.run(gunzip_gff_command, shell = True)
        subprocess.run(gunzip_fasta_command, shell = True)

    os.chdir("..")
    return files

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text


def clean_gffs(headers, input_dir, output_dir):
    import re 
    
    gff = input_dir + '/' + headers + '.gff'
    o = open(gff, 'r')
    stored = str(o.read()) 
    
    d = {"Dbxref=": "dbxref=", "Notes=": "note="}
    stored = replace_all(stored, d) #Ensure GFF annotation format is consistent with prokka output
    
    gene_list = stored.splitlines() 
    title = gene_list[7]
    gene_list = gene_list[9:-1] #remove non-essential content
    
    genes = []
    positions = []
    for y in range(len(gene_list)):
        try:
            if re.split(r'\t+', gene_list[y])[2] == 'gene':
                genes.append(gene_list[y])
        except:
            print(gene_list[y])
    
    cleaned_gffs = "\n".join(str(z) for z in genes)
    cleaned_gffs = title + '\n' + cleaned_gffs

    #Concatenate downloaded fasta and reformatted GFF
    fasta_file = input_dir +'/' + headers + '.fna'
    cleaned_gffs = cleaned_gffs + '\n' + '##FASTA' + '\n' + str((open(fasta_file,'r').read())) 
    filename_cleaned = output_dir + '/' + headers + '_cleaned.gff'
    outfile = open(filename_cleaned,'w')          
    outfile.write(cleaned_gffs)
    outfile.close()

    return

def main():
    start = time.time()
    args = get_options()
    
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    current_dir = os.getcwd() + '/'
    os.makedirs(str(current_dir + 'raw_files'))
    os.makedirs(str(current_dir + args.output_dir))
    accessions = args.accessions.replace('"', "")
    accessions = args.accessions.split(',')
    
    if isinstance(accessions, str):
        headers = genome_downloader(args.email, accessions, int(args.number))
    elif isinstance(accessions, list): 
         headers = []
         for accession in accessions:
             headers.append(str(genome_downloader(args.email, accession, (int(args.number)+1))))#, args.complete)
   
    headers = headers[0].split(',')
    print(headers)
    print(('found {} ids').format(len(headers)))
    
    for header in headers:
        header = header.replace("[","")
        header = header.replace("]","")
        header = header.replace("'","")
        header = header.replace(" ","")
        clean_gffs(header, 'raw_files', args.output_dir)
    
    end = time.time()
    print(end - start)
    
    return

if __name__ == '__main__':
    
    main()
