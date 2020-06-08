#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script downloads fasta and gff files in fasta format for accessions of interest. These are then joined and cleaned for direct input into panaroo.
"""
import time
import os

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
                         help="maximum number of sequences to retrieve",
                         default=100)
    
    args = parser.parse_args()
    
    return (args)

def genome_downloader(email, accession, number):
    
   """Searches for available GFF files and retrieves fastas only for these"""

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

        print('Writing:{}'.format(label))
        
        urllib.request.urlretrieve(gff_link, gff_file)
        urllib.request.urlretrieve(fasta_link, fasta_file)
        gunzip_gff_command = 'gunzip ' + gff_file
        gunzip_fasta_command = 'gunzip ' + fasta_file
        subprocess.run(gunzip_gff_command, shell = True)
        subprocess.run(gunzip_fasta_command, shell = True)

    os.chdir("..")
    return files


def split_contigs(headers, input_dir, output_dir):
    
    """Most annotations are incomplete and consist of multiple contigs. Panaroo only accepts single genomic regions."""
    
    for header in headers:
        d = {"[": "", "]": "", " ": "", "'": ""}
        header = replace_all(header, d)
        fasta = input_dir + '/' + header + '.fna'
        gff = input_dir + '/' + header + '.gff'
        
        with open(fasta, 'r') as f:
            stored_fasta = str(f.read()) 
            
        with open(gff, 'r') as g:
            stored_gff = str(g.read()) 
            
        stored_fasta = stored_fasta.split('>')
        stored_fasta = stored_fasta[1:]
        
        for fasta_region in stored_fasta:
            filename_split = 'contigs' + '/' + fasta_region.split(" ")[0] + ".fna"
            outfile = open(filename_split,'w')          
            outfile.write(">" + fasta_region)
            outfile.close()
        
        stored_gff = stored_gff.split("##sequence-region ")
        stored_gff = stored_gff[1:]
        
        for gff_region in stored_gff:
               filename_split = 'contigs' + '/' + gff_region.split(" ")[0] + ".gff"
               outfile = open(filename_split,'w')          
               outfile.write("##sequence-region " + gff_region)
               outfile.close()
    
    return

def clean_gffs(input_dir, output_dir):
    
    """Remove non-CDS annotations and CDSs that are not divisible by 3"""
    import re 
    import glob 
    import os
    
    #Import and process annotations
    path_gff = 'contigs' + "/*.gff"
    
    files_gff = glob.glob(path_gff)
    
    files_complete = 0
    
    for file in files_gff:
        with open(file,"r") as gff:
            stored_gff = gff.read()
       
        d = {"Dbxref=": "dbxref=", "Notes=": "note="}
        stored_gff = replace_all(stored_gff, d) #Ensure GFF annotation format is consistent with prokka output
            
        gene_list = stored_gff.splitlines() 
        
        for line in range(len(gene_list)):
            if '##sequence-region' in gene_list[line]:
                title = gene_list[line]
            elif '###' in gene_list[line]:
                gene_list = gene_list[:line]
        
        gene_list = gene_list[2:]
        
        genes = []
        for y in range(len(gene_list)):
            split = re.split(r'\t+', gene_list[y])
            if split[2] == 'CDS' and ((int(split[4]) - int(split[3]) + 1) % 3) == 0: #process_pokka input only looks for CDSs and returns duplicate error when the exon is split.
                genes.append(gene_list[y])
            else:
                pass
        
        Ids= []
        for gene in range(len(genes)):
            Id = re.search('ID=(.*?);', genes[gene]).group(1)
            Ids.append(Id)
        
        output_cds = []
        for Id in range(len(Ids)):
            if not Ids.count(Ids[Id]) == 1:
                print('Removing split CDS {}'.format(Ids[Id]))
            else:
                output_cds.append(genes[Id])
        
        if len(output_cds) == 0:
            continue
        
        cleaned_gffs = "\n".join(str(z) for z in output_cds)
        cleaned_gffs = title + '\n' + cleaned_gffs
        
        #Concatenate downloaded fasta and reformatted GFF
        header = os.path.basename(file).split(".gff")[-2]
        fasta_file = 'contigs/' + header + '.fna'
        cleaned_gffs = cleaned_gffs + '\n' + '##FASTA' + '\n' + str((open(fasta_file,'r').read())) 
        filename_cleaned = output_dir + '/' + header + '.gff'
        outfile = open(filename_cleaned,'w')          
        outfile.write(cleaned_gffs)
        outfile.close()
        
        files_complete += 1
        print(str(files_complete) + "/" + str(len(files_gff)) + " complete")
    return


def main():
    start = time.time()
    args = get_options()
    
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    current_dir = os.getcwd() + '/'
    os.makedirs(str(current_dir + 'raw_files'))
    os.makedirs(str(current_dir + 'contigs'))
    os.makedirs(str(current_dir + args.output_dir))
    accessions = args.accessions.replace('"', "")
    accessions = args.accessions.split(',')
    
    print("Retrieving Sequences")
    
    if isinstance(accessions, str):
        headers = genome_downloader(args.email, accessions, int(args.number))
    elif isinstance(accessions, list): 
         headers = []
         for accession in accessions:
             headers.append(str(genome_downloader(args.email, accession, (int(args.number)+1))))#, args.complete)
   
    headers = headers[0].split(',')
    print(('found {} ids').format(len(headers)))
    
    print("Separating contigs")
    
    split_contigs(headers,'raw_files', 'contigs')
    
    print("Processing GFF files")
    
    clean_gffs('contigs', args.output_dir)
    
    print("Done")
    end = time.time()
    print(end - start)
    
    return

if __name__ == '__main__':
    
    main()
