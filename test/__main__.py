#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script downloads fasta and gff files in fasta format for accessions of interest. These are then joined and cleaned for direct input into panaroo.
"""

  
import os

from Entrez import fasta_downloader, gff_downloader
from clean_gffs import clean_gffs

email = "danielanderson1@hotmail.com" #email for entrez access
#Accessions of interest
Accessions = ['NZ_CP027540.1', 'NC_003098.1', 'NC_017592.1', 'NC_012469.1', 'NC_011900.1', 'NZ_AKVY01000001.1']
NZ_CP027540.1,NC_003098.1,NC_017592.1,NC_012469.1,NC_011900.1,NZ_AKVY01000001.1
current_dir = os.getcwd() + '/'
os.makedirs(str(current_dir + 'genomes'))
os.makedirs(str(current_dir + 'gffs_cleaned'))
os.makedirs(str(current_dir + 'gffs'))

#Specify species for all genomes or accessions for a limited number 
#Specify if only complete genomes 
#Specify email for entrez
#Specify ouytput directory for cleaned GFFs. 
#Store fasta option.

def main():
    
    args = get_options()
    
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    
    
    headers = fasta_downloader(email, Accessions)#, args.complete)
        
    headers['gff_filenames'] = gff_downloader(email, Accessions)
    
    headers.apply(lambda row: clean_gffs(row["genome_ids"], row["gff_filenames"], 'genomes', 'gffs', 'gffs_cleaned'), axis = 1) #Count the number of guide RNA off targets
    
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
                        help="sequences to download, specify a species for all accessions or specific accessions of interest.",
                        type=str) #Specify the search term for entrez
     
    io_opts.add_argument("--complete",
                          dest="complete",
                          help="only retrieve complete genomes",
                          action='complete',
                          default=None) #Option to only retrieve complete genomes  
    
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
    
    io_opts.add_argument("--store",
                          dest="store",
                          help="store fastas in output directory",
                          action='store',
                          default=None) #Option to only retrieve complete genomes  
    
    args = parser.parse_args()
    
    return (args)

if __name__ == '__main__':
    
    main()
