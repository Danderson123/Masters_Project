# This script downloads genome files in fasta format for accessions of interest. 
# Fastas are then inputted into Prokka to generate GFF files and used as an input into panaroo.
# This script specifies 5 pneumnococcus reference genomes and currently is set works on mac but both of these can be adapted easily  
  
import subprocess
import os
from Entrez import genome_downloader

email = "danielanderson1@hotmail.com" #email for entrez access
#Accessions of interest
Accessions = ['NZ_CP027540.1', 'NC_017592.1', 'NC_012469.1', 'NC_011900.1', 'NZ_AKVY01000001.1', 'NC_003098.1']

current_dir = os.getcwd() + '/'
os.makedirs(str(current_dir + 'genomes'))

headers = genome_downloader(email,Accessions)

os.makedirs(str(current_dir + 'reannotated'))
os.makedirs(str(current_dir + 'gffs'))

#Run Prokka to generate GFF files and store in GFF folder
for x in headers:
    file = " " + current_dir +'genomes/' + str(x) + '.fasta'
    command = 'prokka --outdir '+ current_dir + 'reannotated/{}'.format(x) + ' --prefix ' + x + ' --force' + file
    subprocess.run(command, shell = True)
    gff = current_dir + 'reannotated/{}'.format(x) + '/' + x + '.gff'
    out_gff = current_dir + 'gffs' + '/' + x + '.gff'
    os.replace(gff, out_gff)

#Run panaroo on the prokka outputted GFFs
os.makedirs(str(current_dir + 'results'))
panaroo_command = 'panaroo -i ' + current_dir + 'gffs/' + '*.gff -o results --clean-mode strict'

subprocess.run(panaroo_command, shell = True)


