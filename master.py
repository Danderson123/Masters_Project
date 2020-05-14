import subprocess
import os
from Entrez import genome_downloader

email = "danielanderson1@hotmail.com"
Accessions = ['NZ_CP027540.1', 'NC_017592.1', 'NC_012469.1', 'NC_011900.1', 'NZ_AKVY01000001.1', 'NC_003098.1']

current_dir = os.getcwd() + '/'
os.makedirs(str(current_dir + 'genomes'))

headers = genome_downloader(email,Accessions)

os.makedirs(str(current_dir + 'reannotated'))
os.makedirs(str(current_dir + 'gffs'))

#Run Prokka to generate GFF files
for x in headers:
    file = " " + current_dir +'genomes/' + str(x) + '.fasta'
    command = 'prokka --outdir '+ current_dir + 'reannotated/{}'.format(x) + ' --force' + file
    subprocess.run(command, shell = True)


