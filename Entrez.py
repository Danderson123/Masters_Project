def genome_downloader(email, Accessions):
    from Bio import Entrez
    #Download the genomes in fasta format  
    Entrez.email = email
    search = " ".join(Accessions)

    #sequences = []
    
    handle = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
    genome_ids = handle['IdList']
    
    for genome_id in genome_ids:
        record = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        filename = '{}.fasta'.format(genome_id)
        directory = 'genomes/' + filename
        print('Writing:{}'.format(filename))
        #sequences.append(record.read())
        with open(directory, 'w') as f:
            f.write(record.read())
    
    return genome_ids #,sequences 

