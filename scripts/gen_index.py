"""
Generate a searcheable index from a panaroo output
"""
import cobs_index as cobs
import sys
import os
from tqdm import tqdm
import tempfile
import shutil
import glob 

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

def translate_alignements(alignement_path):
    """"Translate all alignement files in path and add to single list"""
    all_proteins = []
    for file in tqdm(alignement_path):
        with open(file, "r") as f:
            aln = f.read()
        gene_name = (file.split("/")[1]).split(".")[0]
        print(gene_name)
        aln = aln.split(">")[1:]
        cleaned = []
        for x in aln:
           sequence = "".join(x.splitlines()[1:]).replace("-", "")
           sequence = ">"+gene_name + "\n" + translate(sequence.upper())
           if not sequence == (">"+gene_name + "\n"):
               cleaned.append(sequence)
           
        cleaned = "\n".join(cleaned)
        all_proteins.append(cleaned)
    
    all_proteins = "\n".join(all_proteins)

    return all_proteins

def long_centroid(graph):
    """extract protein sequences associated with each node and add to a single list"""
    all_proteins = []
    G = nx.read_gml(graph)
    for node in G._node:
        y = G._node[node]
        seq = y["protein"].split(";")
        for x in seq:
            all_proteins.append(">" + y["name"] + "\n" + x)
            
    all_proteins = "\n".join(all_proteins)
    return all_proteins

def get_options(): #options for downloading and cleaning
    
    import argparse
    description = 'Generate a searcheable index from panaroo-outputted alignements'
    parser = argparse.ArgumentParser(description=description,
                                     prog='generate-index')
    io_opts = parser.add_argument_group('Inputs')
    
    io_opts.add_argument("-i",
                        "--input_dir",
                        dest="input_dir",
                        required=True,
                        help='path of panaroo-output',
                        type=str) #Specify the search term for entrez
    io_opts.add_argument("-o",
                        "--output",
                        dest="output_dir",
                        required=True,  
                        help="output directory for generated index",
                        type=str)
    io_opts.add_argument("--kmer_length",
                        dest="kmer_length",
                        help="specify kmer length for the constructed index (default = 10)",
                        default=10,
                        type=int)
     io_opts.add_argument("--source",
                      dest="srce",
                      help="build index from alignments for CD-HIT representative sequences",
                      choices=['alignment', 'centroid']
                      default="centroid",
                      type=int)
    io_opts.add_argument("--false_positive_rate",
                      dest="fpr",
                      help="false positive rate for index. Greater fpr means smaller index (default = 0.01).",
                      default=0.01,
                      type=int)
    io_opts.add_argument("--dirty",
                      dest="dirty",
                      help="write file containing translated genes",
                      action='store_true',
                      default=False)
        
    args = parser.parse_args()
    return (args)

def write_files(titles, cleaned, temp_dir):
    """Separate translated proteins and write out to individual files"""
    duplicate_set = set()
    duplicate_count = 0
    for title in range(len(titles)):
        info = cleaned[title]
        if not titles[title] in duplicate_set:
            prot_file = open(temp_dir + titles[title] + ".txt", "w")
            prot_file.write(info)
            prot_file.close()
            duplicate_set.add(titles[title])
        else:
            prot_file = open(temp_dir + titles[title] + ";" + str(duplicate_count) + ".txt", "w")
            prot_file.write(info)
            prot_file.close()
            duplicate_count += 1
    sys.exit(0)

def create_index(input_dir, output_dir, kmer_length, fpr):
    """"Create cobs index"""
    params = cobs.CompactIndexParameters()
    params.term_size = kmer_length
    params.clobber = True               # overwrite output and temporary files
    params.false_positive_rate = fpr    # higher false positive rate -> smaller index
    cobs.compact_construct(os.path.join(input_dir), os.path.join(output_dir) + str(kmer_length) + "_index_index.cobs_compact", index_params=params)
    sys.exit(0)

def main():
    args = get_options()
    
    alignement_path = glob.glob(args.input_dir + "/aligned_gene_sequences/*.aln.fas")
    alignement_path += glob.glob(args.input_dir + "/aligned_gene_sequences/*.fasta")
    
    if args.srce == "alignment":
        all_proteins = translate_alignements(alignement_path)
    if args.srce == "centroid":
        all_proteins = long_centroid(args.input_dir + "/final_graph.gml")   
    
    if args.dirty:
        out_file = open(args.output_dir + "/all_proteins.txt", "w")
        out_file.write(all_proteins)
        out_file.close()
        
    titles = []
    cleaned = []
    for x in all_proteins.split(">")[1:]:
       	cleaned.append("".join(x.splitlines()[1:]))
       	titles.append(x.splitlines()[0])
           
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    
    write_files(titles, cleaned, temp_dir)
    
    create_index(args.input_dir, args.output_dir, args.kmer_length, args.fpr)
    shutil.rmtree(temp_dir)
    
    sys.exit(0)
    
if __name__ == '__main__':
    main()
