#!/usr/bin/env python
'''
Put sequences of the same regions into a same file. 
Gene loss is also detected in this script, and for these groups it is not needed to run alignment.
'''

from os import listdir
import os
from os.path import join
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import argparse
from os import path
  
def check_dir(d):
    if not os.path.exists(d):
	os.makedirs(d)
def cluster(strain, f, output_dir):
#    strain= re.search('.+\/([^\/]+)\.fa', f).group(1)

    gene_gain= []
    for record in SeqIO.parse(f.strip(), 'fasta'):
	group= record.id
	group_f= path.join(output_dir, group+'.fa')
	seq= record.seq
	if re.search('^\w+$', str(seq)) is None:## gene loss (no reads in this region)
	    continue
	else:
	    gene_gain.append(group_f)
	    new_seq_record= SeqRecord(Seq(str(seq), generic_dna), id= strain, name= '', description= '')
	    ## print the sequence
	    with open(group_f, 'a') as handle:
	        SeqIO.write(new_seq_record, handle, 'fasta')
    return(gene_gain)

def run(files_list, out_d, com_gene_f):
    fh= open(files_list, 'r')
    files= {l.strip().split('\t')[0]:l.strip().split('\t')[1]  for l in fh.readlines()}
    #files= [l.strip()  for l in fh.readlines()]
    common_genes= []
    for strain in files:
	gene_gain= cluster(strain, files[strain], out_d)
	common_genes= (gene_gain if len(common_genes)==0 else list(set(common_genes) & set(gene_gain))) 
    fh.close()
    # write the common genes, which are not nessacerily core genes
    ch= open(com_gene_f, 'w') 
    map(lambda g: ch.write(g+"\n"), common_genes)
    ch.close()

if __name__=='__main__':
    parser= argparse.ArgumentParser(description='Collect sequences of same gene from gene files of strains')
    parser.add_argument('--l', dest= 'f_list', required= True, help= 'list of sequence files of each strain')
    parser.add_argument('--d', dest= 'g_d', required= True, help= 'folder of output fasta files')
    parser.add_argument('--c', dest= 'common_genes', required= True, help= 'output file to record common genes')
    args= parser.parse_args()

    files_list= args.f_list
    out_d= args.g_d
    check_dir(out_d)
    com_gene_f= args.common_genes
    run(files_list, out_d, com_gene_f)

