#!/usr/bin/env python
from os import listdir
import os
from os import path
import sys
import re
from Bio import SeqIO
import argparse

def check_files(l_f):
    if not os.path.exists(l_f):
	sys.exit('ERROR: {} not existing\n'.format(l_f)) 
    else:
	fasta_files= [l.strip() for l in open(l_f, 'r')]
	prob_files= filter(lambda f: not os.path.exists(f), fasta_files)
	if len(prob_files) > 0:
	    sys.exit('ERROR: files not existing\n'+','.join(prob_file)+'\n') 
	else:
	    return(fasta_files)

def run (args):
    fasta_files=check_files(args.f)
    for fasta_f in fasta_files:
	group= re.search('\/([^\/]+)\.fa', fasta_f).group(1)
	aln_f= path.join(args.o, group+'.aln')
	mafft_alg= '--'+args.a
	mafft_cpu= str(args.n)
	mafft_mol= '' if args.p == '1' else '--nuc'

	print('mafft --quiet --thread {} --maxiterate 10 {} {} {} > {}'.format(mafft_cpu, mafft_mol, mafft_alg, fasta_f, aln_f))

if __name__=='__main__':
    parser= argparse.ArgumentParser(description= 'Create commands to run mafft and align sequences')
    parser.add_argument('-f', required= True, help= 'list of the fasta files to align')
    parser.add_argument('-o', required= True, help= 'the folder of output alignment results')
    parser.add_argument('-n', help= 'CPU number for mafft', type= int, default= 1)
    parser.add_argument('-p', required= True, help= 'protein (1) or nucleotide (0) sequences', default= '1', choices= ['1', '0'])
    parser.add_argument('-a', required= True, help= 'alignment algorithm for mafft', type= str, choices= ['globalpair', 'genafpair', 'localpair'], default= 'globalpair')
    args= parser.parse_args()
    run(args)

