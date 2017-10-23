#!/usr/bin/env python
import re
import os
import argparse
import sys
import pandas as pd
import pprint

#def create_command(reads_type, fq1, fq2, insertsize, insertsd, ref, sam_f):
def create_command(info, ref, sam_f):
    reads_type= info['type']
    if reads_type == 'mate-pair':
	print('stampy.py -g {0} -h {1} --insertsize2=-{2} --insertsd2={3} -M {4} {5} > {6}'.format(ref, ref, info['insert_size'], info['insert_std'], info['reads_file1'], info['reads_file2'], sam_f))
    elif reads_type == 'paired-end':
	print('stampy.py -g {0} -h {1} --insertsize=-{2} --insertsd={3} -M {4} {5} > {6}'.format(ref, ref, info['insert_size'], info['insert_std'], info['reads_file1'], info['reads_file2'], sam_f))
    elif reads_type == 'single':
	print('stampy.py -g {0} -h {1} -M {2} > {3}'.format(ref, ref,  info['reads_file1'], sam_f))
    else:
	sys.exit('Reads type {} unrecognized'.format(reads_type))


def check_reference(ref):
    #Build a genome (.stidx) file:
    ref_genome= ref+'.stidx'
    if not os.path.exists(ref_genome):
	os.system('stampy.py -G {} {}'.format(ref,ref))
    #Build a hash (.sthash) file:
    ref_hash= ref+'.sthash'
    if not os.path.exists(ref_hash):
	os.system('stampy.py -g {} -H {}'.format(ref, ref))

    return(ref_genome, ref_hash)


if __name__== '__main__':
    parser= argparse.ArgumentParser(description='Create commands to map reads to a reference genome')
    parser.add_argument('--l', dest= 'l', required= True, help= 'list of fastq files')
    parser.add_argument('--r', dest= 'reference', required= True, help= 'reference genome (fasta format)')
    parser.add_argument('--out', dest= 'out', required= True, help= 'a directory of output sam files')

    args= parser.parse_args()
    reads_l_f= args.l
    ref_f= args.reference
    sam_d= args.out

    ref_genome, ref_hash= check_reference(ref_f)

    reads_df= pd.read_csv(reads_l_f, sep='\t', header= 0, index_col= None)
    for n in range(reads_df.shape[0]):
	info= reads_df.iloc[n, :].to_dict()
	sam_f= '{}/{}.sam'.format(sam_d, info['name'])
	create_command(info, ref_f, sam_f)

