#!/usr/bin/env python
import re
import argparse
import os
import os.path
import sys

def create(script, ref_f, genes_l_f, out_f, vcf_gz_f):
    return('{} -f {} -r {} -out {} -v {}'.format(script, ref_f, genes_l_f, out_f, vcf_gz_f))
def check_dir(d):
    if not os.path.exists(d):
	os.makedirs(d)

if __name__=='__main__':
    parser= argparse.ArgumentParser(description='Create commands to compute consensus gene sequences')
    
    parser.add_argument('--s', dest= 'script', required= True, help= 'the script to compute consensus sequence')
    parser.add_argument('--r', dest= 'reference', required= True, help= 'reference genome (fasta format)')
    parser.add_argument('--v', dest= 'vcf_l', required= True, help= 'list of compressed vcf files (.vcf.gz)')
    parser.add_argument('--g', dest= 'genes_l', required= True, help= 'list of gene regions')
    parser.add_argument('--o', dest= 'out', required= True, help= 'folder of output fasta file')
    args= parser.parse_args()
    core_script= args.script
    vcf_l_f= args.vcf_l
    ref_f= args.reference
    genes_l_f= args.genes_l
    out_d= args.out

    files= [f.strip()  for f in open(vcf_l_f, 'r') if not (re.search('vcf\.gz$', f) is None)]
    if len(files) < 1:# check if sam files are available
	sys.exit('no .vcf.gz file found')

    check_dir(out_d)
    all_commands= []
    for f in files:
	strain=re.sub('\.vcf\.gz$', '', f)
	strain=re.sub('.+\/', '', strain)
	out_f= out_d+'/'+strain+'.fa'
	com= create(core_script, ref_f, genes_l_f, out_f, f)
	print(com)
