#!/usr/bin/env python
import os
from os import path
import sys
import re
import argparse
import threading
from subprocess import call

def run(args):
    cwd=os.getcwd()
    ref= args.f if re.search('^/', args.f) else path.join(cwd, args.f)
    vcf_f= args.v if re.search('^/', args.v) else path.join(cwd, args.v)
    out_f= args.out if re.search('^/', args.out) else path.join(cwd, args.out)

    
    ## prepare the environment
    vcf_index_f= vcf_f+'.csi'
    if os.path.exists(vcf_index_f):
        os.system('rm {}'.format(vcf_index_f))
    os.system('bcftools index -c {}'.format(vcf_f))

    if os.path.exists(out_f):
	os.system('rm ' + args.out)

    regions= [l.strip() for l in open(args.r, 'r')]
    for region in regions:
	call('samtools faidx {} {} | vcf-consensus {} >> {}'.format(ref, region, vcf_f, out_f), shell=True)
#	print('samtools faidx {} {} | vcf-consensus {} >> {}'.format(ref, region, vcf_f, out_f))
	#print('samtools faidx {} {} | vcf-consensus {} >> {}'.format(args.f, region, args.v, args.out))
    
if __name__=='__main__':
    parser= argparse.ArgumentParser(
	description= 'make consensus sequences with variant calling results')
    parser.add_argument('-f', required= True, help= 'the fasta file of reference genome')
    parser.add_argument('-r', required= True, help= 'the file listing all regions to make consensus sequences')
    parser.add_argument('-out', required= True, help= 'output file')
    parser.add_argument('-v', required= True, help= 'a vcf.gz file')
    args= parser.parse_args()
    run( args)
