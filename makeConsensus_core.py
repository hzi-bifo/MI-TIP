#!/usr/bin/env python
import os
from os.path import join
import sys
import re
import argparse
import threading
from subprocess import call

def run(args):
    ## prepare the environment
    if not os.path.exists(args.v+'.csi'):
	os.system('bcftools index -c {}'.format(args.v))

    if os.path.exists(args.out):
	os.system('rm ' + args.out)

    regions= [l.strip() for l in open(args.r, 'r')]
    for region in regions:
	call('samtools faidx {} {} | vcf-consensus {} >> {}'.format(args.f, region, args.v, args.out), shell=True)
    
if __name__=='__main__':
    parser= argparse.ArgumentParser(
	description= 'make consensus sequences with variant calling results')
    parser.add_argument('-f', required= True, help= 'the fasta file of reference genome')
    parser.add_argument('-r', required= True, help= 'the file listing all regions to make consensus sequences')
    parser.add_argument('-out', required= True, help= 'output file')
    parser.add_argument('-v', required= True, help= 'a vcf.gz file')
    args= parser.parse_args()
    run( args)
