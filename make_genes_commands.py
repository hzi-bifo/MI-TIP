#!/usr/bin/env python
import re
import os
import subprocess
import argparse
import pandas as pd
import sys

def make_commands(tmp_d, ref_f, genes_f, strain, fq_f1, fq_f2):
    strain_d= os.path.join(tmp_d,strain)
    if not os.path.exists(strain_d):
	os.makedirs(strain_d)

    bwa_sam_f= os.path.join(strain_d, 'bwa.sam')
    bwa_bam_f= os.path.join(strain_d, 'bwa.bam')
    bwa_sortedbam_f=os.path.join(strain_d, 'bwa.sorted.bam')
    stampy_sam_f= os.path.join(strain_d, 'stampy_with_bwa.sam')
    stampy_bam_f= os.path.join(strain_d, 'stampy_with_bwa.bam')
    stampy_sortedbam_f=os.path.join(strain_d, 'stampy_with_bwa.sorted.bam')
    stampy_vcf_f= os.path.join(strain_d, 'stampy_with_bwa.vcf')
    stampy_vcfgz_f= os.path.join(strain_d, 'stampy_with_bwa.vcf'+'.gz')
    out_f= os.path.join(strain_d, 'stampy_with_bwa.fa')

    commands= []
    commands.append('bwa mem  -v 2 -M -R \'@RG\\tID:snps\\tSM:snps\' -t 1 {} {} {}  > {}'.format(ref_f, fq_f1, fq_f2, bwa_sam_f))
    commands.append('samtools view -b -S {} -o {} -@ 1 '.format(bwa_sam_f, bwa_bam_f))
    commands.append('bamtools sort -in {} -out {}'.format(bwa_bam_f, bwa_sortedbam_f))
    commands.append('stampy.py -g {} -h {} --bamkeepgoodreads -M {} > {} '.format(ref_f, ref_f, bwa_sortedbam_f, stampy_sam_f))
    commands.append('samtools view -b -S {} -o {} -@ 1'.format(stampy_sam_f, stampy_bam_f))
    commands.append('bamtools sort -in {} -out {}'.format(stampy_bam_f, stampy_sortedbam_f))
    commands.append('samtools index {}'.format(stampy_sortedbam_f))
    commands.append('freebayes-parallel <(fasta_generate_regions.py {} 100000) 1 -f {} {} > {}'.format(ref_f+'.fai', ref_f, stampy_sortedbam_f, stampy_vcf_f))
    commands.append('bgzip -c {} > {}'.format(stampy_vcf_f,stampy_vcfgz_f))
    commands.append('bcftools index -c {}'.format(stampy_vcfgz_f))
    commands.append('makeConsensus_core.py -f {} -r {} -out {}  -v {}'.format(ref_f, genes_f, out_f, stampy_vcfgz_f))
    commands_pipe= ';'.join(commands)

    print(commands_pipe)

if __name__=='__main__':
    
    parser= argparse.ArgumentParser(description='Create commands to map reads to a reference genome')
    parser.add_argument('--l', dest= 'l', required= True, help= 'list of fastq files')
    parser.add_argument('--r', dest= 'reference', required= True, help= 'reference genome (fasta format)')
    parser.add_argument('--g', dest= 'genes', required= True, help= 'list of gene regions')
    parser.add_argument('--t', dest= 'tmp', default= 'mi-tip.out', required= False, help= 'directory to store MI-TIP outputs')
    args= parser.parse_args()

    # environement
    ref_f=args.reference
    ref_f= ref_f if re.search('^/', ref_f) else os.path.join(os.getcwd(),ref_f)
    l_f= args.l
    genes_f= args.genes
    tmp_d= args.tmp
    if not os.path.exists(tmp_d):
        os.makedirs(tmp_d)

    subprocess.call(['samtools','faidx',ref_f])
    subprocess.call(['bwa', 'index',ref_f])
    subprocess.call(['stampy.py', '-G', ref_f, ref_f])
    subprocess.call(['stampy.py', '-g', ref_f, '-H', ref_f])

    for l in open(l_f, 'r'):
	if len( l.strip().split('\t'))  == 3:
	    strain, fq_f1, fq_f2= l.strip().split('\t')
	elif len( l.strip().split('\t'))  == 2:
	    strain, fq_f1= l.strip().split('\t')
	    fq_f2= ''
	else:
	    sys.exit('{} unrecognizable'.format(l_f))
        make_commands(tmp_d, ref_f, genes_f, strain, fq_f1, fq_f2)
