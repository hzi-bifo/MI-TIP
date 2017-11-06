#!/usr/bin/env python
import os
import os.path
import re
import argparse
import sys

def check_dir(d):
    if not os.path.exists(d):
	os.makedirs(d)
def check_reference_index(ref_f):
    ref_index_f= ref_f+'.fai'
    if not os.path.exists(ref_index_f):
        os.system('samtools faidx {}'.format(ref_f))
    return(ref_index_f)
def sam2bam_command(sam_f, bam_f):
    return('samtools view -b -S {} -o {} -@ 1'.format(sam_f, bam_f))
def sort_bam_command(bam_f, sorted_bam_f):
    return('bamtools sort -in {} -out {}'.format(bam_f, sorted_bam_f))
def bam_index_command(sorted_bam_f):
    return('samtools index {}'.format(sorted_bam_f))
def call_variance_command(ref_index_f, ref_f, sorted_bam_f, vcf_f):
    return('freebayes-parallel <(fasta_generate_regions.py {0} 100000) 1 -f {1} {2} > {3}'.format(ref_index_f, ref_f, sorted_bam_f, vcf_f))
def compress_vcf_command(vcf_f, vcf_gz_f):
    return('bgzip -c {0} > {1}'.format(vcf_f, vcf_gz_f))

def create (sam_f, ref_index_f, bam_d, vcf_d):
    re_out=re.search('(.+)\/([^\/]+)\.sam', sam_f)
    strain= ''
    sam_d= ''
    if re_out is None:
	strain= re.sub('\.sam', '', sam_f)
	sam_d= '.'
    else:
        strain= re_out.group(2)
	sam_d= re_out.group(1)
    bam_f= bam_d+'/'+strain+'.bam'
    sorted_bam_f= bam_d+'/'+strain+'.bam.sorted'
    vcf_f= vcf_d+'/'+strain+'.vcf'
    vcf_gz_f= vcf_f+'.gz'

    commands_dict= {'sam2bam':sam2bam_command(sam_f, bam_f),
	'sort_bam':sort_bam_command(bam_f, sorted_bam_f),
	'bam_index':bam_index_command(sorted_bam_f),
	'call_variance':call_variance_command(ref_index_f, ref_f, sorted_bam_f, vcf_f),
	'compress_vcf':compress_vcf_command(vcf_f, vcf_gz_f)}
    return(commands_dict)


if __name__=='__main__':
    parser= argparse.ArgumentParser(description='Create commands for detecting variants')
    parser.add_argument('--r', dest= 'reference', required= True, help= 'reference genome (fasta format)')
    parser.add_argument('--s', dest= 'sam', required= True, help= 'list of sam files (.sam)')
    parser.add_argument('--b', dest= 'bam', required= True, help= 'folder of output bam files (.bam)')
    parser.add_argument('--v', dest= 'vcf', required= True, help= 'folder of output vcf files (.vcf)')
    parser.add_argument('--n', dest= 'cpu', required= True, help= 'cpu')

    args= parser.parse_args()
    ref_f= args.reference
    sam_l_f= args.sam 
    bam_d= args.bam
    vcf_d= args.vcf
    cpu= args.cpu

    map(lambda d: check_dir(d), [bam_d, vcf_d])

    ref_index_f=check_reference_index(ref_f) # the index files required by variant calling
    
    files= [f.strip()  for f in open(sam_l_f, 'r')]
    if len(files) < 1:# check if sam files are available
	sys.exit('no .sam file found')

    all_commands= []
    for f in files:
	all_commands.append(create(f, ref_index_f, bam_d, vcf_d))
    commands_order= ['sam2bam', 'sort_bam', 'bam_index', 'call_variance', 'compress_vcf']
    commands_files= {group: group+'.commands' for group in  commands_order}
    conductor_file= 'sam2vcf.conductor.commands'
    conductor_fh= open(conductor_file, 'w')
    for command_group in commands_order: # processes with same tool are put together
	commands_fh= open(commands_files[command_group], 'w')
	for n in range(len(all_commands)):
	    commands_fh.write(all_commands[n][command_group]+'\n')
	commands_fh.close()
	conductor_fh.write('parallel -j {} --joblog {}.log --retries 3 < {}\n'.format(cpu, command_group, commands_files[command_group]))
    conductor_fh.close()
