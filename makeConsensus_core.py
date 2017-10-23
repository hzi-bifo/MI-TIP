#!/usr/bin/env python
import os
from os.path import join
import sys
import re
import argparse
import threading
import subprocess

class thr(threading.Thread):
    def __init__(self, sam_args):
	threading.Thread.__init__(self)
	self.ref= sam_args['ref']
	self.region= sam_args['region']
	self.vcf= sam_args['vcf_f']
	self.out= sam_args['output']
	self.tmp_out= self.vcf+'-'+self.region+'.tmp'
	#self.append= sam_args['append']
	
    def run(self):
	trial=0
	outcome=1
	max_trial= 10
	while (trial < max_trial)&(outcome != 0):
	    outcome=os.system('/bin/bash -c \"samtools faidx {} {} | bcftools consensus {} > {}\"'.format(self.ref, self.region, self.vcf, self.tmp_out))
	    trial+= 1

	thr_lock.acquire()
	os.system('cat {} >> {}'.format(self.tmp_out, self.out))
	os.system('rm {}'.format(self.tmp_out))
	thr_lock.release()

def make_threading (threads, materials):
    n= 0
    while threads[n].isAlive():
	n= (n+1)%len(threads)
    new_thr= thr(materials)
    new_thr.start()
    threads[n]= new_thr
    
class blank_thr(threading.Thread):
    def __init__(self):
	threading.Thread.__init__(self)
    def run(self):
	return

def init_threads (cpu_n):
    thr_stack= []
    for n in range(cpu_n):
	new= blank_thr()
	new.setName('thread '+str(n))
	new.start()
	thr_stack.append(new)
    return (thr_stack)

def run(threads, args):
    ## prepare the environment
    if not os.path.exists(args.v+'.csi'):
	os.system('bcftools index -c {}'.format(args.v))

    if os.path.exists(args.out):
	os.system('rm ' + args.out)
    print(os.getcwd())
    print(args.r)
    regions= open(args.r, 'r')
    for l in regions.readlines():
	#fields= re.split('\t', l)
	#target_region= fields[0]
	#samtools_args= {'ref': args.f, 'region': target_region, 'output': args.out, 'vcf_f': args.v, 'append': append}
	#append= '1'
	target_region= l.strip()
	samtools_args= {'ref': args.f, 'region': target_region, 'output': args.out, 'vcf_f': args.v}
	## threading here
	make_threading(threads, samtools_args)

    for t in threads:
	t.join()
    
if __name__=='__main__':
    parser= argparse.ArgumentParser(
	description= 'make consensus sequences with variant calling results')
    parser.add_argument('-f', required= True, help= 'the fasta file of reference genome')
    parser.add_argument('-r', required= True, help= 'the file listing all regions to make consensus sequences')
    parser.add_argument('-out', required= True, help= 'output file')
    parser.add_argument('-v', required= True, help= 'a vcf.gz file')
    parser.add_argument('-n', help= 'cpu/threads number', type= int, default= 1)
    args= parser.parse_args()
    threads_array= init_threads(args.n)
    thr_lock = threading.Lock()
    run(threads_array, args)
