## MI-TIP: MIcrobial Tree Inference Pipeline
#### Introduction
MI-TIP is a pipeline to compute a tree of bacterial population without precomputed genomic sequences.
#### Dependencies
- <a href="#introduction">Introduction</a>
- <a href="#installation">Installation</a>
- <a href="#dependencies">Dependencies</a>
- <a href="#usage">Usage</a>
- <a href="#processes">Principle processes</a>
- <a href="#troubleshooting">What to do when the pipeline doesn't work as expected?</a>
#### Introduction<a name="introduction"></a>
MI-TIP is a pipeline to compute a tree of bacterial population without precomputed genomic sequences. To conduct the pipeline, only one single command is needed. All the options required by the pipeline can be edited in another file, helping users to review and reproduce results with the same method and data. Furthermore, the main script MI-TIP can be cloned and edited to conduct specific processes, allowing users to continue the works without rerunning all the pipeline when any unexpected result is generated. 
#### Installation<a name="installation"></a>
- step 1: Click "Clone or download", which should be found on the upper-right of github main page, and copy the URL.
- step 2: Open a terminal, go to the folder of installation, and clone the repository by the command
```
git clone https://github.com/hzi-bifo/MI-TIP
```
- step 3: Add the path to the environmental variables. If the installation directory is ```~/bin/MI-TIP```, the enviromental variable ```$PATH``` can be updated by the command
```
export PATH='~/bin/MI-TIP':$PATH
```
This command can also be inserted to the ```~/.profile``` to make the change be done automatically. 
#### Dependencies<a name="dependencies"></a>
MI-TIP, like most of other tree inference workflow, involves in a list of software. Considering the stability, a version same as listed here is strongly suggested. 
- samtools (1.3.1),bcftools (1.3.1), and htslib (1.3.1) https://github.com/samtools
Li H, A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627]
- bamtools (2.3.0) https://github.com/pezmaster31/bamtools
Derek W. Barnett et al., BamTools: a C++ API and toolkit for analyzing and managing BAM files, Bioinformatics, Volume 27, Issue 12, 15 June 2011, Pages 1691-1692
- freebayes-parallel in freebayes (v1.1.0) https://github.com/ekg/freebayes
Erik Garrison and Gabor Marth, Haplotype-based variant detection from short-read sequencing, 20 Jul 2012
- GNU parallel (20161222) http://www.gnu.org/software/parallel
O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
    ;login: The USENIX Magazine, February 2011:42-47.
- stampy.py (v1.0.31) http://www.well.ox.ac.uk/project-stampy
G. Lunter and M. Goodson.  Stampy: A statistical algorithm for sensitive and fast mapping of Illumina sequence reads. Genome Res. 2011 21:936-939.
- mafft (v7.305b) https://mafft.cbrc.jp/alignment/software/
Katoh, Misawa, Kuma, Miyata 2002 (Nucleic Acids Res. 30:3059-3066) 
MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. 
- trimal (1.2rev59) http://trimal.cgenomics.org/
trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses.
Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon. Bioinformatics 2009 25: 1972-1973.
- FastTreeMP (2.1.10 Double precision, No SSE3, OpenMP) 
Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. doi:10.1371/journal.pone.0009490.
#### Usage<a name="usage"></a>
##### 1. Check the materials
- fq list (see FQ_LIST.FORMAT for details)
- gene regions list (see GENE_REGIONS.FORMAT for details)
- reference genome (fasta format)
- fastq files
##### 2. Edit the environment and specify the material in the config file
- MI-TIP.config (copy and modify before running MI-TIP)
##### 3. Run MI-TIP
```
MI-TIP <path/of/your/MI-TIP.config>
```
An usage with ```nohup``` is recommended:
```
nohup MI-TIP <path/of/your/MI-TIP.config> &
```
#### Principle processes<a name="processes"></a>
##### 1. detect variant
```
# Create commands to run stampy
make_mapping_commands.py --l $FQ --r $REF_FASTA --out $SAM_DIR > $make_sam_commands
# Run commands parallely
parallel --retries 3 -j 30 --joblog $sam_log < $make_sam_commands
# generating six commands files, including five of running different processes and one for ordering them
make_sam2vcf_commands.py --r $REF_FASTA --s $sam_list --b $BAM_DIR --v $VCF_DIR --n $THR_NUM
bash sam2vcf.conductor.commands
```
##### 2. compute coding sequences
```
# Create commands to run makeConsensus_core.py
makeConsensus_commands.py --s makeConsensus_core.py --r $REF_FASTA --v $vcf_list --g $GENE_REGIONS --o $GENE_SEQ_DIR > $make_concensus_commands
# Run commands parallely
parallel --retries 3 -j $THR_NUM --joblog $consensus_log < $make_concensus_commands
```
##### 3. align gene sequences and conduct tree inference
```
# Cluster the consensus sequences of coding region by genes
makeGroupFasta.py --l $seqfiles_list --d $GENE_FAMILY_SEQ_DIR --c $core_genes
# Filter gene families by alignment quality
Compute_avIdent.sh $aln_list >  $aln_ident_list
cat $aln_ident_list | sort -rnk 2 | awk -v c=$aln_ident_list_CUT '$2 > c'| awk '{print $1}' > $good_genes_list # quality control
# Make a concatenated alignment 
concatenateAln.py --l  $good_genes_list --o $FINAL_ALN # make a concatenated alignment
# Compute a tree with the concatenated alignment
FastTreeMP -nt -gtr -gamma $FINAL_ALN > $FINAL_TREE
```
#### What to do when the pipeline doesn't work as expected?<a name="troubleshooting"></a>
MI-TIP is a bash script. By copying and editing it, processes can be easily conducted again. 
##### step 1: check the log file (default: tmp/MI-TIP.log)
A log file is written to help people track problems. The file includes two columns: the time stamp, and the message. For example, if the last line in the log was:
```
17:39:14        Computing consensus sequences...
```
it means that the pipeline stopped when computing the consensus sequences, suggesting that some problem happened after the message was written by MI-TIP.
##### step 2: find possible problems in the MI-TIP script
By checking the log file, it is known that some problems happened after the message was written. In this example, subsequent processes which didn't't run correctly after the last message was written are
```
makeConsensus_commands.py --s makeConsensus_core.py --r $REF_FASTA --v $vcf_list --g $GENE_REGIONS --o $GENE_SEQ_DIR > $make_concensus_commands
```
and
```
parallel --retries 3 -j $THR_NUM --joblog $consensus_log < $make_concensus_commands
```
It is reasonable to firstly check whether the commands in ```$make_concensus_commands``` are correct. Then, if there was nothing problematic, the ```$consensus_log``` can then be checked to see whether any command exit with unusual status. it can be suggested to test the single line command again to find the putative problem. A common problem can be caused by incorrect path of files. 
##### step 3: solve the problem and continue the piepline
Once the problem is solved, the pipeline can be continued to finish the rest parts. Please copy the```$MITIP_HOME/MI-TIP``` to ```./MI-TIP.copy```, like:
```
cp $MITIP_HOME/MI-TIP ./MI-TIP.copy
```
With any prefered editor, such as vim, users can disable the parts which have been done. Those parts can be deleted or encapsulated into functions. 

_Please note that commands before_
```
# Workflow begins #
```
_are strongly suggested to be unchanged._
##### step 4: ask for help if the problem cannot be solved
Please open an issue in this repository.

