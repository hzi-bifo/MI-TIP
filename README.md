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
- step 3: Add the path to the environmental variables. If the installation path is ```~/bin/MI-TIP```, the enviromental variable ```$PATH``` can be updated by the command
```
export PATH='~/bin/MI-TIP':$PATH
```
This command can also be inserted to the ```~/.profile``` to make the change be done automatically. 
#### Dependencies<a name="dependencies"></a>
MI-TIP, like most of other tree inference workflow, involves in a list of software. Versions listed below were tested. 
- samtools (1.3.1),bcftools (1.3.1), and htslib (1.3.1) https://github.com/samtools
- bamtools (2.3.0) https://github.com/pezmaster31/bamtools
- freebayes-parallel in freebayes (v1.1.0) https://github.com/ekg/freebayes
- GNU parallel (20161222) http://www.gnu.org/software/parallel
- bwa (Version: 0.7.15-r1140) http://bio-bwa.sourceforge.net/
- stampy.py (v1.0.31) http://www.well.ox.ac.uk/project-stampy
- VCFtools (0.1.15) https://vcftools.github.io/downloads.html
- mafft (v7.305b) https://mafft.cbrc.jp/alignment/software/
- trimal (1.2rev59) http://trimal.cgenomics.org/
- FastTreeMP (2.1.10 Double precision, No SSE3, OpenMP) 
#### Usage<a name="usage"></a>
##### 1. Check the materials
- fq list (see FQ_LIST.FORMAT for details)
- reference genome (.fasta and .gff)
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
##### 1. detect variant and compute consensus sequences
```
make_GeneRegions.py --g $REF_GFF --f gene > $GENE_REGIONS
##### 2. detect variant and compute consensus sequences
```
make_genes_commands.py --t $TMP_DIR --l $FQ --r $REF_FASTA --g $GENE_REGIONS  > $mapping_commands
parallel --retries 3 -j $THR_NUM --joblog $mapping_log < $mapping_commands
```
##### 3. sort gene sequences by gene family
```
ls $TMP_DIR/*/stampy_with_bwa.fa | awk -F'/' '{print $(NF-1)"\t"$0}' > $seqfiles_list 
makeGroupFasta.py --l $seqfiles_list --d $GENE_FAMILY_SEQ_DIR --c $core_genes
```
##### 4. align gene sequences and conduct tree inference
```
# Cluster the consensus sequences of coding region by genes
makeGroupAln.py -f $core_genes -o $GENE_FAMILY_ALN_DIR -p 0 -a globalpair > $mafft_commands 
parallel -j $THR_NUM --retries 5 --joblog $mafft_log < $mafft_commands
ls $GENE_FAMILY_ALN_DIR/* > $aln_list
# Filter gene families by alignment quality
Compute_avIdent.sh $aln_list >  $aln_ident_list
cat $aln_ident_list | sort -rnk 2 | awk -v c=$aln_ident_list_CUT '$2 > c'| awk '{print $1}' > $good_genes_list # quality control
# Make a concatenated alignment 
concatenateAln.py --l  $good_genes_list --o $FINAL_ALN # make a concatenated alignment
# Compute a tree with the concatenated alignment
FastTreeMP -nt -gtr -gamma $FINAL_ALN > $FINAL_TREE
```
#### What to do when the pipeline doesn't work as expected?<a name="troubleshooting"></a>
##### step 1: check the log file (default: $TMP_DIR/MI-TIP.log)
MI-TIP writes a general log file as well as files for each part. When the output looks unusual, please go to the output folder (eg. $TMP_DIR) and check MI-TIP.log or the other log files to find putative causal comments.
##### step 2: find possible problems and edit the MI-TIP script or MI-TIP.config
MI-TIP is a bash script. By copying and editing it, the pipeline can be easily resumed. 
If a problem is found and solved, users can copy the script and remove commands alreaady done. The edited MI-TIP script can then be run with the the same or edited MI-TIP.config file. Please ensure the correct verion of both MI-TIP script and options in MI-TIP.config before running. 
##### step 3: report and feedback 
Please open an issue in the github repository and we will try to help you to solve problems and provide technical supports. Because software can always be improved, we also encourage users to give us any constructive comment.

