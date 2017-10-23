## MI-TIP: MIcrobial Tree Inference Pipeline
#### Introduction
MI-TIP is a pipeline to compute a tree of bacterial population without precomputed genomic sequences. To conduct the pipeline, only one single command is needed. All the options required by the pipeline can be edited in another file, helping to review and reproduce results with the same method and data. 
#### Dependencies
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
#### Usage 
##### 1. Check required files
- MI-TIP.config (copy and modify before running MI-TIP)
- fq list (see FQ_LIST.FORMAT for details)
- gene regions list (see GENE_REGIONS.FORMAT for details)
- reference genome (fasta format)
- fastq files
##### 2. Run the command
```
MI-TIP <MI-TIP.config>
```

#### Principle processes
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
