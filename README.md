# MetaSeq-Toolkit
![MetaSeq_Flowchart](images/MetaSeq.png)


In general this pipeline is intended to be a convenient way to work though large sets of metagenomic or metatranscriptiomic datasets while also retaining high analytical flexibility due to retained intermediate results that might be useful outside of the intended purpose. Downstream analysis such as statistics etc. are not included here.

## Scope
1.	Quality trimming and adapter removal of Illumina reads
2.	Assembly of Metagenomes using metaSPAdes
  2.1.	Metaquast evaluation of assembly
3.	Quality filtering of assembled contigs & predicted genes
4.	Generation of a non-redundant gene catalogue
5.	Classification of genes
  5.1.	Taxonomy (Kaiju): RefSeq, Progenomes
  5.2.	Function: KEGG (GHOSTKOLA), CAZYmes, Uniref90
6.	MAG binning (Metabat2)
7.	rRNA depletion of transcriptomes (sortmerna)
8.	Read mapping (bowtie2) 

Since the assembly is a very memory intensive process, it may be required to run this in a HPC environment if local computational power is insufficient (such as [CARL](https://uol.de/fk5/wr/hochleistungsrechnen/hpc-facilities) and requires the transfer of high quality/trimmed reads. More information on the [HPC of the Carl-von-Ossietzky University of Oldenburg](https://uol.de/fk5/wr/hochleistungsrechnen/faq-frequently-asked-questions). metaseq is modular and can be run completely with one command or using individual modules for purposes outside the intended scope of this wrapper – or simply to reanalyse data with different parameters. 
Setup (… if used outside CvO University)
In order to work properly, you will need to install dependencies and/or define paths to executables in your environment. When working on the ICBM server Rhea everything should be correct. If not, you will need to install the following third-party software:

- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Usearch (64bit)](https://www.drive5.com/usearch/)
- [Kaiju](https://kaiju.binf.ku.dk/)
- [Diamond](https://github.com/bbuchfink/diamond)
- [cdhit](https://github.com/weizhongli/cdhit)
- [Metabat2](https://kbase.us/applist/apps/metabat/run_metabat/release)
- [samtools](http://www.htslib.org/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [(Meta)SPAdes](https://cab.spbu.ru/software/meta-spades/)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [Tab2fasta](https://github.com/shenwei356/bio_scripts/blob/master/sequence/tab2fasta)
- [Removesmalls](https://github.com/burgsdorf/removesmalls/blob/master/removesmalls.pl)
- [Fasx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)

## metaseq usage
```
metaseq [module] –i [input directory] [options]
MODULES
QC			Quality and adapter clipping of illumine reads using trimmomatic [-minQ, 			-minR]
assembly		Assembly using metaSPAdes [-k, -mem]
predict_genes		Predict genes using Prodigal
filter_genes	 	Filter genes according to completeness, length and coverage thresholds 				[-mincov, -minlen, -complete_genes]
cluster_genes		Clustering of gene sequences to generate a non-redundant gene catalogue 			[-cluster_method, -id, -no_cluster]
classification	Taxonomic and functional classification of sequences [-eval, -prot_id, 	         -no_tax, -no_kegg, -no_cazy, -no_uniref90]
map	bowtie2 mapping of (HQ) reads to nucleotide database [-db]
bin	binning, evaluation and classification of MAGs
rna_depletion	Insilico depletion of rRNA reads from (meta)transcriptomic reads

OPTIONS
-i 		[PATH] directory containing input files. Varies with module.
-o 			[PATH] Path to output directory. Subdirectories will be created.
-t 			[INT] Number of available threads; default: 16
-minR 	[INT] [QC] Minimal read length (R1 & R2 read) after adapter and quality clipping. Shorter reads will be discarded; default = 100
-minQ 		[INT] [QC] Minimal average Q-score within a 4bp sliding window.
default: 20
-mem		[INT] [assembly] memory available for metaSPAdes assembly. Default 120GB
-k			[INT,INT,INT…] [assembly] k-mer size for assembly, default: 21,33,55
-complete_genes	[0/1] 
-minlen	[INT] [filter_genes] Minimal length of genes in potential amino acid sequence. default: 210 (=70 amino acids)
-mincov 		[INT] [filter_genes] Minimal coverage for genes to be kept. default: 3
-cluster_method	[usearch/cdhit] [cluster_genes] Cluster method used for clustering of nucleotide sequences. default: usearch
-id 	[INT] [cluster_genes] Threshold for nucleotide identity (%) clustering.  Ignored if –no_cluster 1. default: 95
-no_cluster		[0/1] [cluster_genes] if set to 1: do not cluster genes. default: 0
-eval	[FLOAT] [classification] E-value threshold for classification default: 0.00001
-no_cazy	[0/1] [classification] if set to 1: do not classify sequences using the CAZyme database. default: 0
-no_kegg	[0/1] [classification] if set to 1: do not split amino acid sequences in parts for GHOSTKoala classification. default: 0
-no_uniref	[0/1] [classification] if set to 1: do not classify sequences using the uniref90 database. default: 0
-no_tax	[0/1] [classification] if set to 1: do not classify sequences taxonomically using kaiju and RefSeq/ProGenomes databases. default: 0
-prot_id	[INT] [module 5] minimal amino acid sequence identity to CAZyme and UniRef90 databse. default: 70
-db	[PATH] [map] Nucleotide sequence file for mapping of reads. Only required if module is run separately.
```

