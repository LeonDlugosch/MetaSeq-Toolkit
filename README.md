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
- [metaquast](http://bioinf.spbau.ru/metaquast)

## metaseq usage
```
metaseq [module] –i [input directory/file] [options]
```
### MODULES
```
qc              Quality and adapter clipping of illumina reads using trimmomatic [-minQ, -minR]
assembly        Assembly using metaSPAdes [-k, -mem]
predict_genes   Predict genes using Prodigal [-minlen]
filter_genes    Filter genes according to completeness, length and coverage thresholds [-mincov, -minlen, -complete_genes]
cluster_genes   Clustering of gene sequences to generate a non-redundant gene catalogue [-cluster_method, -id, -no_cluster]
classify        Taxonomic and functional classification of sequences [-eval, -prot_id, -no_tax, -no_kegg, -no_cazy, -no_uniref90]
map             bowtie2 mapping of (HQ) reads to nucleotide database [-db]
bin             binning, evaluation and classification of MAGs
rna_depletion   in silico depletion of rRNA reads from (meta)transcriptomic reads
complete        starts with illimina reads and includes QC, assembly, predict_genes, filter_genes, cluster_genes, classification and map
post_assembly   starts with assembled contigs and includes predict_genes, filter_genes, cluster_genes, classification and map
```

### OPTIONS
```
-i              [PATH] directory containing input files. Varies with module.
-o              [PATH] Path to output directory. Subdirectories will be created.
-t              [INT] Number of available threads; default: 16
-minR           [INT] [QC] Minimal read length (R1 & R2 read) after adapter and quality clipping. Shorter reads will be discarded; default = 100
-minQ           [INT] [QC] Minimal average Q-score within a 4bp sliding window. default: 20
-mem            [INT] [assembly] memory available for metaSPAdes assembly. Default 180GB
-k              [INT,INT,INT…] [assembly] k-mer size for assembly, default: 21,33,55
-complete_genes [0/1] 
-minlen         [INT] [filter_genes] Minimal length of genes in potential amino acid sequence. default: 210 (=70 amino acids)
-mincov         [INT] [filter_genes] Minimal coverage for genes to be kept. default: 3
-cluster_method [usearch/cdhit] [cluster_genes] Cluster method used for clustering of nucleotide sequences. default: usearch
-id             [INT] [cluster_genes] Threshold for nucleotide identity (%) clustering.  Ignored if –no_cluster 1. default: 95
-no_cluster     [0/1] [cluster_genes] if set to 1: do not cluster genes. default: 0
-eval           [FLOAT] [classification] E-value threshold for classification default: 0.00001
-no_cazy        [0/1] [classification] if set to 1: do not classify sequences using the CAZyme database. default: 0
-no_kegg        [0/1] [classification] if set to 1: do not split amino acid sequences in parts for GHOSTKoala classification. default: 0
-no_uniref      [0/1] [classification] if set to 1: do not classify sequences using the uniref90 database. default: 0
-no_tax         [0/1] [classification] if set to 1: do not classify sequences taxonomically using kaiju and RefSeq/ProGenomes databases. default: 0
-prot_id        [INT] [module 5] minimal amino acid sequence identity to CAZyme and UniRef90 database. default: 70
-db             [PATH] [map] Nucleotide sequence file for mapping of reads. Only required if module is run separately.
```

The pipeline or rather certain modules (assembly, cluster_genes, classify, map) may take a considerable amout of time to complete. I advise you to execute this in a _virtual terminal_ using e.g.
```
tmux
```
### Module usage examples
#### complete
This executes modules QC, assembly, predict_genes, filter_genes, cluster_genes, classification and map in sequence and is probably the easiest way to get results quick. For detailed module descriptions see below. 

```
metaseq complete –i /path/to/paired/illumina_reads –o output/directory –minR 100 –minQ 20 –k 21,33,55,77 –mem 180 -minlen 210 -mincov 3 -cluster_method usearch -id 95 -prot_id 70 -eval 0.00001 –t 16
```
When finished, gene abundance tables for each sample, gene classification, gene catalogue (fasta) can be found in 08_results. Binning is not part of the standart analysis and has to be run seperately. 

#### qc
Although sequences obtained from modern Illumina machines are quite good, some form of quality control is still advised as erroneous sequences and sequences containing adapters may introduce errors in the assembly. Sequences for widely used adapers is provided [here](https://github.com/LeonDlugosch/MetaSeq-Toolkit/blob/main/data/Adapter.fna). For identification of forward and reverse reads samples should be named *Filename_R1.fastq* and *Filename_R2.fastq*

#### Usage: 
```
metaseq qc –i /path/to/paired/illumina/reads –o output/directory –minR 100 –minQ 20 –t 16
```
This will use Trimmomatic to cut remaining adapter sequence fragments from reads and cuts read ends with low quality. Read pairs, in which at least one read is shorter (after quality trimming and adapter removel) than minR are discared. Don’t overdo it here, metaSPADes uses [BaysHammer](https://link.springer.com/article/10.1186/1471-2164-14-S1-S7) read error correction tool prior to assembly to further mitigate sequencing errors. Trimmed reads will be stored in outputDirectory/01_QC/Paired.These need to be transferred to the HPC CARL for assembly.

#### assembly
If HPC clusters are required for assembly trimmed reads need to be transferred to the HPC environment. [See here for details](). Using the HPC can improve speed significantly due to parallelization. 

metaSPAdes uses deBrujin-graphs and variable kmer sizes for assembly. Depending on number of samples, sequencing depth and sample heterogeinety, this may take a considerable amount of time. Get a coffe... no, not from your office coffee machine but maybe in Brittany - it's pretty nice all year long :)
In addition metaquast statistics are calculated for each assembly.  

##### Usage: 
```
metaseq assembly –i /path/to/paired/illumina/trimmed/reads –o output/directory –k 21,33,55,77 –mem 180 –t 16
```
Results are saved in outDir/02_contigs and outDir/02_1_Quast_500.

#### predict_genes
Prodigal uses a combination of GC content, ribisomal binding sites, hexamer statistics and Start/Stop codons to dertermine where a gene starts, ends and what the correct frame is. Prodigal is run in meta-mode and generates files for amino acid and nucleotide sequence. Contigs smaller than the -minlen threshold are discarded. 

##### Usage: 
```
metaseq predict_genes –i /path/to/contigs –o output/directory -minlen 210 –t 16
```
Results are saved in outDir/03_filtered_contigs and outDir/04_genes/fna for nucleotides sequences and outDir/04_genes/faa for amino acid sequences.

#### filter_genes
Gene sequences are filtered according to –minlen and –mincov thresholds and renamed according to sample ID, contig and gene number. Even though this is not entirely necessary for the analysis, it may come in helpful in analysis outside the intentioned functionality (e.g. contig taxonomy, operons…). This modules an [R-skript](https://github.com/LeonDlugosch/MetaSeq-Toolkit/blob/main/scripts/genefilter.R). … and yes, I‘m sure there is a more efficient way to do that, but it’s not excruciatingly slow and it works fine :) 
##### Usage: 
```
metaseq filter_genes -i path/to/04_genes –o output/directory -minlen 210 -mincov 3 –t 16
```
Results are saved in outDir/03_filtered_contigs and outDir/04_genes/nuc_filtered for nucleotides sequences and outDir/04_genes/prot_filtered for amino acid sequences.

#### cluster_genes
Depending on sample number, sequencing depth and sample heterogeneity, metagenomes may contain from hundreds of thousands to some million sequences some of which will be duplicates or at least very close relatives of the same template sequences. However, therein lies a bit of a problem: clustering a few hundred thousand or even a few million sequences may take some time (hours to days), but in really, _really_ large datasets (>20M sequences), computation gets a bit out of hand. Each new sequences will have to be search against an ever growing database to a point where this process soft-cappes the clustering speed. In extreme cases it might be neccessary to skip clustering and use dereplicated sequences (-no_cluster option). 

##### Usage: 
```
metaseq cluster_genes -i path/to/04_genes/nuc_filtered –o output/directory -cluster_method usearch -id 95 –t 16
```
or
```
metaseq cluster_genes -i path/to/04_genes/nuc_filtered –o output/directory -no_cluster 1 -t 16
```
Dereplicated or clustered sequences are saved in outDir/05_nr.

#### classify
Of course, metagenome are most useful if function and taxonomy of sequences are known. Nucleotide based taxonomic sequence classification is done using Kaiju in combination with RefSeq and the ProGenomes databases. Functional classification using the GHOSTKoala online tool* and uniref90 and CAZyme database are optional (diamond blastp _--more-sensitive_).
```
metaseq classify -i path/to/GeneCatalogue.fasta –o output/directory -prot_id 70 -eval 0.00001 -t 16
```
or (if no functional classification is required)
```
metaseq classify -i path/to/GeneCatalogue.fasta –o output/directory -no_kegg 1 -no_uniref90 1 -no_cayz -eval 0.00001 -t 16
```

_*requires upload of partitioned amino (< 300 mb) acid sequence data (outDir/06_classification/GHOSTKOALA_splits) to https://www.kegg.jp/ghostkoala/ . Select 'genus_prokaryotes + family_eukaryotes + viruses' database. Only one job is allowed per registered e-mail address. Use multiple to queue up some jobs. This takes about 24h per sample._ 


#### map
To estimate gene abundance, HQ reads are mapped to the gene catalogue using bowtie2 in _--very-sensitive-local_ mode and converted to gene abundance tables using samtools.
```
metaseq map -i path/to/trimmed_reads -db pat/to/gene_catalogue –o output/directory -t 16
```
or (if no functional classification is required)
```
metaseq classify -i path/to/GeneCatalogue.fasta –o output/directory -no_kegg 1 -no_uniref90 1 -no_cayz -eval 0.00001 -t 16
```

#### bin
This module uses metabat2 to create metagenomic bins of representative organisms based on contig coverage and g+c content (... and a few more things. [See here](https://peerj.com/articles/1165/)). Resulting bins are evaluated using CheckM, 16S sequences are extracted using usearch _-search_16s_ and classified with gtbd-tk _-classify_wf_.

```
metaseq bin -i path/to/trimmed_reads -c path/to/corresponding/contigs -o output/directory -t 16
```
Resulting bins, taxonomic classification and bin-statistics are saved in output/directory/11_bins and output/directory/12_bin_stats
