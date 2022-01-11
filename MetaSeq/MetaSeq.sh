################################################################################################################################
#                                                   INPUT REQUIRED                                                             #
#                                                 General information                                                          #
################################################################################################################################

################################################################################################################################
#                                          Paths to third-party software and databases                                         #
################################################################################################################################
usearch=/usr/local/bin/usearch64
trimmo=/bioinf/home/leon.dlugosch/Resources/Trimmomatic/trimmomatic-0.36.jar        
rm_smalls=/bioinf/home/leon.dlugosch/Resources/RemoveSmalls/RemoveSmalls.pl
tab2fasta=/bioinf/home/leon.dlugosch/Resources/Tab2Fasta/Tab2Fasta.py

Kaiju_RefSeq_Nodes=/bioinf/home/leon.dlugosch/Resources/Kaiju_nr_2021_03/nodes.dmp
Kaiju_RefSeq_Names=/bioinf/home/leon.dlugosch/Resources/Kaiju_nr_2021_03/names.dmp
Kaiju_RefSeq=/bioinf/home/leon.dlugosch/Resources/Kaiju_nr_2021_03/nr/kaiju_db_nr.fmi

Kaiju_ProGenomes_Nodes=/bioinf/home/leon.dlugosch/Resources/Kaiju_Progenomes_2021_03/nodes.dmp
Kaiju_ProGenomes_Names=/bioinf/home/leon.dlugosch/Resources/Kaiju_Progenomes_2021_03/names.dmp
Kaiju_ProGenomes=/bioinf/home/leon.dlugosch/Resources/Kaiju_Progenomes_2021_03/progenomes/kaiju_db_progenomes.fmi

jgi_summarize=/bioinf/home/leon.dlugosch/Resources/metabat/jgi_summarize_bam_contig_depths

metabat2=/nfs/data/kila0393/metabat/metabat2

adapter=/bioinf/home/leon.dlugosch/Resources/Adapter/Adapter_new.fna
################################################################################################################################
#                                                    Analysis defaults                                                         #
################################################################################################################################
mode=0

complete_genes=0				# should only complete genes (Prodigal) be used for analysis?
cluster_method=usearch 			# Clustering algorythm. can be "usearch" or "cdhit"
no_cazy=0						# Skip cazyme calssification
no_uniref90=0					# Skip uniref calssification
no_kegg=0						# Skip splitting of aa data for GhostKOALA
no_cluster=0					# dereplicated sequences will not be clustered, if set to 1: id will be ignored
id=95							# identity threshold for clustering
no_tax=0						# Skip taxonomic calssification using kaiju and RefSeq/ProGenomes databases
prot_id=70						# Identity threshold for diamond protein similatity searches (uniref90 and cazyme)

threads=20  					# threads used for computation                                                  
mem=120							# memory for SPAdes assembly in GB
minR=100						# minimal read length in quality control of raw illumina reads
minlen=210						# minimal contig and gene length
mincov=3						# minimal average coverage of contigs
k=11,21,33,55					# kmer size for assembly

db=0							# database for mapping of reads and read abundance tables
rDir=0							# input directory
oDir=.							# output directory

################################################################################################################################
#                                                          Options                                                             #
################################################################################################################################
while : ; do
	case $1 in
	-i)
		shift
		rDir=$1
		shift
		;;
	-o)
		shift
		oDir=$1
		shift
		;;
	-t)
		shift
		threads=$1
		shift
		;;
	-minR)
		shift
		minR=$1
		shift
		;;
	-cluster_method)
		shift
		cluster_method=$1
		shift
		;;
	-no_cluster)
		shift
		no_cluster=$1
		shift
		;;
	-no_tax)
		shift
		no_tax=$1
		shift
		;;
	-no_kegg)
		shift
		no_kegg=$1
		shift
		;;
	-no_uniref90)
		shift
		no_uniref90=$1
		shift
		;;
	-minlen)
		shift
		minlen=$1
		shift
		;;
	-id)
		shift
		id=$1
		shift
		;;
	-prot_id)
		shift
		prot_id=$1
		shift
		;;
	-complete_genes)
		shift
		complete_genes=$1
		shift
		;;
	-mincov)
		shift
		mincov=$1
		shift
		;;
	-k)
		shift
		k=$1
		shift
		;;
	-mem)
		shift
		mem=$1
		shift
		;;
	-h|--help)
		echo ""
		echo "USAGE:"
		echo "        path/to/Samples/directory"
		echo "        -t [INT] Number of available threads for calculation; default: 8"
		echo "        -minR [INT] Minimal read length (R1 & R2 read). Shorter reads will be discarded; default = 100 (QC)"
		echo "        -o [PATH] Path to output directory. Subdirectories will be created."
		echo "        -h this help"
		echo ""
		echo "        path/to/Samples.fastq, path to directory containing samples."
		echo "        [INT]: integer (e.g. 1, 54)"
		echo "        [NUM]: numeric (e.g. 0.5, 0.1)"
		echo "        [option]: one of pre-defined options."
		echo ""
		echo "        For additional support, contact: leon.dlugosch>uni-oldenburg.de"
		exit
		;;
	*)  if [ -z "$1" ]; then break; fi
		mode=$1
		shift
		;;
	esac
done

minAA=$(($minlen/3))

if [[ "$mode" == 0 ]]; then
	echo "Error: no mode defined"
	exit
fi

if [[ "$rDir" == 0 ]]; then
	echo "No input directory selected."
	echo "Exiting skript."
	exit
fi

if [[ "$oDir" == . ]]; then
	echo "No output directory selected."
	echo "Output will be written in your current working directory."
	exit
fi

rename "s/-//g" $inDir/*.fastq* ### deletes "-" from strings, because they make problems
if [[ "${inDir: -1}" == "/" ]]; then
	rDir=${inDir::-1}
fi

if [[ "${oDir: -1}" == "/" ]]; then
	oDir=${oDir::-1}
fi

###################################################################################################################################
# MODE: gzip_d & gzip_c 												                                                          #
###################################################################################################################################
if [[ "$mode" == "gzip_d" ]]; then
	( cd $rDir && gzip -d *.gz )
fi

if [[ "$mode" == "gzip_c" ]]; then
	( cd $rDir && gzip *.f* )
fi
###################################################################################################################################
# MODE: QC														                                                           		  #
###################################################################################################################################
if [[ "$mode" == "QC"  || "$mode" == "complete" ]]; then

	mkdir -p $oDir/01_QC/Paired
	mkdir -p $oDir/01_QC/Unpaired
	mkdir -p $oDir/tmp

	( cd $inDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/tmp/trim_files.txt

	for s in $(cat $oDir/tmp/trim_files.txt); do
		java -jar $trimmo PE \
		-phred33 \
		-threads $threads \
		$rDir/${s}_*R1*.fastq \
		$rDir/${s}_*R2*.fastq \
		$oDir/01_QC/Paired/${s}_R1.fastq \
		$oDir/01_QC/Unpaired/${s}_SE_R1.fastq \
		$oDir/01_QC/Paired/${s}_R2.fastq \
		$oDir/01_QC/Unpaired/${s}_SE_R2.fastq \
		ILLUMINACLIP:$adapter:2:30:10:2:true \
		SLIDINGWINDOW:4:20 \
		LEADING:20 \
		MINLEN:$minR 
	done
	rm -rf $oDir/tmp/
fi

###################################################################################################################################
# MODE: assembly (only run this if you have sufficient memory - requirements vary by sample heterogeinety and sequencing depth)   #
###################################################################################################################################
if [[ "$mode" == "assembly"  || "$mode" == "complete"  ]]; then

if [[ "$mode" == "complete" ]]; then
	rDir=$oDir/01_QC/Paired
fi
	mkdir -p $oDir/02_contigs
	mkdir -p $oDir/tmp/SPAdes
	( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/tmp/Files.txt
	for s in $(cat Files_${n}.txt); do

		if [ -d $oDir/tmp/SPAdes ]; then
			rm -rf $oDir/tmp/SPAdes
			mkdir -p $oDir/tmp/SPAdes
		fi

		spades.py -1 $rDir/${s}_R1.fastq \
		-2 $rDir/${s}_R2.fastq \
		--meta \
		--memory $mem \
		-k $k \
		-o $oDir/tmp/SPAdes \
		-t $threads

		mv $oDir/tmp/SPAdes/contigs.fasta $oDir/02_contigs/${s}_contigs.fasta
	done

	( cd $oDir/02_contigs/ && ls *.f* ) > $oDir/tmp/metaquast_files.txt
	for s in $(cat $oDir/tmp/metaquast_files.txt); do
		mkdir -p $oDir/02_1_Quast_500/${s}
		metaquast.py -t $threads --max-ref-number 0 -m 500 -o $oDir/02_1_Quast_m${m}/${s} $oDir/02_contigs/${s}
	done
	rm -rf $oDir/tmp/SPAdes
fi

###################################################################################################################################
# MODE: predict														                                                           		  #
###################################################################################################################################
if [[ "$mode" == "predict" || "$mode" == "postassembly" || "$mode" == "complete" ]]; then
	if [[ "$mode" == "complete"  || "$mode" == "postassembly" ]]; then
		rDir=$oDir/02_contigs
	fi
	mkdir -p $oDir/tmp/faa
	mkdir -p $oDir/tmp/fna
	mkdir -p $oDir/03_Filtered_Contigs
	
	( cd $rDir && ls *.fasta ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $oDir/tmp/RemoveSmalls.txt
	echo "Removing Contigs smaller than "$minlen"bp..."
	for s in $(cat $oDir/tmp/RemoveSmalls.txt); do
		echo ${s}"..."
		perl $rm_smalls $minlen $rDir/${s}_contigs.fasta > $oDir/03_Filtered_Contigs/${s}_f.fasta
	done
	
	 mkdir -p $oDir/04_Genes/fna 	# Nucleotide sequences
	 mkdir -p $oDir/04_Genes/faa 	# Proteins sequneces

	 ( cd $oDir/03_Filtered_Contigs && ls *.fasta ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $oDir/tmp/Prodigal_Files.txt
	 for s in $(cat $oDir/tmp/Prodigal_Files.txt); do
	 	(FILE=$s
	 		echo $FILE" ..." 
	 		prodigal  -p meta -q -i $oDir/03_Filtered_Contigs/${FILE}_f.fasta -d $oDir/tmp/fna/${FILE}_genes.fna -a $oDir/tmp/faa/${s}_aas.faa
	 		) &

	 	if [[ $(jobs -r -p | wc -l) -gt $threads ]]; then
	 		wait -n
	 	fi
	 	sleep 1s
	 done

	 ( cd $oDir/tmp/fna/ && ls *.f* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $oDir/tmp/format_files.txt
	 for s in $(cat $oDir/tmp/format_files.txt); do
	 	fasta_formatter -i $oDir/tmp/fna/${s}.fna -o $oDir/04_genes/fna/${s}.fna
	 	fasta_formatter -i $oDir/tmp/faa/${s}.faa -o $oDir/04_genes/faa/${s}.faa
	 done
	 rm -rf $oDir/tmp/
	fi

###################################################################################################################################
# MODE: filter genes														                                               		  #
###################################################################################################################################

if [[ "$mode" == "filter_genes" || "$mode" == "postassembly" || "$mode" == "complete" ]]; then
		if [[ "$mode" == "complete"  || "$mode" == "postassembly" ]]; then
			rDir=$oDir/04_genes
		fi
	mkdir -p $oDir/tmp
	faaDir=$rDir/faa/
	fnaDir=$rDir/fna/

	if [[ "${complete_genes}" == 1 ]]; then
		mkdir $rDir/complete_fna
		mkdir $rDir/complete_faa  
		
		( cd $rDir/fna && ls *.f* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $oDir/tmp/filter_files.txt
		for s in $(cat $oDir/tmp/filter_files.txt); do
			cat $rDir/fna/${s}.fna | grep "partial=00" -A1 > $oDir/04_genes/complete_fna/${s}_complete.fna
			cat $rDir/faa/${s}.faa | grep "partial=00" -A1 > $oDir/04_genes/complete_faa/${s}_complete.faa
		done
		faaDir=$oDir/04_Genes/complete_faa
		fnaDir=$oDir/04_Genes/complete_fna
	fi

	if [[ "${complete_genes}" == 1 ]]; then
		mkdir $oDir/04_genes/complete_fna_filtered
		mkdir $oDir/04_genes/complete_faa_filtered
		faaOut=$oDir/04_genes/complete_faa_filtered
		fnaOut=$oDir/04_genes/complete_fna_filtered
	else
		mkdir $oDir/04_genes/fna_filtered
		mkdir $oDir/04_genes/faa_filtered
		faaOut=$oDir/04_genes/faa_filtered
		fnaOut=$oDir/04_genes/fna_filtered
	fi

	export faaOut
	export fnaOut
	export mincov
	export minlen

	( cd $fnaOut && ls *.fna ) > $oDir/tmp/FilterGenes_files.txt
	for s in $(cat $oDir/tmp/FilterGenes_files.txt); do
		(
			export $s
			Rscript /bioinf/home/leon.dlugosch/Resources/R_functions/GeneFilter.R
			) &

		if [[ $(jobs -r -p | wc -l) -gt $threads ]]; then
			wait -n
		fi
		sleep 1s
	done

	( cd $rDir/fna && ls *.fna ) | cut -f 1 -d '.' > $oDir/tmp/Format_Files.txt 
	for s in $(cat $oDir/tmp/Format_Files.txt ); do
		s_new=$(echo "$s" | sed 's/_gf//g')
		python $tab2fasta $fnaOut/${s}.fna 2 1 > $fnaOut/${s_new}_filtered.fna
		python $tab2fasta $faaOut/${s}.faa 2 1 > $faaOut/${s_new}_filtered.faa
	done

	rm $faaOut/*_gf*
	rm $fnaOut/*_gf*
	mkdir -p 04_genes/prot_filtered
	mkdir -p 04_genes/nuc_filtered
	mv $faaOut/*filtered* 04_genes/nuc_filtered 
	mv $fnaOut/*filtered* 04_genes/prot_filtered 
	rm -rf $fnaOut 
	rm -rf $faaOut 
	rm -rf $oDir/tmp/
fi

###################################################################################################################################
# MODE: Cluster_genes														                                               		  #
###################################################################################################################################
if [[ "$mode" == "cluster_genes" || "$mode" == "postassembly" || "$mode" == "complete" ]]; then
		if [[ "$mode" == "complete"  || "$mode" == "postassembly" ]]; then
			rDir=$oDir/04_genes/nuc_filtered
		fi

	mkdir -p $oDir/04_genes/	
	cat $rDir/*.fna > $oDir/04_genes/all_genes_nuc.fna
	
	usearch -fastx_uniques $oDir/04_genes/all_genes_nuc.fna \
	-sizeout \
	-threads $threads \
	-fastaout $oDir/04_genes/all_genes_nuc_derep.fna
	
	if [[ "$no_cluster" == "0" ]]; then
		mkdir -p 05_nr/
		if [[ "$cluster_method" == "usearch" ]]; then
			usearch -cluster_fast $oDir/05_genes/all_genes_nuc_derep.fna \
			-id 0.$id \
			-centroids $oDir/06_nr/Genes_nr${id}.fna
			gzip $oDir/04_genes/all_genes_nuc_derep.fna
		fi
		
		if [[ "$cluster_method" == "cdhit" ]]; then
			cdhit -i $oDir/05_genes/all_genes_nuc_derep.fna \
			-o $oDir/05_nr/Genes_nr${id}.fna \
			-c 0.$id \
			-T $threads \
			-M 0

			gzip $oDir/04_genes/all_genes_nuc_derep.fna
		fi
	fi
fi

###################################################################################################################################
# MODE: classification													                                                 		  #
###################################################################################################################################
if [[ "$mode" == "classification" || "$mode" == "postassembly" || "$mode" == "complete" ]]; then
		if [[ "$mode" == "complete"  || "$mode" == "postassembly" ]]; then
			if [[ "$no_cluster" == 1 ]]; then
				rdir=$oDir/05_genes/all_genes_nuc_derep.fna
			else
				rDir=$oDir/05_nr/Genes_nr${id}.fna
			fi
			
		fi	
	mkdir $oDir/06_classification/Kaiju_taxonomy/RefSeq
	mkdir $oDir/06_classification/Kaiju_taxonomy/ProGenomes
	if [[ "$no_tax" == "0" ]]; then
		echo "Starting Kaiju: RefSeq NR"
		kaiju -t $Kaiju_RefSeq_Nodes -f $Kaiju_RefSeq -z $threads -a greedy -e 5 -E 0.00001 -i $rDir -o $oDir/06_classification/Kaiju_taxonomy/RefSeq/RefSeq_Taxonomy.txt
		kaiju-addTaxonNames -t $Kaiju_RefSeq_Nodes -n $Kaiju_RefSeq_Names -i $oDir/06_classification/Kaiju_taxonomy/RefSeq/RefSeq_Taxonomy.txt -o $oDir/06_classification/Kaiju_taxonomy/RefSeq/RefSeq_names.txt -r superkingdom,phylum,class,order,family,genus,species -u
		
		echo "Starting Kaiju: ProGenomes"
		kaiju -t $Kaiju_Progenomes_Nodes -f $Kaiju_Progenomes -z $threads -a greedy -e 5 -E 0.00001 -i $rDir -o $oDir/06_classification/Kaiju_taxonomy/ProGenomes/ProGenomes_Taxonomy.txt
		kaiju-addTaxonNames -t $Kaiju_Progenomes_Nodes -n $Kaiju_Progenomes_Names -i $oDir/06_classification/Kaiju_taxonomy/ProGenomes/ProGenomes_Taxonomy.txt -o $oDir/06_classification/Kaiju_taxonomy/ProGenomes/Progenomes_names.txt -r superkingdom,phylum,class,order,family,genus,species -u
	fi
	
	if [[ "$no_kegg" == "0" || "$no_cazy" == "0" || "$no_uniref90" == "0" ]]; then
		transeq -sequence $rDir -outseq $oDir/06_classification/Genes_aa.fna -frame 1
	fi
	
	if [[ "$no_cazy" == "0" ]]; then
		mkdir $oDir/06_classification/CAZy/
		diamond blastp --more-sensitive -p $threads --id 70 -e 0.0000000001 -k 1 -d /bioinf/home/leon.dlugosch/Resources/DiamondDB_CAZy/diamond_cazy_db.dmnd -q $oDir/06_classification/Genes_aa.fna -o $oDir/06_classification/CAZy/CAZy_IDs.txt
	fi
	
	if [[ "$no_uniref90" == "0" ]]; then
		mkdir $oDir/06_classification/uniref90/
		diamond blastp --more-sensitive -p $threads --id 70 -e 0.0000000001 -k 1 -d /bioinf/home/leon.dlugosch/Resources/UniRef/diamond_db/uniref90_2021_11.dmnd -q $oDir/06_classification/Genes_aa.fna -o $oDir/06_classification/uniref90/niref90_ids.txt
	fi
	
	if [[ "$no_kegg" == "0" ]]; then
		mkdir $oDir/07_classification/GhostKOALA_parts
		perl /bioinf/home/leon.dlugosch/Resources/FastaSplitter/fasta-splitter.pl $oDir/06_classification/Genes_aa.fna --part-size 300000000 --line-length 0 --out-dir $oDir/06_classification/GhostKOALA_parts
		gzip $oDir/06_classification/GhostKOALA_parts/*
	fi
fi

###################################################################################################################################
# MODE: MG_mapping     													                                                 		  #
###################################################################################################################################
if [[ "$mode" == "MG_mapping" || "$mode" == "postassembly" || "$mode" == "complete" ]]; then

	if [[ "$mode" == "postassembly" || "$mode" == "complete" ]]; then
		rDir=$oDir/01_QC/Paired
		db=$oDir/05_nr/Genes_nr${id}.fna
	fi

	mkdir -p $oDir/tmp/bam
	mkdir -p $oDir/tmp/sam
	mkdir -p $oDir/tmp/db
	mkdir -p $oDir/07_map
	
	bowtie2-build ${db} $oDir/tmp/db/Bowtie2.db
	DB=$oDir/tmp/db/Bowtie2.db
	
	( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/tmp/files.txt
	for s in $(cat $oDir/tmp/files.txt); do
		bowtie2 --very-sensitive-local \
		-x $DB \
		-1 $rDir/${s}_R1.fastq \
		-2 $rDir/${s}_R2.fastq \
		-p $threads \
		-S $oDir/tmp/sam/${s}.sam
		
		samtools view -b -S $oDir/tmp/sam/${s}.sam > $oDir/tmp/bam/${s}.bam
		rm $oDir/tmp/${s}.sam
		
		samtools sort $oDir/tmp/bam/${s}.bam > $oDir/tmp/bam/${s}.sorted.bam 
		rm $oDir/tmp/bam/${s}.bam
		
		samtools index $oDir/tmp/bam/${s}.sorted.bam 
		samtools idxstats $oDir/tmp/bam/${s}.sorted.bam > $oDir/07_map/${s}_mapped.txt
		
		rm $oDir/tmp/bam/*.bam
		rm $oDir/tmp/bam/*.bai
	done 
fi
###################################################################################################################################
# MODE: gathering results --- only in complete and postassembly mode	                                                 		  #
###################################################################################################################################
if [[ "$mode" == "postassembly" || "$mode" == "complete" ]]; then
	mkdir -p 08_results/mapped
	mkdir -p 08_results/classification
	
	mv $oDir/07_map/*.txt 08_results/mapped
	if [[ "$no_cazy" == "0" ]]; then
		mkdir -p 08_results/classification
		mv $oDir/06_classification/CAZy/CAZy_IDs.txt 08_results/classification
	fi

	if [[ "$no_uniref90" == "0" ]]; then	
		mkdir -p 08_results/classification
		mv $oDir/06_classification/uniref90/uniref90_ids.txt 08_results/classification
	fi

	if [[ "$no_tax" == "0" ]]; then
		mkdir -p 08_results/classification
		mv $oDir/06_classification/Kaiju_taxonomy/ProGenomes/*.txt 08_results/classification
		mv $oDir/06_classification/Kaiju_taxonomy/RefSeq/*.txt 08_results/classification
	fi
	
	if [[ "$no_kegg" == "0" ]]; then
		mkdir -p 08_results/classification
		mkdir -p 08_results/GhostKOALA_parts
		mv $oDir/06_classification/GhostKOALA_parts/* 08_results/GhostKOALA_parts
	fi

	if [[ "$no_cluster" == "0" ]]; then
		mv $oDir/05_nr/Genes_nr${id}.fna 08_results
	else
		mv $oDir/05_genes/all_genes_nuc_derep.fna 08_results
	fi


fi

###################################################################################################################################
# MODE: binning     													                                                 		  #
###################################################################################################################################
if [[ "$mode" == "binning" ]]; then

	mkdir -p $oDir/tmp 
	mkdir -p $oDir/10_contig_depth
	mkdir -p $oDir/11_bins

	( cd $rDir && ls *.f* ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq > $oDir/tmp/read_files.txt 
	for s in $(cat $oDir/tmp/read_files.txt ); do
		bbmap.sh in1=$rDir/${s}_R1.fastq in2=$rDir/${s}_R2.fastq ref=$cDir/${s}_contigs.fasta out=$oDir/tmp/mapped.sam nodisk threads=${threads}
		samtools view -S -b $oDir/tmp/mapped.sam > $oDir/tmp/mapped.bam
		rm $oDir/tmp/mapped.sam
		samtools sort $oDir/tmp/mapped.bam -o $oDir/tmp/mapped_sorted.bam
		rm $oDir/tmp/mapped.bam
		$jgi_summazize $out/tmp/mapped_sorted.bam --outputDepth $oDir/10_contig_depth/${s}_depth.txt
		rm $oDir/tmp/mapped_sorted.bam
	done

	mkdir $oDir/11_bins/all_bins
	mkdir $oDir/12_bins_stats/

	( cd $cDir && ls *.fasta ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $oDir/tmp/metabat2_files.txt
	for s in $(cat $oDir/tmp/metabat2_files.txt); do
		mkdir -p $oDir/11_bins/${s}/${s}
		$metabat2 -i $cDir/${s}_contigs.fasta -a $oDir/10_contig_depth/${s}_depth.txt -o $oDir/11_bins/${s}/${s}_bin -m 1500
		mv $oDir/11_bins/${s}/${s}/* $oDir/11_bins/all_bins/
		rm -rf $oDir/11_bins/${s}/
	done

	mv $oDir/11_bins/all_bins/* $oDir/11_bins/
	rm -rf $oDir/11_bins/all_bins/
	mkdir -p $oDir/tmp/MAG_check/CheckM/02_tree
	checkm tree $oDir/11_bins/ $oDir/tmp/MAG_check/CheckM/02_tree -x fa -t $t
	mkdir -p $oDir/tmp/MAG_check/CheckM/03_markerfile
	checkm checkm lineage_set $oDir/tmp/MAG_check/CheckM/02_tree $oDir/tmp/MAG_check/CheckM/03_markerfile
	checkm taxon_set domain Bacteria $oDir/tmp/MAG_check/CheckM/03_markerfile
	mkdir -p $oDir/tmp/MAG_check/CheckM/04_output
	checkm analyze $oDir/tmp/MAG_check/CheckM/03_markerfile $oDir/11_bins/ $oDir/tmp/MAG_check/CheckM/04_output -t $threads -x fa
	checkm qa $oDir/tmp/MAG_check/CheckM/03_markerfile $oDir/tmp/MAG_check/CheckM/04_output > $oDir/12_bins_stats/CheckM_results.txt

	mkdir -p $oDir/12_bins_stats/ssu/fasta
	mkdir -p $oDir/12_bins_stats/ssu/tab
	wget https://github.com/biocore/qiime-default-reference/raw/master/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta.gz
	gzip -d 97_otus.fasta.gz
	usearch -makeudb_usearch 97_otus.fasta -wordlength 13 -output gg97.udb
	rm 97_otus.fasta
	usearch -udb2bitvec gg97.udb -output $oDir/tmp/gg97.bitvec
	rm gg97.udb
	( cd $oDir/11_bins/ && ls *.fa ) > $oDir/tmp/bin_list.txt
	for s in $(cat $oDir/tmp/bin_list.txt); do
		usearch -search_16s $oDir/11_bins/${s} -bitvec $oDir/tmp/gg97.bitvec -fastaout $oDir/12_bins_stats/ssu/fasta/16s_${s} -tabbedout $oDir/12_bins_stats/ssu/tab/${s}.txt
	done
	find $out/MAG_check/ssu/fasta/ -size 0 -print -delete

	mkdir -p $out/tmp/gtdbtk/
	gtdbtk classify_wf --genome_dir $oDir/11_bins/ --out_dir $out/tmp/gtdbtk/ -x fa --cpus $threads
	mv $out/tmp/gtdbtk/classify/*.tsv $oDir/12_bins_stats/

	for s in $(cat $oDir/tmp/bin_list.txt); do
		cat $oDir/11_bins/$s | grep -v ">" | awk 'BEGIN{a=0; c=0; g=0; t=0;} {a+=gsub("A",""); c+=gsub("C",""); g+=gsub("G",""); t+=gsub("T","");} END{print a,t,g,c}' >> $oDir/tmp/ATCG_Counts.txt
		grep ">" $b/$s | wc -l  >> $oDir/tmp/Totalseqs.txt
		grep "^>" -v $b/$s | wc -m >> $oDir/tmp/total_length.txt
	done 

	echo -e "bin\tContigs\tA\tT\tG\tC\tTotal_length" > $oDir/12_bins_stats/bin_sequence_Stats.txt
	paste -d "\t" $oDir/tmp/bin_list.txt $oDir/tmp/Totalseqs.txt  $oDir/tmp/ATCG_Counts.txt $oDir/tmp/total_length.txt >> $oDir/12_bins_stats/bin_sequence_Stats.txt
	rm -rf $oDir/tmp
fi

###################################################################################################################################
# MODE: rna_depletion     													                                              		  #
###################################################################################################################################
if [[ "$mode" == "rna_depletion" ]]; then
	mkdir -p $oDir/01_SortmeRNA/Out/ 
	mkdir -p $oDir/01_SortmeRNA/Unaligned/ 
	
	( cd $rDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $oDir/tmp/sortME_files.txt
	for s in $(cat $oDir/tmp/sortME_files.txt); do
		sortmerna --ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/bac/silva-bac-16s-id90.fasta \
		--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/bac/silva-bac-23s-id98.fasta \
		--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/euk/silva-euk-18s-id95.fasta \
		--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/euk/silva-euk-28s-id98.fasta \
		--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/5s/rfam-5.8s-database-id98.fasta \
		--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/5s/rfam-5s-database-id98.fasta \
		--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/arc/silva-arc-16s-id95.fasta \
		--ref /bioinf/home/leon.dlugosch/Resources/SortMeRNA_DBs/arc/silva-arc-23s-id98.fasta \
		--reads $rDir/${s}_R1.fastq \
		--reads $rDir/${s}_R2.fastq \
		--workdir $oDir/01_SortmeRNA/Out/ \
		--other $oDir/01_SortmeRNA/Unaligned/ \
		--threads 1:1:$threads \
		--paired_out \
		--fastx \
		-e 0.00001 \
		-v
	done
fi
