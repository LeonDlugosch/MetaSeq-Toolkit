################################################################################################################################
#                                                   INPUT REQUIRED                                                             #
#                                                 General information                                                          #
################################################################################################################################
Version=0.2
# IMPORTANT
# Sample-file naming convention:
# Names for the file must not contain any whitespaces. Use '_' instead. Use a common identifier for all Sample files in the pipeine!
# Example: Foreward: Sample_x_R1.fastq & Sample_y_R1.fastq Reverse: Sample_x_R2.fastq & Sample_y_R2.fastq

# If the directory-structure stays as it is, only "WD" has to be defined in the general part.

Tools=/home/leon/Desktop/Tools                                                   # Path to thrid-party tools directory
usearch=/usr/local/bin/usearch64                                                 # Define path to the Usearch executable
trimmo=/bioinf/home/leon.dlugosch/Resources/Trimmomatic/trimmomatic-0.36.jar     # Path to Trimmomatic          
rm_smalls=/bioinf/home/leon.dlugosch/Resources/RemoveSmalls/RemoveSmalls.pl
tab2fasta=/bioinf/home/leon.dlugosch/Resources/Tab2Fasta/Tab2Fasta.py

Kaiju_RefSeq_Nodes=/bioinf/home/leon.dlugosch/Resources/Kaiju_nr_2021_03/nodes.dmp
Kaiju_RefSeq_Names=/bioinf/home/leon.dlugosch/Resources/Kaiju_nr_2021_03/names.dmp
Kaiju_RefSeq=/bioinf/home/leon.dlugosch/Resources/Kaiju_nr_2021_03/nr/kaiju_db_nr.fmi

Kaiju_ProGenomes_Nodes=/bioinf/home/leon.dlugosch/Resources/Kaiju_Progenomes_2021_03/nodes.dmp
Kaiju_ProGenomes_Names=/bioinf/home/leon.dlugosch/Resources/Kaiju_Progenomes_2021_03/names.dmp
Kaiju_ProGenomes=/bioinf/home/leon.dlugosch/Resources/Kaiju_Progenomes_2021_03/progenomes/kaiju_db_progenomes.fmi 

echo Modules of the Metagenome assembly pipeline:
echo 1: pre-assembly - Quality and adapter trimming of illumina reads
echo 2: post-assembly - contig filtering and gene prediction
echo 3: post-assembly - Gene filtering
echo 4: post-assembly - Dereplication and clustering
echo 5: post-assembly - Classification (Kaiju: ProGenomes/RefSeq nr and Diamond: CAZymes)
echo 10: steps 2-5

echo " "
echo Assembly requires the HPC cluster
echo Run module: 
read STEP

################################################################################################################################
#                                                    Library preparation                                                       #
################################################################################################################################
#Defining Pipeline defaults
threads=30                                                     # Available threads
minR=100                                                        # Minimal readlength of single reads for QC
inDir=0
outdir=.
id=95
relabel=Uniq_
minAA=70
minNuc=$(($minAA*3))
minCov=2
# COMMANDLINE OPTIONS
while : ; do
    case $1 in
    	-i)
            shift
            inDir=$1
            shift
            ;;
        -t)
            shift
            threads=$1
            shift
            ;;
        -a)
            shift
            ADAPTER=$1
            shift
            ;;
        -o)
            shift
            outDir=$1
            shift
            ;;
        -id)
            shift
            id=$1
            shift
            ;;
        -minR)
            shift
            minR=$1
            shift
            ;;
        -relabel)
            shift
            relabel=$1
            shift
            ;;
        -minF)
            shift
            minF=$1
            shift
            ;;
        -minCov)
            shift
            minCov=$1
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
            inDir=$1
            shift
            ;;
    esac
done

if [[ "$inDir" == 0 ]]; then
echo "No input directory selected."
echo "Exiting skript."
fi
if [[ "$outDir" == . ]]; then
echo "No output directory selected."
echo "Output will be written in your current working directory."
fi

rename "s/-//g" $inDir/*.fastq* ### deletes "-" from strings, because they make problems
if [[ "${inDir: -1}" == "/" ]]; then
inDir=${inDir::-1}
fi

if [[ "${outDir: -1}" == "/" ]]; then
outDir=${outDir::-1}
fi



  if [[ "$STEP" == 1 ]]; then
mkdir -p $outDir/01_QC/Paired
mkdir -p $outDir/01_QC/Unpaired
mkdir -p $outDir/Temp

     ( cd $inDir && ls *.fastq ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' | uniq -d > $outDir/Temp/TrimFiles.txt
     cat $outDir/Temp/Trimfiles.txt
       for s in $(cat $outDir/Temp/TrimFiles.txt); do
         java -jar $trimmo PE \
          -phred33 \
          -threads $threads \
          $inDir/${s}_*R1*.fastq \
          $inDir/${s}_*R2*.fastq \
          $outDir/01_QC/Paired/${s}_R1.fastq \
          $outDir/01_QC/Unpaired/${s}_SE_R1.fastq \
          $outDir/01_QC/Paired/${s}_R2.fastq \
          $outDir/01_QC/Unpaired/${s}_SE_R2.fastq \
          ILLUMINACLIP:$ADAPTER:2:30:10:2:true \
          SLIDINGWINDOW:4:15 \
          LEADING:3 \
          MINLEN:$minR 
         done
     rm $outDir/Temp/TrimFiles.txt

    echo "Finished Trimming - compressing original data"
    gzip $inDir/*.fastq
  fi

if [[ "$STEP" == 2 || "$STEP" == 10 ]]; then

################################################################################################################################
#                                                    Removing small contigs                                                    #
################################################################################################################################
 mkdir -p $outDir/Temp
 mkdir -p $outDir/04_Filtered_Contigs
( cd $inDir && ls *.fasta ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $outDir/Temp/RemoveSmalls.txt
 echo "Removing Contigs smaller than "$minNuc"bp..."
   for s in $(cat $outDir/Temp/RemoveSmalls.txt); do
      echo ${s}"..."
      perl $rm_smalls $minNuc $inDir/${s}_contigs.fasta > $outDir/04_Filtered_Contigs/${s}_f.fasta
   done

################################################################################################################################
#                                                    Prodigal ORF-prediction                                                   #
################################################################################################################################

 mkdir -p $outDir/05_Genes/fna 	# Nucleotide sequences
 mkdir -p $outDir/05_Genes/faa 	# Proteins sequneces
 mkdir -p $outDir/Temp

( cd $outDir/04_Filtered_Contigs && ls *.fasta ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $outDir/Temp/Prodigal_Files.txt
   for s in $(cat $outDir/Temp/Prodigal_Files.txt); do
      (FILE=$s
      echo $FILE" ..." 
      prodigal  -p meta -q -i $outDir/04_Filtered_Contigs/${FILE}_f.fasta -d $outDir/05_Genes/fna/${FILE}_genes.fna -a $outDir/05_Genes/faa/${s}_aas.faa
      ) &
     
      if [[ $(jobs -r -p | wc -l) -gt $threads ]]; then
         wait -n
     fi
     sleep 2s
     done
   fi 

if [[ "$STEP" == 3  || "$STEP" == 10 ]]; then
mkdir -p $outDir/05_Genes/fna_filtered  # Nucleotide sequences
mkdir -p $outDir/05_Genes/faa_filtered  # Proteins sequneces
mkdir -p $outDir/Temp
( cd $outDir/05_Genes/fna && ls *.fna ) | awk 'BEGIN{FS=OFS="_"}{NF--; print}' > $outDir/Temp/Files.txt

echo "Removing Genes and Proteins smaller than "$minNuc"bp..."
    for s in $(cat $outDir/Temp/Files.txt); do
     echo $s" ..."
     perl $rm_smalls $minNuc $outDir/05_Genes/fna/${s}_genes.fna > $outDir/05_Genes/fna_filtered/${s}_${minNuc}_na.fna
     perl $rm_smalls $minAA $outDir/05_Genes/faa/${s}_aas.faa > $outDir/05_Genes/faa_filtered/${s}_${minAA}_AAs.faa
    done
  
mkdir -p $outDir/06_TabularFiles/fna
mkdir -p $outDir/06_TabularFiles/faa
echo "Transforming .fasta files to tabular format..."

( cd $outDir/05_Genes/fna_filtered/ && ls *.fna ) | cut -f 1 -d '.' | cut -f 1 -d "_"> $outDir/Temp/Format_Files.txt
  for s in $(cat $outDir/Temp/Format_Files.txt); do
    echo ${s}" nucleotide sequences..." 
    fasta_formatter -i $outDir/05_Genes/fna_filtered/${s}_${minNuc}_na.fna -o $outDir/06_TabularFiles/fna/${s}_tab.fna -t
    echo ${s}" amino acid sequences..."
    fasta_formatter -i $outDir/05_Genes/faa_filtered/${s}_${minAA}_AAs.faa -o $outDir/06_TabularFiles/faa/${s}_tab.faa -t
    echo "=================================="
  done  

export outDir
export minCov
echo "Running R-Script, to rename Sequences in .fna files ..."
Rscript /bioinf/home/leon.dlugosch/Resources/R_functions/FASTA_renamerNuc.R
echo "Running R-Script, to rename Sequences in .faa files ..."
Rscript /bioinf/home/leon.dlugosch/Resources/R_functions/FASTA_renamerAA.R

##### Executed R-Script
# library(stringr)
# path = Sys.getenv("outDir")
# path = paste(path, "/06_TabularFiles/fna/")
# setwd(path)
# file.names = dir(path, pattern = ".fna")
# minCov = Sys.getenv("minCov")
# 
# for (i in 1:length(file.names)) {
#   file.n = str_split(file.names[i], "_")
#   Cruise = unlist(lapply(file.n, "[[", 1))
#   Station = unlist(lapply(file.n, "[[", 2))
#   MG_Name = paste(Cruise, Station, sep = "_")
#   print(MG_Name)
#   
#   fasta.tab = read.delim(file.names[i],
#                          header = F,
#                          stringsAsFactors = F)
#   names(fasta.tab) = c("ID", "Sequence")
#   fasta.tab$ID = str_replace(fasta.tab$ID, " ", "_")
#   split = str_split(fasta.tab[, 1], "_")
#   ContigNr = unlist(lapply(split, "[[", 2))
#   GeneNR = unlist(lapply(split, "[[", 7))
#   
#   fasta.tab$Coverage = as.numeric(unlist(lapply(split, "[[", 6)))
#   fasta.tab$ID_short = paste(MGName, ContigNr, GeneNR, sep = "_")
#   
#   fasta.tab.filtered = fasta.tab[which(fasta.tab$Coverage >= minCov), ]
#   fasta.tab.filtered = fasta.tab.filtered[, -4]
#   fasta.tab.filtered$Seq = fasta.tab.filtered$Sequence
#   fasta.tab.filtered = fasta.tab.filtered[, -c(1:2)]
#   
#   file = paste(MG_Name, "_FR.fna", sep = "")
#   write.table(
#     file = file,
#     fasta.tab.filtered,
#     sep = "\t",
#     row.names = F,
#     col.names = F,
#     quote = F
#   )
# }

mkdir -p $outDir/07_RenamedGenes/fna/
mkdir -p $outDir/07_RenamedGenes/faa/
mkdir -p $outDir/07_RenamedGenes/fna_fasta/
mkdir -p $outDir/07_RenamedGenes/faa_fasta/
mkdir -p $outDir/07_RenamedGenes/fna_fasta_c/
mkdir -p $outDir/07_RenamedGenes/faa_fasta_c/

mkdir -p $outDir/07_RenamedGenes/fna_complete/
mkdir -p $outDir/07_RenamedGenes/faa_complete/

mv $outDir/06_TabularFiles/fna/*_FR.fna $outDir/07_RenamedGenes/fna
mv $outDir/06_TabularFiles/faa/*_FR.faa $outDir/07_RenamedGenes/faa

mv $outDir/06_TabularFiles/fna/*_FR_complete.fna $outDir/07_RenamedGenes/fna_complete
mv $outDir/06_TabularFiles/faa/*_FR_complete.faa $outDir/07_RenamedGenes/faa_complete


( cd $outDir/07_RenamedGenes/fna && ls *.fna ) | cut -f 1 -d '.' > $outDir/Temp/Format_Files.txt 
for s in $(cat $outDir/Temp/Format_Files.txt ); do
echo "Reformatting (Tabular to fasta) "${s}"..."  
python $tab2fasta $outDir/07_RenamedGenes/fna/${s}.fna 2 1 >  $outDir/07_RenamedGenes/fna_fasta/${s}_f.fna
python $tab2fasta $outDir/07_RenamedGenes/faa/${s}.faa 2 1 >  $outDir/07_RenamedGenes/faa_fasta/${s}_f.faa
rm $outDir/07_RenamedGenes/fna/${s}.fna
rm $outDir/07_RenamedGenes/faa/${s}.faa
done

( cd $outDir/07_RenamedGenes/fna_complete && ls *.fna ) | cut -f 1 -d '.' > $outDir/Temp/Format_Files.txt 
for s in $(cat $outDir/Temp/Format_Files.txt ); do
echo "Reformatting (Tabular to fasta) "${s}"..."  
python $tab2fasta $outDir/07_RenamedGenes/fna_complete/${s}.fna 2 1 >  $outDir/07_RenamedGenes/fna_fasta_c/${s}_f.fna
python $tab2fasta $outDir/07_RenamedGenes/faa_complete/${s}.faa 2 1 >  $outDir/07_RenamedGenes/faa_fasta_c/${s}_f.faa
rm $outDir/07_RenamedGenes/fna_complete/${s}.fna
rm $outDir/07_RenamedGenes/faa_complete/${s}.faa
done

fi

if [[ "$STEP" == 4  || "$STEP" == 10 ]]; then
echo "=============================================================="
echo "Merging all Genes into one file (AllGenes_na.fasta)"
echo "=============================================================="
echo " "
mkdir -p $outDir/08_AllGenes/

cat $inDir/*fna > $outDir/08_AllGenes/AllGenes_na.fasta

echo "=============================================================="
echo "Dereplicating genes in (AllGenes_na.fasta)"
echo "=============================================================="
echo " "
usearch -fastx_uniques $outDir/08_AllGenes/AllGenes_na.fasta \
        -sizeout \
        -relabel $relabel \
        -threads $threads \
        -fastaout $outDir/08_AllGenes/AllGenes_na_derep.fasta

echo "=============================================================="
echo "Clustering dereplicated gene sequences at " $id"% identity"
echo "=============================================================="
echo " "

mkdir -p $outDir/09_nr
usearch -cluster_fast $outDir/08_AllGenes/AllGenes_na_derep.fasta \
        -id 0.$id \
        -centroids $outDir/09_nr/AllGenes_na_nr${id}.fasta
fi

if [[ "$STEP" == 5  || "$STEP" == 10 ]]; then
echo "=============================================================="
echo "Kaiju taxonomy prediction"
echo "=============================================================="
echo " "

mkdir $outDir/08_Kaiju_RefSeq/
mkdir $outDir/08_Kaiju_Progenomes/
echo "Starting Kaiju: RefSeq NR"
    kaiju -t $Kaiju_RefSeq_Nodes -f $Kaiju_RefSeq -z $threads -a greedy -e 5 -E 0.00001 -i $outDir/09_nr/AllGenes_na_nr${id}.fasta -o $outDir/08_Kaiju_RefSeq/RefSeq_Taxonomy.txt
    kaiju-addTaxonNames -t $Kaiju_RefSeq_Nodes -n $Kaiju_RefSeq_Names -i $outDir/08_Kaiju_RefSeq/RefSeq_Taxonomy.txt -o $outDir/08_Kaiju_RefSeq/RefSeq_names.txt -r superkingdom,phylum,class,order,family,genus,species -u

echo "Starting Kaiju: ProGenomes"
  kaiju -t $Kaiju_ProGenomes_Nodes -f $Kaiju_ProGenomes -z $threads -a greedy -e 5 -E 0.00001 -i $outDir/09_nr/PGC_complete_Genes_na_nr${id}.fasta -o $outDir/08_Kaiju_Progenomes/ProGenomes_Taxonomy.txt
  kaiju-addTaxonNames -t $Kaiju_ProGenomes_Nodes -n $Kaiju_ProGenomes_Names -i $outDir/08_Kaiju_Progenomes/ProGenomes_Taxonomy.txt -o $outDir/08_Kaiju_Progenomes/ProGenomes_names.txt -r superkingdom,phylum,class,order,family,genus,species -u

echo "========================================"
echo "Predicting CAZy functions of nr genes..."
echo "========================================"
echo " "

transeq -sequence $outDir/09_nr/AllGenes_na_nr${id}.fasta -outseq $outDir/09_nr/AllGenes_aa_nr${id}.fasta -frame 1
diamond blastx --more-sensitive -p $threads --id 70 -e 0.0000000001 -k 1 -d /bioinf/home/leon.dlugosch/Resources/DiamondDB_CAZy/diamond_cazy_db.dmnd -q $outDir/09_nr/AllGenes_na_nr${id}.fasta -o $outDir/09_nr/CAZy_IDs.txt

fi
