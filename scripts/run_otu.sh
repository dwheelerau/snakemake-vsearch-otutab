#!/bin/bash

# This is an example of a pipeline
# will perform initial paired-end
# read merging, quality filtering, chimera removal and OTU clustering.
# based on the pipeline here

# heavily inspired by https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline

# RUN PARAMs
THREADS=24
PERL=$(which perl)
VSEARCH=$(which vsearch)

RAW_READS="./projectData"
TMP="./tmp"
REF_DB="/home/dwheeler/databases/gg_13_8_otus/rep_set/97_otus.fasta"
LOGDIR="./logs"
mkdir -p $TMP
mkdir -p $LOGDIR

# MERGE PARAMs
MINOVLEN="--fastq_minovlen 60" #200
MAXDIFF="--fastq_maxdiffs 3" #15

# FILTER PARAMs
MAXEE="--fastq_maxee 0.5"
MINLEN="--fastq_minlen 400"
MAXLEN="--fastq_maxlen 500"
MAXNS="--fastq_maxns 0"

# CLUSTER PARAMs
CLUSTERMODE="--cluster_size"
CLUSTERID="--id 0.97"

# cleanup otherwise these will be incorprated due to >>

rm -f $TMP/*
rm -f $LOGDIR/*
rm -f *fastq
rm -f *.uc
rm -f *.fasta
rm -f *.stats
rm -f *.txt
rm -f all*

# directory for read stats
mkdir -p readStats
rm -f readStats/*

date

# not used (I think you need a licence)
#echo Obtaining Gold reference database for chimera detection
#
#REF_DB="./refdb/gold.fasta"
#mkdir -p ${REF_DB%/*} # removes /filename to catch just the base directory
#if [ ! -e $REF_DB ]; then
#
#    if [ ! -e Silva.gold.bacteria.zip ]; then
#        wget https://www.mothur.org/w/images/f/f1/Silva.gold.bacteria.zip
#    fi
#
#    echo Decompressing and reformatting...
#    unzip -p Silva.gold.bacteria.zip silva.gold.align | \
#        sed -e "s/[.-]//g" > $REF_DB
#
#fi
#
## extract reads if required and copy to tmp for processing
#echo
#
#echo "Moving reads to tmp and unzipping if required"

for f in $RAW_READS/*
do
  echo "Processing $f"
  if [[ "$f" == *.gz ]] ; then
    gunzip -c $f > $TMP/$(basename ${f%.gz})
  else
    cp $f $TMP
  fi
done

echo
echo Checking FASTQ format version for one file

$VSEARCH --threads $THREADS \
  --fastq_chars $(ls -1 $TMP/*.fastq | head -1)

# Process samples

for f in $TMP/*_R1_*.fastq; do
  fbname=$(basename $f)
  r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
  s=$(cut -d_ -f1 <<< "$fbname")

  echo
  echo ====================================
  echo Processing sample $s
  echo ====================================

  $VSEARCH --threads $THREADS \
      --fastq_mergepairs $f \
      --reverse $r \
      $MINOVLEN \
      $MAXDIFF \
      --fastqout $s.merged.fastq \
      --fastq_eeout

  # Commands to demultiplex and remove tags and primers
  # using e.g. cutadapt may be added here.

  echo
  echo Calculate quality statistics

  $VSEARCH --threads $THREADS \
      --fastq_eestats $s.merged.fastq \
      --output $s.stats
  echo
  echo Quality filtering

  $VSEARCH --threads $THREADS \
      --fastq_filter $s.merged.fastq \
      $MAXEE \
      $MINLEN \
      $MAXLEN \
      $MAXNS \
      --fastaout $s.filtered.fasta \
      --relabel $s. \
      --log $s.filter.log \
      --fasta_width 0

  echo
  # this part is to get ready to create OTUS
  echo Dereplicate at sample level and relabel with sample_n

  $VSEARCH --threads $THREADS \
      --derep_fulllength $s.filtered.fasta \
      --strand plus \
      --output $s.derep.fasta \
      --sizeout \
      --uc $s.derep.uc \
      --relabel $s. \
      --fasta_width 0

  done

echo Sum of unique sequences in each sample: $(cat *.derep.fasta | grep -c "^>")

# At this point there should be one fasta file for each sample
# It should be quality filtered and dereplicated.

echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

cat *.derep.fasta > all.fasta
# these will be used for mapping
cat *.filtered.fasta > sample_seqs_filtered.fasta

# move sample files to tmp
mv ./*.merged.fastq $TMP

# move files that are not needed
mv *.stats ./readStats
mv *.log $LOGDIR

# excludes all* sample_seqs_filtered that is need for mapping
shopt -s extglob
mv !(all*|sample_seqs_filtered).fasta $TMP
mv !(all*).uc $TMP
shopt -u extglob

echo
echo Dereplicate across samples and remove singletons

$VSEARCH --threads $THREADS \
  --derep_fulllength all.fasta \
  --minuniquesize 2 \
  --sizein \
  --sizeout \
  --fasta_width 0 \
  --uc all.derep.uc \
  --output all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" all.derep.fasta)

echo
echo Precluster at 98% before chimera detection

$VSEARCH --threads $THREADS \
  --cluster_size all.derep.fasta \
  --id 0.98 \
  --strand plus \
  --sizein \
  --sizeout \
  --fasta_width 0 \
  --uc all.preclustered.uc \
  --centroids all.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" all.preclustered.fasta)

echo
echo De novo chimera detection

$VSEARCH --threads $THREADS \
  --uchime_denovo all.preclustered.fasta \
  --sizein \
  --sizeout \
  --fasta_width 0 \
  --nonchimeras all.denovo.nonchimeras.fasta \

echo Unique sequences after de novo chimera detection: $(grep -c "^>" all.denovo.nonchimeras.fasta)

echo
echo Reference chimera detection

$VSEARCH --threads $THREADS \
  --uchime_ref all.denovo.nonchimeras.fasta \
  --db $REF_DB \
  --sizein \
  --sizeout \
  --fasta_width 0 \
  --nonchimeras all.ref.nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" all.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

$PERL ./map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

$PERL ./map.pl all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

echo
echo Cluster at 97% and relabel with OTU_n, generate OTU table

# removed --sizeout becuase I think this crashes fastree, I checked and it does
#removed --otutabout all.otutab.txt
# not effect the OTU table, just removes size=x form the fasta OUT file.
$VSEARCH --threads $THREADS \
  $CLUSTERMODE all.nonchimeras.fasta \
  $CLUSTERID \
  --strand plus \
  --sizein \
  --fasta_width 0 \
  --uc all.clustered.uc \
  --relabel OTU_ \
  --centroids all.otus.fasta \
  --otutabout all.otutab.txt

echo
echo Number of OTUs: $(grep -c "^>" all.otus.fasta)

echo
echo Mapping sample_seqs_filtered.fasta to all.otus.fasta

# note the use of ALL reads for this step
vsearch -usearch_global ./sample_seqs_filtered.fasta -db ./all.otus.fasta \
  $CLUSTERID -strand plus -uc readmap.uc --notrunclabels \
  --otutabout sample.otu_table.txt --log map_reads_to_clusters.log \

echo
echo Done

date
