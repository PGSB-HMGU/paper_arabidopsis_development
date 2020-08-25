#!/bin/bash
#$ -q queue
#$ -o /path/to/Trimmomatic/err_out
#$ -e /path/to/Trimmomatic/err_out
#$ -pe serial 6
#$ -l job_mem=4G


INPATH=/path/to/RNAseq

R1=${INPATH}/sample_R1.fastq.gz
R2=${INPATH}/sample_R2.fastq.gz

OUTDIR=/path/to/Trimmomatic

LOG="./sample.trimmomatic.log"

r1=$(basename $R1)
r1=$OUTDIR/${r1/.fastq.gz}

r2=$(basename $R2)
r2=$OUTDIR/${r2/.fastq.gz}


java -jar  trimmomatic-0.35.jar PE -threads 6 -validatePairs $R1 $R2 ${r1}_paired.fastq.gz ${r1}_unpaired.fastq.gz ${r2}_paired.fastq.gz ${r2}_unpaired.fastq.gz ILLUMINACLIP:Illumina-PE.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 > $LOG 2>&1

