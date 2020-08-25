#!/bin/bash
#$ -q queue
#$ -o /path/to/Kallisto/err_out
#$ -e /path/to/Kallisto/err_out
#$ -pe serial 6
#$ -l job_mem=4G


OUTDIR=/path/to/Kallisto

INPATH=/path/to/Trimmomatic

INDEX=/path/to/Kallisto/index/Araport11_genes.201606.cdna.kallisto.index

R1=sample_R1_paired.fastq.gz
R2=sample_R2_paired.fastq.gz

name=sample

mkdir $OUTDIR/$name


kallisto quant --index $INDEX --output $OUTDIR/$name --threads=6 --rf-stranded -b 100 --bias --seed 27042012 $R1 $R2 
 



