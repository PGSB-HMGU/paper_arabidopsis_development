#!/bin/bash
#$ -q queue
#$ -o /path/to/FastQC_raw/err_out
#$ -e /path/to/FastQC_raw/err_out
#$ -pe serial 4
#$ -l job_mem=2G


INDIR=/path/to/RNAseq

OUTDIR=/path/to/FastQC_raw

IN="$INDIR/sample_R%s.fastq.gz"

raw=( $(printf $IN "1")  $(printf $IN "2")  )


fastqc  --threads 4 --outdir $OUTDIR -f fastq ${raw[*]}


