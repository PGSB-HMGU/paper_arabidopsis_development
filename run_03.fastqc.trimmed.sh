#!/bin/bash
#$ -q queue
#$ -o /path/to/FastQC_trimmed/err_out
#$ -e /path/to/FastQC_trimmed/err_out
#$ -pe serial 4
#$ -l job_mem=2G


INDIR=/path/to/Trimmomatic

OUTDIR=/path/to/FastQC_trimmed

IN="$INDIR/sample_R%s_paired.fastq.gz"

raw=( $(printf $IN "1")  $(printf $IN "2")  )


fastqc --threads 4 --outdir $OUTDIR -f fastq ${raw[*]}


