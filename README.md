# paper_arabidopsis_development
Sources code used in the analysis of Arabidopsis development 

The shell and R scripts can be used to reproduce the processing of the RNAseq and proteomics data.

1. `run_01.fastqc.raw.sh` (quality control of raw data)
2. `run_02.trimmomatic.sh` (trimming of sequencing adapters and low quality reads)
3. `run_03.fastqc.trimmed.sh` (quality control of trimmed data)
4. `run_04.kallisto.sh` (mapping of trimmed reads to `Araport11_genes.201606.cdna.fasta.gz`, can be downloaded from [here](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FSequences%2FAraport11_blastsets) and indexed with `kallisto index [arguments] FASTA-files`)
5. `run_05.txi2gene_kallisto.R` (R script to collapse transcripts to gene level)
6. `run_06.preprocessing_rna.DV.R` (R script to normalize and filter transcript data for development datasets)
7. `run_07.preprocessing_protein_DVFL.R` (R script to normalize and filter proteinGroups.txt and phospho (STY)site.txt MaxQuant output data files for flower dataset)
8. `run_08.preprocessing_protein_DVEB.R` (R script to normalize and filter proteinGroups.txt and phospho (STY)site.txt MaxQuant output data files for seed & silique dataset)
9. `run_09.preprocessing_protein_DVLF.R` (R script to normalize and filter proteinGroups.txt and phospho (STY)site.txt MaxQuant output data files for leaf dataset)
