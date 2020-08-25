###########################
### Kallisto - txi2gene ###
###########################

### collapse to gene level ###

library(tximport)
library(readr)


IDs<-read.delim("sample/abundance.tsv", stringsAsFactors=F)
IDs<-transform(IDs, gene_id=vapply(strsplit(IDs$target_id, "[.]"), '[', 1, FUN.VALUE=character(1)))

tx2gene<-IDs[,c(1,6)]
colnames(tx2gene)[1]<-"transcript_id"

Experimental_Design<-read.delim("/path/to/sample_sheet.tsv")

files.kallisto<-file.path("/path/to/Kallisto",samples,"abundance.h5")
names(files.kallisto)<-samples
all(file.exists(files.kallisto))

txi2<-tximport(files.kallisto, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

write.table(txi2$counts, "RNAseq.gene.counts.txt", sep="\t", quote=F, row.names=T)
write.table(txi2$abundance, "RNAseq.gene.abundance.txt", sep="\t", quote=F, row.names=T)


