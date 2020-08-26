###########################
### Preprocessing - RNA ###
###########################

library(data.table)
library(na.tools)

mediancenter <- function(x, refrows=NULL) {
  if (!is.null(refrows)) {
    gm <- apply(x[refrows,], 2, median, na.rm = T)
    x <- sweep(x, 2, gm, `-`)+median(x[refrows,], na.rm=T)
  } else {
    gm <- apply(x, 2, median, na.rm=T)
    x <- sweep(x, 2, gm, `-`)+median(x, na.rm=T)
  }
  return(x)
}

# Kallisto output file - normalization and filter

## load RNA dataset
rna <- fread("/path/to/RNAseq.gene.abundance.txt", stringsAsFactors = F, integer64 = 'double')

tmpmat <- as.matrix(rna[,grep('AGI', colnames(rna), invert = T),with=F])
cnames <- colnames(tmpmat)
tmpmat[tmpmat==0] <- NA
tmpmat <- log2(tmpmat)
tmpmat <- mediancenter(tmpmat)
tmpmat[tmpmat<=0] <- NA

rna[,eval(cnames):=as.data.table(tmpmat)]
rna_DV <- rna[,c(1,29:57)]
rna_DV[,RowSum:=rowSums(cbind(rna_DV[,c(2:30)]),na.rm=T)]
rna_DV <- rna_DV[rna_DV$RowSum > 0]

fwrite(rna_DV, paste0("/path/to/", "preproc_rna_DV.txt"), sep = '\t', col.names = T, row.names = F)