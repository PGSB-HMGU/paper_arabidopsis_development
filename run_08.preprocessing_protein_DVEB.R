########################################################
### Preprocessing - Developmental set silique & seed ###
#######################################################

library(data.table)
library(stringr)


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

# MaxQuant output tables - filter and median centering

## load MaxQuant proteinGroups.txt Developmental set seed & silique
tmtprot <- fread("/path/to/proteinGroups.txt", stringsAsFactors = F, integer64 = 'double')
setkey(tmtprot, Reverse)
tmtprot <- tmtprot[.(''),]
setkey(tmtprot, `Potential contaminant`)
tmtprot <- tmtprot[.(''),]
setkey(tmtprot, `Only identified by site`)
tmtprot <- tmtprot[.(''),]

tmt_embryo <- as.matrix(tmtprot[,grep('Reporter intensity [0-9]{1,2} DV_', colnames(tmtprot)),with=F])
cnames <- colnames(tmt_embryo)
tmt_embryo[tmt_embryo==0] <- NA
tmt_embryo_log2 <- log2(tmt_embryo)
tmt_embryo_log2 <- mediancenter(tmt_embryo_log2)
tmtprot[,eval(cnames):=as.data.table(tmt_embryo_log2)]

tmtprot[,AGI_code:=gsub('\\.[0-9]*;', ';', `Protein IDs`)]
tmtprot[,AGI_code:=gsub('\\.[0-9]*$', '', AGI_code)]
tmtprot[,AGI_code:=gsub('REV__AT[0-9|C|M]G[0-9]*.[0-9];', '', tmtprot$AGI_code)]
tmtprot[,AGI_code:=gsub('REV__AT[0-9|C|M]G[0-9]*.[0-9]', '', tmtprot$AGI_code)]
tmtprot[,AGI_code:=gsub('CON__.*;', '', tmtprot$AGI_code)]
tmtprot[,AGI_code:=gsub('CON__.*', '', tmtprot$AGI_code)]
tmtprot[,AGI_code:=sapply(strsplit(tmtprot$AGI_code,";"), function(x) paste(unique(x), collapse = ";"))]
tmtprot[,'Number of genes':=str_count(tmtprot$AGI_code,pattern="AT")]
tmtprot <- tmtprot[order(tmtprot$AGI_code, -abs(tmtprot$`Razor + unique peptides`)),]
tmtprot <- tmtprot[!duplicated(tmtprot$AGI_code),]

fwrite(tmtprot, paste0("/path/to/", "preproc_proteinGroups_DV_EB.txt"), sep = '\t', col.names = T, row.names = F)


## load MaxQuant phospho (STY)sites.txt Developmental set seed & silique
tmtphos <- fread("/path/to/Phospho (STY)Sites.txt", stringsAsFactors = F, integer64 = 'double')
setkey(tmtphos, Reverse)
tmtphos <- tmtphos[.(''),]
setkey(tmtphos, `Potential contaminant`)
tmtphos <- tmtphos[.(''),]

tmtphos[,AGI_code:=gsub('\\.[0-9]*;', ';', `Proteins`)]
tmtphos[,AGI_code:=gsub('\\.[0-9]*$', '', AGI_code)]
tmtphos[,AGI_code:=gsub('REV__AT[0-9|C|M]G[0-9]*.[0-9];', '', tmtphos$AGI_code)]
tmtphos[,AGI_code:=gsub('REV__AT[0-9|C|M]G[0-9]*.[0-9]', '', tmtphos$AGI_code)]
tmtphos[,AGI_code:=gsub('CON__.*;', '', tmtphos$AGI_code)]
tmtphos[,AGI_code:=gsub('CON__.*', '', tmtphos$AGI_code)]
tmtphos[,AGI_code:=sapply(strsplit(tmtphos$AGI_code,";"), function(x) paste(unique(x), collapse = ";"))]
tmtphos[,'Number of genes':=str_count(tmtphos$AGI_code, pattern="AT")]
tmtphos$'class I'<- ifelse(tmtphos$`Localization prob`>0.75, "+", "")

tmtp_embryo <- as.matrix(tmtphos[,grep('Reporter intensity [0-9]{1,2} P_', colnames(tmtphos)),with=F])
cnames <- colnames(tmtp_embryo)
tmtp_embryo[tmtp_embryo==0] <- NA
tmtp_embryo_log2 <- log2(tmtp_embryo)
tmtphos[,eval(cnames):=as.data.table(tmtp_embryo_log2)]

fwrite(tmtphos, paste0("/path/to/", "preproc_phosphoSTYsites_DV_EB.txt"),sep = '\t', col.names = T, row.names = F)


