######################################################
### Preprocessing - Developmental set rosette leaf ###
######################################################

library(data.table)
library(stringr)
library(sva)
library(reshape2)

mediancenter <- function(x, refrows=rownames(x), islog) {
  if (islog) {
    gm <- apply(x[refrows,], 2, median, na.rm = T)
    x <- sweep(x, 2, gm, `-`)+median(x[refrows,], na.rm=T)
  } else {
    gm <- apply(x[refrows,], 2, median, na.rm = T)
    x <- sweep(x, 2, gm, `/`)*median(x[refrows,], na.rm=T)
  }
}


# MaxQuant output table - combat normalization and filter - Full proteome

## load MaxQuant proteinGroups.txt Developmental set rosette leaf
tmtprot <- fread("/path/to/proteinGroups.txt", stringsAsFactors = F, integer64 = 'double')
naming <- fread("/path/to/naming_fullproteomes.csv", stringsAsFactors = F, integer64 = 'double')
naming[,Sample:=paste0(gsub('[0-9] ', '', Name),'_',Sample)]
setkey(naming, Name)
tmt <- as.matrix(tmtprot[,grep('Reporter intensity [0-9]{1,2} L', colnames(tmtprot)),with=F])
oldnames <- colnames(tmt)
colnames(tmt) <- naming[gsub('Reporter intensity ', '', colnames(tmt)),Sample]
rownames(tmt) <- tmtprot$id
tmt[tmt==0] <- NA
bat <- as.factor(grepl('LF2', colnames(tmt)))
tmt <- do.call(cbind, lapply(unique(bat), function(i) mediancenter(tmt[,bat==i], islog = F)))

tmt_combat <- log2(tmt)
colnames(tmt_combat) <- gsub('LF[12]_', '', colnames(tmt_combat))
idx <- which(rowSums(!is.na(tmt_combat[,bat==F]))>=2 & rowSums(!is.na(tmt_combat[,bat==T]))>=2)
tmt_combat <- ComBat(dat = tmt_combat[idx,], batch = bat)
tmt_combat <- rbind(tmt_combat,log2(tmt[setdiff(rownames(tmt),rownames(tmt_combat)),]))
tmt_combat <- tmt_combat[order(as.integer(rownames(tmt_combat))),]

tmt_combat <- as.data.table(tmt_combat)
setnames(tmt_combat, oldnames)
setkey(tmtprot, id)
tmt <- as.matrix(tmtprot[,grep('Reporter intensity [0-9]{1,2} L', colnames(tmtprot)),with=F])
tmtprot[,c(grep('Reporter intensity [0-9]{1,2} L', colnames(tmtprot), value = T)):=tmt_combat[,grep('Reporter intensity [0-9]{1,2} L', colnames(tmt_combat)),with=F]]
fwrite(tmtprot, paste0("/path/to/", "ComBat_normalized.txt"), sep = '\t', col.names = T, row.names = F)


## load Combat normalized file - filter and combine Leaf set 1 & 2 Full proteome
tmtprot <- fread("/path/to/ComBat_normalized.txt", stringsAsFactors = F, integer64 = 'double')
setkey(tmtprot, Reverse)
tmtprot <- tmtprot[.(''),]
setkey(tmtprot, `Potential contaminant`)
tmtprot <- tmtprot[.(''),]
setkey(tmtprot, `Only identified by site`)
tmtprot <- tmtprot[.(''),]

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

naming <- fread("/path/to/naming_fullproteomes.csv", stringsAsFactors = F, integer64 = 'double')
naming[,Sample:=paste0(gsub('[0-9] ', '', Name),'_',Sample)]
setkey(naming, Name)
tmt_leaf <- as.matrix(tmtprot[,grep('Reporter intensity [0-9]{1,2} L', colnames(tmtprot)),with=F])
oldnames <- colnames(tmt_leaf)
colnames(tmt_leaf) <- naming[gsub('Reporter intensity ', '', colnames(tmt_leaf)),Sample]
rownames(tmt_leaf) <- tmtprot[,AGI_code]

tmp <- as.data.table(melt(data.table(tmt_leaf, keep.rownames = T), id.vars = 'rn', variable.factor = F))
tmp[,variable:=gsub('LF[0-9]{1}_', 'LF_', variable)]
tmp <- tmp[,list(value=mean(value, na.rm = T)),by=list(rn, variable)]
tmt_leaf_log2 <- acast(data = tmp, formula = rn~variable, value.var = 'value')
setkey(tmtprot, AGI_code)
tmtprot[rownames(tmt_leaf_log2),paste0('Average Reporter intensity ', colnames(tmt_leaf_log2)):=data.table(tmt_leaf_log2, check.names = F, keep.rownames = F)]
fwrite(tmtprot, paste0("/path/to/", "ComBat_normalized_average.txt"), sep = '\t', col.names = T, row.names = F)


# MaxQuant output table - combat normalization and filter - phosphorylation sites

## load MaxQuant phospho (STY)sites.txt Developmental set rosette leaf
tmtphos <- fread("/path/to/Phospho (STY)Sites.txt", stringsAsFactors = F, integer64 = 'double')
naming <- fread("/path/to/naming_fullproteomes.csv", stringsAsFactors = F, integer64 = 'double')
naming[,Sample:=paste0(gsub('[0-9] ', '', Name),'_',Sample)]
naming[,Name:=gsub(' LF', ' P_LF',Name)]
setkey(naming, Name)
tmt <- as.matrix(tmtphos[,grep('Reporter intensity [0-9]{1,2} P_L', colnames(tmtphos)),with=F])
oldnames <- colnames(tmt)
colnames(tmt) <- naming[gsub('Reporter intensity ', '', colnames(tmt)),Sample]
rownames(tmt) <- tmtphos$id
tmt[tmt==0] <- NA
bat <- as.factor(grepl('LF2', colnames(tmt)))
tmt <- do.call(cbind, lapply(unique(bat), function(i) mediancenter(tmt[,bat==i], islog = F)))

tmt_combat <- log2(tmt)
colnames(tmt_combat) <- gsub('LF[12]_', '', colnames(tmt_combat))
idx <- which(rowSums(!is.na(tmt_combat[,bat==F]))>=2 & rowSums(!is.na(tmt_combat[,bat==T]))>=2)
tmt_combat <- ComBat(dat = tmt_combat[idx,], batch = bat)
tmt_combat <- rbind(tmt_combat,log2(tmt[setdiff(rownames(tmt),rownames(tmt_combat)),]))
tmt_combat <- tmt_combat[order(as.integer(rownames(tmt_combat))),] 

tmt_combat <- as.data.table(tmt_combat)
setnames(tmt_combat, oldnames)
setkey(tmtphos, id)
tmt <- as.matrix(tmtphos[,grep('Reporter intensity [0-9]{1,2} L', colnames(tmtphos)),with=F])
tmtphos[,c(grep('Reporter intensity [0-9]{1,2} P_L', colnames(tmtphos), value = T)):=tmt_combat[,grep('Reporter intensity [0-9]{1,2} P_L', colnames(tmt_combat)),with=F]]
fwrite(tmtphos, paste0("/path/to/", "ComBat_normalized_phospho.txt"), sep = '\t', col.names = T, row.names = F)


## Combat normalized file - filter and combine Leaf set 1 & 2 phosphorylation sites
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
 
naming <- fread("/path/to/naming_fullproteomes.csv", stringsAsFactors = F, integer64 = 'double')
naming[,Sample:=paste0(gsub('[0-9] ', '', Name),'_',Sample)]
naming[,Name:=gsub(' LF', ' P_LF',Name)]
setkey(naming, Name)
tmt_leaf <- as.matrix(tmtphos[,grep('Reporter intensity [0-9]{1,2} P_L', colnames(tmtphos)),with=F])
oldnames <- colnames(tmt_leaf)
colnames(tmt_leaf) <- naming[gsub('Reporter intensity ', '', colnames(tmt_leaf)),Sample]
rownames(tmt_leaf) <- tmtphos[,id]

tmp <- as.data.table(melt(data.table(tmt_leaf, keep.rownames = T), id.vars = 'rn', variable.factor = F))
tmp[,variable:=gsub('LF[0-9]{1}_', 'LF_', variable)]
tmp <- tmp[,list(value=mean(value, na.rm = T)),by=list(rn, variable)]
tmt_leaf_log2 <- acast(data = tmp, formula = rn~variable, value.var = 'value')
setkey(tmtphos, id)
tmtphos[.(as.integer(rownames(tmt_leaf_log2))),paste0('Average Reporter intensity ', colnames(tmt_leaf_log2)):=data.table(tmt_leaf_log2, check.names = F, keep.rownames = F)]
fwrite(tmtphos, paste0("/path/to/", "ComBat_normalized_phospho_average.txt"), sep = '\t', col.names = T, row.names = F)







