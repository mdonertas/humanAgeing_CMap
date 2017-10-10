source('../shared/functions/functions.R')
library(biomaRt)
gtex_cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
brainmicro_cormat=readRDS(file='./data/processed/humanBrainMicroarray/cormat.rds')
genes=union(rownames(gtex_cormat),rownames(brainmicro_cormat))
cormat=matrix(NA,ncol=ncol(gtex_cormat)+ncol(brainmicro_cormat),nrow=length(genes),dimnames=list(genes,c(paste('GTEx',colnames(gtex_cormat),sep='_'),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep='')))))
cormat[rownames(gtex_cormat),paste('GTEx',colnames(gtex_cormat),sep='_')]=gtex_cormat
cormat[rownames(brainmicro_cormat),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep=''))]=brainmicro_cormat
brtis=grep('Brain',colnames(cormat),v=T)
cormat=cormat[,brtis]
ups=names(which(rowMeans(cormat>0)==1))
downs=names(which(rowMeans(cormat<0)==1))
humanmart=useMart('ensembl','hsapiens_gene_ensembl')
up_affy=as.character(getBM(attributes ='affy_hg_u133a' ,filters = 'ensembl_gene_id',values = ups,mart = humanmart)[,1])
down_affy=as.character(getBM(attributes ='affy_hg_u133a' ,filters = 'ensembl_gene_id',values = downs,mart = humanmart)[,1])
system('mkdir ./data/processed/allBrain')
write.table(up_affy,file='./data/processed/allBrain/up_affy.grp',quote = F,row.names = F,col.names = F)
write.table(down_affy,file='./data/processed/allBrain/down_affy.grp',quote = F,row.names = F,col.names = F)

cmap=read.csv('./data/processed/allBrain/CMap_results/allbrain.csv')
ps=as.numeric(as.character(cmap$p))
cmap$padj=p.adjust(ps,method='fdr')
head(cmap)
cmap_sig=cmap[which(cmap$padj<0.05),]
write.table(cmap_sig,file='./data/processed/allBrain/CMap_results/significantres.csv',quote = F,row.names = F,col.names = T,sep = ',')



bg=rownames(cormat[complete.cases(cormat),])
genelist=setNames(rep(0,length(bg)),bg)
genelist[ups]=1
genelist[downs]=-1
go_up=go_enrich.test(genelist = genelist,selection = function(x)x==1,nodesize = 5,ontologyx = 'BP')
go_down=go_enrich.test(genelist = genelist,selection = function(x)x==-1,nodesize = 5,ontologyx = 'BP')

save(list=ls(),file='./data/processed/allBrain/all.RData')

###############################

rm(list=ls())
source('../shared/functions/functions.R')
library(biomaRt)
gtex_cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
brainmicro_cormat=readRDS(file='./data/processed/humanBrainMicroarray/cormat.rds')
genes=union(rownames(gtex_cormat),rownames(brainmicro_cormat))
cormat=matrix(NA,ncol=ncol(gtex_cormat)+ncol(brainmicro_cormat),nrow=length(genes),dimnames=list(genes,c(paste('GTEx',colnames(gtex_cormat),sep='_'),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep='')))))
cormat[rownames(gtex_cormat),paste('GTEx',colnames(gtex_cormat),sep='_')]=gtex_cormat
cormat[rownames(brainmicro_cormat),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep=''))]=brainmicro_cormat
brtis=grep('Lu2004',grep('Brain',colnames(cormat),v=T),v=T,invert = T)
cormat=cormat[,brtis]
ups=names(which(rowMeans(cormat>0)==1))
downs=names(which(rowMeans(cormat<0)==1))
humanmart=useMart('ensembl','hsapiens_gene_ensembl')
up_affy=as.character(getBM(attributes ='affy_hg_u133a' ,filters = 'ensembl_gene_id',values = ups,mart = humanmart)[,1])
down_affy=as.character(getBM(attributes ='affy_hg_u133a' ,filters = 'ensembl_gene_id',values = downs,mart = humanmart)[,1])
write.table(up_affy,file='./data/processed/allBrain/up_affy_noLu.grp',quote = F,row.names = F,col.names = F)
write.table(down_affy,file='./data/processed/allBrain/down_affy_noLu.grp',quote = F,row.names = F,col.names = F)

cmap=read.csv('./data/processed/allBrain/CMap_results/allbrain_noLu.csv')
ps=as.numeric(as.character(cmap$p))
cmap$padj=p.adjust(ps,method='fdr')
head(cmap)
cmap_sig=cmap[which(cmap$padj<0.05),]
write.table(cmap_sig,file='./data/processed/allBrain/CMap_results/significantres_noLu.csv',quote = F,row.names = F,col.names = T,sep = ',')



bg=rownames(cormat[complete.cases(cormat),])
genelist=setNames(rep(0,length(bg)),bg)
genelist[ups]=1
genelist[downs]=-1
go_up=go_enrich.test(genelist = genelist,selection = function(x)x==1,nodesize = 5,ontologyx = 'BP')
go_down=go_enrich.test(genelist = genelist,selection = function(x)x==-1,nodesize = 5,ontologyx = 'BP')

save(list=ls(),file='./data/processed/allBrain/all_noLu.RData')

