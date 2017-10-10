source('../shared/functions/functions.R')
library(biomaRt)
gtex_cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
brainmicro_cormat=readRDS(file='./data/processed/humanBrainMicroarray/cormat.rds')
genes=union(rownames(gtex_cormat),rownames(brainmicro_cormat))
cormat=matrix(NA,ncol=ncol(gtex_cormat)+ncol(brainmicro_cormat),nrow=length(genes),dimnames=list(genes,c(paste('GTEx',colnames(gtex_cormat),sep='_'),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep='')))))
cormat[rownames(gtex_cormat),paste('GTEx',colnames(gtex_cormat),sep='_')]=gtex_cormat
cormat[rownames(brainmicro_cormat),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep=''))]=brainmicro_cormat
brtis=grep('GTEx',grep('Brain',colnames(cormat),v=T),v=T)
cormat=cormat[,brtis]
ups=sort(rowMeans(cormat)[names(which(rowMeans(cormat>0)==1))],dec=T)
downs=sort(rowMeans(cormat)[names(which(rowMeans(cormat<0)==1))],dec=F)
humanmart=useMart('ensembl','hsapiens_gene_ensembl')
up_affy=getBM(attributes =c('affy_hg_u133a','ensembl_gene_id') ,filters = 'ensembl_gene_id',values = names(ups),mart = humanmart)
down_affy=getBM(attributes =c('affy_hg_u133a','ensembl_gene_id') ,filters = 'ensembl_gene_id',values = names(downs),mart = humanmart)
up_affy=up_affy[!rowSums(up_affy=='')>0,]
up_affy=setNames(up_affy$ensembl_gene_id,up_affy$affy_hg_u133a)
down_affy=down_affy[!rowSums(down_affy=='')>0,]
down_affy=setNames(down_affy$ensembl_gene_id,down_affy$affy_hg_u133a)
library(reshape2)
up_affy2=sort(unlist(lapply(unique(names(up_affy)),function(x){
  a=unname(ups[names(ups)%in%up_affy[names(up_affy)%in%x]])
  setNames(a,rep(x,length(a)))
})),dec=T)
down_affy2=sort(unlist(lapply(unique(names(down_affy)),function(x){
  a=unname(downs[names(downs)%in%down_affy[names(down_affy)%in%x]])
  setNames(a,rep(x,length(a)))
})),dec=F)
int=intersect(names(up_affy2),names(down_affy2))
up_affy3=up_affy2[!names(up_affy2)%in%int]
down_affy3=down_affy2[!names(down_affy2)%in%int]

# system('mkdir ./data/processed/brainGTEx')

up_affy4=unique(names(up_affy3))[1:500]
down_affy4=unique(names(down_affy3))[1:500]

write.table(up_affy4,file='./data/processed/brainGTEx/up_affy.grp',quote = F,row.names = F,col.names = F)
write.table(down_affy4,file='./data/processed/brainGTEx/down_affy.grp',quote = F,row.names = F,col.names = F)

cmap=read.csv('./data/processed/brainGTEx/CMap_results/GTEx_brain.csv')
ps=as.numeric(as.character(cmap$p))
cmap$padj=p.adjust(ps,method='fdr')
head(cmap)
cmap_sig=cmap[which(cmap$padj<0.05),]
write.table(cmap_sig,file='./data/processed/brainGTEx/CMap_results/significantres.csv',quote = F,row.names = F,col.names = T,sep = ',')

bg=rownames(cormat[complete.cases(cormat),])
genelist=setNames(rep(0,length(bg)),bg)
genelist[names(ups)]=1
genelist[names(downs)]=-1
go_up=go_enrich.test(genelist = genelist,selection = function(x)x==1,nodesize = 5,ontologyx = 'BP')
go_down=go_enrich.test(genelist = genelist,selection = function(x)x==-1,nodesize = 5,ontologyx = 'BP')

go_up[go_up$p.adjusted<0.05,2]
go_down[go_down$p.adjusted<0.05,2]

save(list=ls(),file='./data/processed/brainGTEx/all.RData')

###################
# correlation matrix: cormat
# genes consistently up- or down-regulated: names(ups), names(downs)
# affy ids used to query CMap: up_affy4, down_affy4

library(VennDiagram)
?venn.diagram
venn.diagram(list(microarray_down=readRDS('./data/processed/humanBrainMicroarray/all_down_affy.rds'),
                  gtex_down=down_affy4,
                  microarray_up=readRDS('./data/processed/humanBrainMicroarray/all_up_affy.rds'),
                  gtex_up=up_affy4),filename='~/Desktop/gtex_vs_micro.png',imagetype = 'png')

venn.diagram(list(microarray_down=readRDS('./data/processed/humanBrainMicroarray/all_down_ensembl.rds'),
                  gtex_down=names(downs),
                  microarray_up=readRDS('./data/processed/humanBrainMicroarray/all_up_ensembl.rds'),
                  gtex_up=names(ups)),filename='~/Desktop/gtex_vs_micro_genes.png',imagetype = 'png')

venn.diagram(list(microarray_down=readRDS('./data/processed/humanBrainMicroarray/all_down_ensembl.rds'),
                  gtex_down=unname(down_affy[names(down_affy)%in%down_affy4]),
                  microarray_up=readRDS('./data/processed/humanBrainMicroarray/all_up_ensembl.rds'),
                  gtex_up=unname(up_affy[names(up_affy)%in%up_affy4])),filename='~/Desktop/gtex_vs_micro_genes2.png',imagetype = 'png')


down_affy4

#################################################

load('./data/processed/brainGTEx/all.RData')
tisx=gsub('GTEx_','',colnames(cormat))
rm(list=setdiff(ls(),'tisx'))
permf=list.files('../GTEx/data/processed/permutations/permres/',full.names = T)
resmat=matrix(NA,ncol=2,nrow=0)
for(permfx in permf){
  load(permfx)
  cormat=cormat[,tisx]
  rm(list=setdiff(ls(),c('tisx','permf','cormat','resmat')))
  resmat=rbind(resmat,c(sum(rowMeans(cormat>0)==1,na.rm=T),sum(rowMeans(cormat<0)==1,na.rm=T)))
  rm(cormat)
}

colnames(resmat)=c('up','down')

null_sharedness=resmat
rm(list=setdiff(ls(),'null_sharedness'))
load('./data/processed/brainGTEx/all.RData')

hist(null_sharedness[,1])
mean(null_sharedness[,1]>=length(names(ups)))
mean(null_sharedness[,2]>=length(names(downs)))

null_sharedness=data.frame(Up=null_sharedness[,1],Down=null_sharedness[,2])
null_sharedness$UpCol=null_sharedness$Up>=length(names(ups))
null_sharedness$DownCol=null_sharedness$Down>=length(names(downs))

library(ggplot2)
p1=ggplot(null_sharedness,aes(x=Up))+
  geom_histogram(color='gray50',bins=200)+
  geom_segment(arrow = arrow(length = unit(0.1,'inches')),x=length(names(ups)),xend=length(names(ups)),y=15,yend=3,color='darkred',lwd=1)+
  geom_text(x=length(names(ups)),y=17,label=paste('p=',mean(null_sharedness$Up>=length(names(ups))),sep=''))+ xlab('')+
  theme_bw()+ggtitle('Up-Regulated Genes')
p2=ggplot(null_sharedness,aes(x=Down))+
  geom_histogram(color='gray50',bins=200)+
  geom_segment(arrow = arrow(length = unit(0.1,'inches')),x=length(names(downs)),xend=length(names(downs)),y=25,yend=3,color='darkred',lwd=1)+
  geom_text(x=length(names(downs)),y=28,label=paste('p=',mean(null_sharedness$Down>=length(names(downs))),sep=''))+
  theme_bw()+xlab('')+ggtitle('Down-Regulated Genes')
# system('mkdir ./results/brainGTEx')
pdf('./results/brainGTEx/no_shared_genes.pdf')
multiplot(p1,p2)
dev.off()

save(list=ls(),file='./data/processed/brainGTEx/all.RData')

