library(corrplot)
library(scales)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(knitr)
library(reshape2)
library(ggplot2)

filesx=list.files('./data/processed/humanBrainMicroarray/exp_cors/',full.names = T)
# filesx=filesx[!grepl('Lu2004',filesx)]
filesx=grep('age.rds',filesx,v=T)
cors=lapply(filesx,readRDS)
nms=gsub('_age.rds','',sapply(strsplit(filesx,'/'),function(x)x[7]))
names(cors)=nms
genelist=unique(unname(unlist(sapply(cors,rownames))))
cormat=matrix(NA,nrow=length(genelist),ncol=length(nms),dimnames = list(genelist,nms))
pmat=matrix(NA,nrow=length(genelist),ncol=length(nms),dimnames = list(genelist,nms))
padjmat=matrix(NA,nrow=length(genelist),ncol=length(nms),dimnames = list(genelist,nms))
for(nm in nms){
  cor=cors[[nm]]
  genes=rownames(cor)
  cormat[genes,nm]=cor[,1]
  pmat[genes,nm]=cor[,2]
  padjmat[genes,nm]=cor[,3]
}
print(cormat[1:10,1:5])
print(pmat[1:10,1:5])
print(padjmat[1:10,1:5])
saveRDS(cormat,file='./data/processed/humanBrainMicroarray/cormat.rds')
saveRDS(pmat,file='./data/processed/humanBrainMicroarray/pmat.rds')
saveRDS(padjmat,file='./data/processed/humanBrainMicroarray/padjmat.rds')
rm(list=ls())


cormat=readRDS('./data/processed/humanBrainMicroarray/cormat.rds')
pmat=readRDS('./data/processed/humanBrainMicroarray/pmat.rds')
padjmat=readRDS('./data/processed/humanBrainMicroarray/padjmat.rds')

# Number of sub-datasets:
ncol(cormat)
# 26
# Number of datasources:
nms=colnames(cormat)
length(unique(sapply(strsplit(nms,'_'),function(x)x[1])))
# 7
# Number of genes detected in total:
nrow(cormat)
# 28237
# Number of genes detected in total:
nrow(cormat[complete.cases(cormat),])
# 5677
dsdist=data.frame(`Number_of_Genes`=colSums(!is.na(cormat)))
dsdist$Dataset=rownames(dsdist)
pdf('./results/humanBrainMicroarray/numberofgenes.pdf',height=5)
ggplot(dsdist,aes(x=Dataset,y=Number_of_Genes))+
  geom_bar(stat='identity')+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=90,hjust=1),
        plot.title = element_text(face = 'bold',size=16))+
  ylab('Number of genes')+ 
  xlab('')+
  ggtitle('Number of genes detected in each of the sub-datasets')
dev.off()

pdf('./results/humanBrainMicroarray/detectedgenedistribution.pdf',height=5)
genedist=melt(rowSums(!is.na(cormat)))
ggplot(genedist,aes(x=value))+
  geom_bar()+
  theme_bw()+
  ylab('Number of genes')+ 
  xlab('Number of datasets')+
  ggtitle('How many genes are detected in how many datasets?')+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))
dev.off()  

signs=data.frame(genes=rownames(cormat),
                 increase=rowSums(sign(cormat)==1,na.rm=T),
                 decrease=rowSums(sign(cormat)==-1,na.rm=T))
signs$total=signs$increase+signs$decrease
signs$consistent=factor(as.numeric(apply(signs,1,function(x)max(x[2],x[3]))),levels=26:1)
head(signs)
signs$color=signs$total==26 & as.numeric(as.character(signs$consistent))==26
pdf('./results/humanBrainMicroarray/genesharednessdistribution.pdf',height=10,width = 12)
ggplot(signs,aes(x=consistent,fill=color))+
  geom_bar(stat='count')+
  facet_wrap(~total,scales='free_x')+
  theme_bw()+
  xlab('Number of Consistent Changes')+
  ggtitle('How many consistent changes?')+
  guides(fill=F)+
  scale_fill_manual(values=c('gray30','darkred'))+
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))
dev.off()

pc <- prcomp(t(cormat[complete.cases(cormat),]), scale = T)
pcx=cbind(pc$x,data.frame(tissue=rownames(pc$x),majtissue=sapply(strsplit(rownames(pc$x),'_'),function(x)x[2])))
foo=pcx
pdf('./results/humanBrainMicroarray/PCA.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Breast','Cells','Liver','Colon','Heart','Prostate','Artery'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 2,nudge_y = c(0,0,5,0,0,2,2,0,-2), hjust=0)+
  guides(color=guide_legend(title='Tissue Type'))
ggplot(pcx,aes(x=PC3,y=PC4,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 3] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 4] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  guides(color=guide_legend(title='Tissue Type'))
dev.off()

rmat2=cormat[complete.cases(cormat),]
rmat2=rmat2[!rowMeans(rmat2>0)==1,]>0

pc <- prcomp(t(rmat2[!apply(rmat2,1,sd)==0,]), scale = T)
pcx=cbind(pc$x,data.frame(tissue=rownames(pc$x),majtissue=sapply(strsplit(rownames(pc$x),'_'),function(x)x[2])))
foo=pcx
pdf('./results/humanBrainMicroarray/PCA_type_of_change.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Liver','Colon'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 0,nudge_y = -5, hjust=0)+
  guides(color=guide_legend(title='Tissue Type'))
ggplot(pcx,aes(x=PC3,y=PC4,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 3] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 4] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  guides(color=guide_legend(title='Tissue Type'))
dev.off()

cormat2=cormat[complete.cases(cormat),]
annotcol=data.frame(TissueType=sapply(strsplit(colnames(cormat2),'_'),function(x)x[2]))
rownames(annotcol) <- colnames(cormat2)
TissueType=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))]
names(TissueType) <-unique(annotcol$TissueType)
anno_colors <- list(TissueType = TissueType)
rmat2=cormat2


pdf('./results/humanBrainMicroarray/expression_heatmap.pdf',width = 10,height = 12)
pheatmap(rmat2,
         scale = 'none',
         annotation_col = annotcol,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(30),
         annotation_names_col = F,kmeans_k = 50,show_rownames = T,
         show_colnames = T, breaks=seq(-0.5,0.5,length.out = 31))
dev.off()

co=cor(cormat,use = 'pairwise',method = 's')
pdf('./results/humanBrainMicroarray/dataset_correlations.pdf',width=9,height=7)
pheatmap(co,
         annotation_col = annotcol,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(50),
         annotation_names_col = F, annotation_names_row = F,
         show_rownames = F,
         breaks=seq(-1,1,length.out = 51),
         display_numbers = T,
         number_format = '%.1f',
         number_color = 'black',
         cutree_rows = 3,
         cutree_cols = 3,annotation_legend = T)
dev.off()

## Consistent Changes

cons_up=rownames(cormat)[which(rowMeans(cormat>0)==1)]
cons_down=rownames(cormat)[which(rowMeans(cormat<0)==1)]
bg=names(which(rowMeans(is.na(cormat))==0))
data.frame(Up=length(cons_up),Down=length(cons_down),Background=length(bg))
# |Up |Down|Background|
# |---|----|----------|
# |100|117 |5677      |


humanmart=useMart('ensembl','hsapiens_gene_ensembl')
upsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = cons_up,mart = humanmart)[,1])
downsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = cons_down,mart = humanmart)[,1])

upsp2=setdiff(upsp,downsp)
downsp2=setdiff(downsp,upsp)

write.table(upsp2,file = './data/processed/humanBrainMicroarray/all_up.grp',quote = F,row.names = F,col.names = F)
write.table(downsp2,file = './data/processed/humanBrainMicroarray/all_down.grp',quote = F,row.names = F,col.names = F)

saveRDS(cons_up,'./data/processed/humanBrainMicroarray/all_up_ensembl.rds')
saveRDS(cons_down,'./data/processed/humanBrainMicroarray/all_down_ensembl.rds')
saveRDS(upsp2,'./data/processed/humanBrainMicroarray/all_up_affy.rds')
saveRDS(downsp2,'./data/processed/humanBrainMicroarray/all_down_affy.rds')

upshugo=unique(getBM(attributes = c('hgnc_symbol'),filters = 'ensembl_gene_id',values = cons_up,mart = humanmart)[,1])
downhugo=unique(getBM(attributes = c('hgnc_symbol'),filters = 'ensembl_gene_id',values = cons_down,mart = humanmart)[,1])

saveRDS(upshugo,'./data/processed/humanBrainMicroarray/all_up_hgnc.rds')
saveRDS(downhugo,'./data/processed/humanBrainMicroarray/all_down_hgnc.rds')

length(upsp2)
length(downsp2)

res=read.csv('data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
psx=as.numeric(as.character(res$p))
# psx[is.na(psx)]=1
res$padj=p.adjust(psx,method='fdr')

as.data.frame(apply(res[which(res$padj<0.05),],2,function(x)paste('|',x)))

#############
# Except Lu #
#############

library(corrplot)
library(scales)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(knitr)
library(reshape2)
library(ggplot2)

filesx=list.files('./data/processed/humanBrainMicroarray/exp_cors/',full.names = T)
filesx=filesx[!grepl('Lu2004',filesx)]
filesx=grep('age.rds',filesx,v=T)
cors=lapply(filesx,readRDS)
nms=gsub('_age.rds','',sapply(strsplit(filesx,'/'),function(x)x[7]))
names(cors)=nms
genelist=unique(unname(unlist(sapply(cors,rownames))))
cormat=matrix(NA,nrow=length(genelist),ncol=length(nms),dimnames = list(genelist,nms))
pmat=matrix(NA,nrow=length(genelist),ncol=length(nms),dimnames = list(genelist,nms))
padjmat=matrix(NA,nrow=length(genelist),ncol=length(nms),dimnames = list(genelist,nms))
for(nm in nms){
  cor=cors[[nm]]
  genes=rownames(cor)
  cormat[genes,nm]=cor[,1]
  pmat[genes,nm]=cor[,2]
  padjmat[genes,nm]=cor[,3]
}
print(cormat[1:10,1:5])
print(pmat[1:10,1:5])
print(padjmat[1:10,1:5])
saveRDS(cormat,file='./data/processed/humanBrainMicroarray/cormat_noLu.rds')
saveRDS(pmat,file='./data/processed/humanBrainMicroarray/pmat_noLu.rds')
saveRDS(padjmat,file='./data/processed/humanBrainMicroarray/padjmat_noLu.rds')
rm(list=ls())


cormat=readRDS('./data/processed/humanBrainMicroarray/cormat_noLu.rds')
pmat=readRDS('./data/processed/humanBrainMicroarray/pmat_noLu.rds')
padjmat=readRDS('./data/processed/humanBrainMicroarray/padjmat_noLu.rds')

# Number of sub-datasets:
ncol(cormat)
# 25
# Number of datasources:
nms=colnames(cormat)
length(unique(sapply(strsplit(nms,'_'),function(x)x[1])))
# 6
# Number of genes detected in total:
nrow(cormat)
# 28200
# Number of genes detected in total:
nrow(cormat[complete.cases(cormat),])
# 10893
dsdist=data.frame(`Number_of_Genes`=colSums(!is.na(cormat)))
dsdist$Dataset=rownames(dsdist)
pdf('./results/humanBrainMicroarray/numberofgenes_noLu.pdf',height=5)
ggplot(dsdist,aes(x=Dataset,y=Number_of_Genes))+
  geom_bar(stat='identity')+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=90,hjust=1),
        plot.title = element_text(face = 'bold',size=16))+
  ylab('Number of genes')+ 
  xlab('')+
  ggtitle('Number of genes detected in each of the sub-datasets')
dev.off()

pdf('./results/humanBrainMicroarray/detectedgenedistribution_noLu.pdf',height=5)
genedist=melt(rowSums(!is.na(cormat)))
ggplot(genedist,aes(x=value))+
  geom_bar()+
  theme_bw()+
  ylab('Number of genes')+ 
  xlab('Number of datasets')+
  ggtitle('How many genes are detected in how many datasets?')+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))
dev.off()  

signs=data.frame(genes=rownames(cormat),
                 increase=rowSums(sign(cormat)==1,na.rm=T),
                 decrease=rowSums(sign(cormat)==-1,na.rm=T))
signs$total=signs$increase+signs$decrease
signs$consistent=factor(as.numeric(apply(signs,1,function(x)max(x[2],x[3]))),levels=25:1)
head(signs)
signs$color=signs$total==25 & as.numeric(as.character(signs$consistent))==25
pdf('./results/humanBrainMicroarray/genesharednessdistribution_noLu.pdf',height=10,width = 12)
ggplot(signs,aes(x=consistent,fill=color))+
  geom_bar(stat='count')+
  facet_wrap(~total,scales='free_x')+
  theme_bw()+
  xlab('Number of Consistent Changes')+
  ggtitle('How many consistent changes?')+
  guides(fill=F)+
  scale_fill_manual(values=c('gray30','darkred'))+
  theme(axis.text.x = element_text(angle=90),
        axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))
dev.off()

pc <- prcomp(t(cormat[complete.cases(cormat),]), scale = T)
pcx=cbind(pc$x,data.frame(tissue=rownames(pc$x),majtissue=sapply(strsplit(rownames(pc$x),'_'),function(x)x[2])))
foo=pcx
pdf('./results/humanBrainMicroarray/PCA_noLu.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Breast','Cells','Liver','Colon','Heart','Prostate','Artery'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 2,nudge_y = c(0,0,5,0,0,2,2,0,-2), hjust=0)+
  guides(color=guide_legend(title='Tissue Type'))
ggplot(pcx,aes(x=PC3,y=PC4,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 3] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 4] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  guides(color=guide_legend(title='Tissue Type'))
dev.off()

rmat2=cormat[complete.cases(cormat),]
rmat2=rmat2[!rowMeans(rmat2>0)==1,]>0

pc <- prcomp(t(rmat2[!apply(rmat2,1,sd)==0,]), scale = T)
pcx=cbind(pc$x,data.frame(tissue=rownames(pc$x),majtissue=sapply(strsplit(rownames(pc$x),'_'),function(x)x[2])))
foo=pcx
pdf('./results/humanBrainMicroarray/PCA_type_of_change_noLu.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Liver','Colon'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 0,nudge_y = -5, hjust=0)+
  guides(color=guide_legend(title='Tissue Type'))
ggplot(pcx,aes(x=PC3,y=PC4,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 3] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 4] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  guides(color=guide_legend(title='Tissue Type'))
dev.off()

cormat2=cormat[complete.cases(cormat),]
annotcol=data.frame(TissueType=sapply(strsplit(colnames(cormat2),'_'),function(x)x[2]))
rownames(annotcol) <- colnames(cormat2)
TissueType=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))]
names(TissueType) <-unique(annotcol$TissueType)
anno_colors <- list(TissueType = TissueType)
rmat2=cormat2


pdf('./results/humanBrainMicroarray/expression_heatmap_noLu.pdf',width = 10,height = 12)
pheatmap(rmat2,
         scale = 'none',
         annotation_col = annotcol,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(30),
         annotation_names_col = F,kmeans_k = 50,show_rownames = T,
         show_colnames = T, breaks=seq(-0.5,0.5,length.out = 31))
dev.off()

co=cor(cormat,use = 'pairwise',method = 's')
pdf('./results/humanBrainMicroarray/dataset_correlations_noLu.pdf',width=9,height=7)
pheatmap(co,
         annotation_col = annotcol,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(50),
         annotation_names_col = F, annotation_names_row = F,
         show_rownames = F,
         breaks=seq(-1,1,length.out = 51),
         display_numbers = T,
         number_format = '%.1f',
         number_color = 'black',
         cutree_rows = 3,
         cutree_cols = 3,annotation_legend = T)
dev.off()

## Consistent Changes

cons_up=rownames(cormat)[which(rowMeans(cormat>0)==1)]
cons_down=rownames(cormat)[which(rowMeans(cormat<0)==1)]
bg=names(which(rowMeans(is.na(cormat))==0))
data.frame(Up=length(cons_up),Down=length(cons_down),Background=length(bg))
# |Up |Down|Background|
# |---|----|----------|
# |192|202 |10893      |


humanmart=useMart('ensembl','hsapiens_gene_ensembl')
upsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = cons_up,mart = humanmart)[,1])
downsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = cons_down,mart = humanmart)[,1])

upsp2=setdiff(upsp,downsp)
downsp2=setdiff(downsp,upsp)

write.table(upsp2,file = './data/processed/humanBrainMicroarray/noLu_up.grp',quote = F,row.names = F,col.names = F)
write.table(downsp2,file = './data/processed/humanBrainMicroarray/noLu_down.grp',quote = F,row.names = F,col.names = F)

saveRDS(cons_up,'./data/processed/humanBrainMicroarray/noLu_up_ensembl.rds')
saveRDS(cons_down,'./data/processed/humanBrainMicroarray/noLu_down_ensembl.rds')
saveRDS(upsp2,'./data/processed/humanBrainMicroarray/noLu_up_affy.rds')
saveRDS(downsp2,'./data/processed/humanBrainMicroarray/noLu_down_affy.rds')

upshugo=unique(getBM(attributes = c('hgnc_symbol'),filters = 'ensembl_gene_id',values = cons_up,mart = humanmart)[,1])
downhugo=unique(getBM(attributes = c('hgnc_symbol'),filters = 'ensembl_gene_id',values = cons_down,mart = humanmart)[,1])

saveRDS(upshugo,'./data/processed/humanBrainMicroarray/noLu_up_hgnc.rds')
saveRDS(downhugo,'./data/processed/humanBrainMicroarray/noLu_down_hgnc.rds')

length(upsp2)
length(downsp2)

res=read.csv('data/processed/humanBrainMicroarray/CMap_results/noLu_consistent_byname.csv')
psx=as.numeric(as.character(res$p))
# psx[is.na(psx)]=1
res$padj=p.adjust(psx,method='fdr')

as.data.frame(apply(res[which(res$padj<0.05),],2,function(x)paste('|',x)))
