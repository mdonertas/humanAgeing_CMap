source('./../shared/functions/functions.R')
setwd('./../GTEx/')
filesx=list.files('./data/processed/cors/',full.names = T)
tisx=gsub('.rds','',list.files('./data/processed/cors/',full.names = F))
tsx=gsub(')','',gsub('[(]','-',tisx))

corss=lapply(filesx,readRDS)
names(corss)=tsx
genelist=unique(unname(unlist(sapply(corss,rownames))))

rmat=matrix(NA,nrow=length(genelist),ncol=length(tsx),dimnames=list(genelist,tsx))
pmat=matrix(NA,nrow=length(genelist),ncol=length(tsx),dimnames=list(genelist,tsx))
padjmat=matrix(NA,nrow=length(genelist),ncol=length(tsx),dimnames=list(genelist,tsx))

for(nm in names(corss)){
  print(nm)
  corx=corss[[nm]]
  genes=rownames(corx)
  rmat[genes,nm]=corx$rho
  pmat[genes,nm]=corx$p
  padjmat[genes,nm]=corx$padj
}

exps=lapply(paste('./data/processed/exp2/',tsx,'.rds',sep=''),readRDS)
setwd('~/GD/humanAgeing_CMap/')

names(exps)=tsx

mytisx=tsx[sapply(exps,ncol)>20]
mytisx=grep('Cells',mytisx,v=T,invert=T)
exps=lapply(mytisx,function(nm)exps[[nm]])
rmat=rmat[,mytisx]
pmat=pmat[,mytisx]
padjmat=padjmat[,mytisx]

library(ggplot2)
library(RColorBrewer)
pc <- prcomp(t(rmat[complete.cases(rmat),]), scale = T)
pcx=cbind(pc$x,data.frame(tissue=rownames(pc$x),majtissue=sapply(strsplit(rownames(pc$x),'-'),function(x)x[1])))
foo=pcx
pdf('./results/GTEx/PCA.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Breast','Cells','Liver','Colon','Heart','Prostate','Artery'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 2,nudge_y = c(0,0,5,0,0,2,2,0,-2), hjust=0)+
  guides(color=guide_legend(title='Tissue Type'))
ggplot(pcx,aes(x=PC3,y=PC4,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 3] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 4] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  geom_text(data=subset(foo,subset = (tissue == 'Brain-Hippocampus')|(majtissue %in% c('Cells','Prostate'))),aes(label=tissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 2,hjust = 0)+
  # geom_text(aes(x=PC3+8,y=PC4,label=tissue),fontface='bold',size=6,show.legend = FALSE)+
  guides(color=guide_legend(title='Tissue Type'))
dev.off()

rmat2=rmat[complete.cases(rmat),]
rmat2=rmat2[!rowMeans(rmat2>0)==1,]>0
pc <- prcomp(t(rmat2), scale = T)
pcx=cbind(pc$x,data.frame(tissue=rownames(pc$x),majtissue=sapply(strsplit(rownames(pc$x),'-'),function(x)x[1])))
foo=pcx
pdf('./results/GTEx/PCA_type_of_change.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Liver','Colon'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 0,nudge_y = -5, hjust=0)+
  guides(color=guide_legend(title='Tissue Type'))
ggplot(pcx,aes(x=PC3,y=PC4,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 3] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 4] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  # geom_text(data=subset(foo,subset = (tissue == 'Brain-Hippocampus')|(majtissue %in% c('Cells','Prostate'))),aes(label=tissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 2,hjust = 0)+
  # geom_text(aes(x=PC3+8,y=PC4,label=tissue),fontface='bold',size=6,show.legend = FALSE)+
  guides(color=guide_legend(title='Tissue Type'))
dev.off()

co=cor(rmat,use = 'pairwise',method = 's')
library(corrplot)
pdf('./results/GTEx/tissue_correlations.pdf',width=15,height=15)
corrplot(co, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 5,cl.pos = "b", tl.col = "indianred4", tl.cex = 1.5, cl.cex = 1.5, addCoef.col = "gray15", number.digits = 1, number.cex = 0.75, col = colorRampPalette(c("darkred","white","midnightblue"))(100))
dev.off()

saveRDS(rmat,file='./data/processed/GTEx/cormat.rds')
saveRDS(pmat,file='./data/processed/GTEx/pmat.rds')
saveRDS(padjmat,file='./data/processed/GTEx/padjmat.rds')
saveRDS(exps,file='./data/processed/GTEx/exps.rds')

cormat2=cormat[complete.cases(cormat),]
annotcol=data.frame(TissueType=sapply(strsplit(colnames(cormat2),'-'),function(x)x[1]))
rownames(annotcol) <- colnames(cormat2)
TissueType=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set3'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(annotcol$TissueType))]
names(TissueType) <-unique(annotcol$TissueType)
anno_colors <- list(TissueType = TissueType)
rmat2=cormat2

pdf('./results/GTEx/expression_heatmap.pdf',width = 10,height = 9)
pheatmap(rmat2,
         scale = 'none',
         annotation_col = annotcol,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(30),
         annotation_names_col = F,kmeans_k = 30,show_rownames = T,
         show_colnames = T, breaks=seq(-0.5,0.5,length.out = 31))
dev.off()

cormat2=cormat2[,grep('Brain',colnames(cormat2),invert=T)]
annotcol=data.frame(TissueType=sapply(strsplit(colnames(cormat2),'-'),function(x)x[1]))
rownames(annotcol) <- colnames(cormat2)
TissueType=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set3'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(annotcol$TissueType))]
names(TissueType) <-unique(annotcol$TissueType)
anno_colors <- list(TissueType = TissueType)
rmat2=cormat2

pdf('./results/GTEx/expression_heatmap_nobrain.pdf',width = 10,height = 9)
pheatmap(rmat2,
         scale = 'none',
         annotation_col = annotcol,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(30),
         annotation_names_col = F,kmeans_k = 30,show_rownames = T,
         show_colnames = T, breaks=seq(-0.4,0.4,length.out = 31))
dev.off()

################# 
# CMap Analysis #
#################

cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
pmat=readRDS(file='./data/processed/GTEx/pmat.rds')
padjmat=readRDS(file='./data/processed/GTEx/padjmat.rds')
exps=readRDS(file='./data/processed/GTEx/exps.rds')

library(biomaRt)
humanmart=useMart('ensembl','hsapiens_gene_ensembl')
getBM(attributes = c('ensembl_gene_id','hgnc_symbol','description'),filters = 'ensembl_gene_id',values = names(which(rowMeans(cormat>0)==1|rowMeans(cormat>0)==0)),mart = humanmart)

cormat2=cormat
cormat=cormat[complete.cases(cormat),]
dim(cormat)
# [1] 19064    35

quantile(rowSums(cormat>0),probs=c(0.005,1-0.01/2))

xx=data.frame(`Datasets`=rowSums(cormat>0))
xx$cons=c('-','Consistent')[1+(xx$Datasets>=32 | xx$Datasets <=4)]
pdf('./results/GTEx/consistency_distribution.pdf',width = 9)
ggplot(xx,aes(x=Datasets,fill=cons))+
  geom_bar(stat='count') +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(size=14)) +
  scale_fill_manual(values=c('gray35','darkred'))+
  xlab('Number of Datasets')+
  ylab('Count')+
  guides(fill=guide_legend(title='')) + 
  ggtitle('For each gene, how many datasets have a positive correlation with age?\n*Irrespective of the effect size')
dev.off()

ups=names(which(rowSums(cormat>0)>=32))
downs=names(which(rowSums(cormat>0)<=4))

ups2=as.character(read.table('./data/processed/humanBrainMicroarray_LuIncluded/ensembl_consistent_up.txt')[,1])
downs2=as.character(read.table('./data/processed/humanBrainMicroarray_LuIncluded/ensembl_consistent_down.txt')[,1])

xx=list(ups,downs,ups2,downs2)
names(xx)=c('GTEx-up','GTEx-down','Brain_micro-up','Brain_micro-down')
round(sapply(xx,function(x){
  sapply(xx,function(y){
    mean(x%in%y)
  })
}),3)*100

library('pheatmap')

Var1        <- c("navy", "darkgreen")
names(Var1) <- c("Exp1", "Exp2")
anno_colors <- list(Var1 = Var1)

annotcol=data.frame(TissueType=sapply(strsplit(colnames(cormat2),'-'),function(x)x[1]))
rownames(annotcol) <- colnames(cormat2)
TissueType=c(brewer.pal(8,'Pastel2'),brewer.pal(8,'Set2'),brewer.pal(8,'Accent'))[1:length(unique(annotcol$TissueType))]
names(TissueType) <-unique(annotcol$TissueType)
anno_colors <- list(TissueType = TissueType)
pdf("./results/GTEx/Brain_microarrayUpsandDowns_inGTEx.pdf")
pheatmap(cormat2[c(ups2,downs2),],
         scale = 'none',
         cluster_rows = F,
         annotation_col = annotcol,
         show_rownames = F,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(50),annotation_names_col = F,gaps_row = length(ups2))
dev.off()

upsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = ups,mart = humanmart)[,1])
downsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = downs,mart = humanmart)[,1])

upsp2=setdiff(upsp,downsp)
downsp2=setdiff(downsp,upsp)

write.table(upsp2,file = './data/processed/GTEx/all_up.grp',quote = F,row.names = F,col.names = F)
write.table(downsp2,file = './data/processed/GTEx/all_down.grp',quote = F,row.names = F,col.names = F)

#############################  
# do the same without brain #
#############################

cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
pmat=readRDS(file='./data/processed/GTEx/pmat.rds')
padjmat=readRDS(file='./data/processed/GTEx/padjmat.rds')
exps=readRDS(file='./data/processed/GTEx/exps.rds')

cormat=cormat[,grep('Brain',colnames(cormat),invert=T,v=T)]

library(biomaRt)
humanmart=useMart('ensembl','hsapiens_gene_ensembl')
getBM(attributes = c('ensembl_gene_id','hgnc_symbol','description'),filters = 'ensembl_gene_id',values = names(which(rowMeans(cormat>0)==1|rowMeans(cormat>0)==0)),mart = humanmart)

cormat2=cormat
cormat=cormat[complete.cases(cormat),]
dim(cormat)
# [1] 19456    22

quantile(rowSums(cormat>0),probs=c(0.005,1-0.01/2))

xx=data.frame(`Datasets`=rowSums(cormat>0))
xx$cons=c('-','Consistent')[1+(xx$Datasets>=20 | xx$Datasets <=3)]
pdf('./results/GTEx/consistency_distribution_noBrain.pdf',width = 9)
ggplot(xx,aes(x=Datasets,fill=cons))+
  geom_bar(stat='count') +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(size=14)) +
  scale_fill_manual(values=c('gray35','darkred'))+
  xlab('Number of Datasets')+
  ylab('Count')+
  guides(fill=guide_legend(title='')) + 
  ggtitle('For each gene, how many datasets have a positive correlation with age?\n*Irrespective of the effect size')
dev.off()

ups=names(which(rowSums(cormat>0)>=20))
downs=names(which(rowSums(cormat>0)<=3))

ups2=as.character(read.table('./data/processed/humanBrainMicroarray_LuIncluded/ensembl_consistent_up.txt')[,1])
downs2=as.character(read.table('./data/processed/humanBrainMicroarray_LuIncluded/ensembl_consistent_down.txt')[,1])

xx=list(ups,downs,ups2,downs2)
names(xx)=c('GTEx-up','GTEx-down','Brain_micro-up','Brain_micro-down')
round(sapply(xx,function(x){
  sapply(xx,function(y){
    mean(x%in%y)
  })
}),3)*100

upsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = ups,mart = humanmart)[,1])
downsp=unique(getBM(attributes = c('affy_hg_u133a'),filters = 'ensembl_gene_id',values = downs,mart = humanmart)[,1])

upsp2=setdiff(upsp,downsp)
downsp2=setdiff(downsp,upsp)

write.table(upsp2,file = './data/processed/GTEx/nobrain_up.grp',quote = F,row.names = F,col.names = F)
write.table(downsp2,file = './data/processed/GTEx/nobrain_down.grp',quote = F,row.names = F,col.names = F)

##################################

gtex=read.csv('data/processed/GTEx/CMap_results/all_consistent_byname.csv')
psx=as.numeric(as.character(gtex$p))
# psx[is.na(psx)]=1
gtex$padj=p.adjust(psx,method='fdr')

as.data.frame(apply(gtex[which(gtex$padj<0.05),],2,function(x)paste('|',x)))



gtex_nobrain=read.csv('data/processed/GTEx/CMap_results/nobrain_consistent_byname.csv')
psx=as.numeric(as.character(gtex_nobrain$p))
# psx[is.na(psx)]=1
gtex_nobrain$padj=p.adjust(psx,method='fdr')

as.data.frame(apply(gtex_nobrain[which(gtex_nobrain$padj<0.05),],2,function(x)paste('|',x)))
