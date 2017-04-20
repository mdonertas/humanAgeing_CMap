gtex_cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
brainmicro_cormat=readRDS(file='./data/processed/humanBrainMicroarray/cormat.rds')
genes=intersect(rownames(gtex_cormat),rownames(brainmicro_cormat))
cormat=cbind(gtex_cormat[genes,],brainmicro_cormat[genes,])
rmat=cormat



library(ggplot2)
library(RColorBrewer)
pc <- prcomp(t(rmat[complete.cases(rmat),]), scale = T)

pcx=cbind(pc$x,data.frame(majtissue=c(sapply(strsplit(rownames(pc$x),'-'),function(x)x[1])[1:35],rep('Brain',ncol(brainmicro_cormat)))))
foo=pcx
pdf('./results/all_PCA.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  # geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Breast','Cells','Liver','Colon','Heart','Prostate','Artery'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 2,nudge_y = c(0,0,5,0,0,2,2,0,-2), hjust=0)+
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

rmat2=rmat[complete.cases(rmat),]
rmat2=rmat2[!rowMeans(rmat2>0)==1,]>0
pc <- prcomp(t(rmat2), scale = T)
pcx=cbind(pc$x,data.frame(majtissue=c(sapply(strsplit(rownames(pc$x),'-'),function(x)x[1])[1:35],rep('Brain',ncol(brainmicro_cormat)))))
foo=pcx
pdf('./results/all_PCA_type_of_change.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Pastel2'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
  xlab(paste("PC 1 (%", round(summary(pc)$importance[2, 1] * 100, 2), ")", sep = ""))+
  ylab(paste("PC 2 (%", round(summary(pc)$importance[2, 2] * 100, 2), ")", sep = ""))+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size=12),
        plot.title = element_text(face = 'bold',size=16))+
  # geom_text(data=subset(foo,subset = (tissue == 'Brain-Cortex')|(majtissue %in% c('Liver','Colon'))),aes(label=majtissue),fontface='bold',size=6,show.legend = FALSE,nudge_x = 0,nudge_y = -5, hjust=0)+
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
pdf('./results/all_correlations.pdf',width=20,height=15)
annotcol=data.frame(TissueType=c(sapply(strsplit(rownames(pc$x),'-'),function(x)x[1])[1:35],rep('Brain',ncol(brainmicro_cormat))))
rownames(annotcol) <- colnames(co)
TissueType=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set3'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(annotcol$TissueType))]
names(TissueType) <-unique(annotcol$TissueType)
anno_colors <- list(TissueType = TissueType)
pheatmap(co,
         scale = 'none',
         annotation_col = annotcol,
         annotation_row = annotcol,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(30),
         annotation_names_col = F,
         annotation_names_row = F,
         cutree_rows = 2,
         cutree_cols = 2,
         show_colnames = F, breaks=seq(-1,1,length.out = 31),
         display_numbers = T,
         number_format = '%.1f',number_color = 'black')
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
