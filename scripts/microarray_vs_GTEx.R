gtex_cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
brainmicro_cormat=readRDS(file='./data/processed/humanBrainMicroarray/cormat.rds')
genes=union(rownames(gtex_cormat),rownames(brainmicro_cormat))
cormat=matrix(NA,ncol=ncol(gtex_cormat)+ncol(brainmicro_cormat),nrow=length(genes),dimnames=list(genes,c(paste('GTEx',colnames(gtex_cormat),sep='_'),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep='')))))
cormat[rownames(gtex_cormat),paste('GTEx',colnames(gtex_cormat),sep='_')]=gtex_cormat
cormat[rownames(brainmicro_cormat),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep=''))]=brainmicro_cormat

rmat=cormat

dset=sapply(strsplit(colnames(cormat),'_'),function(x)x[1])
majtis=sapply(strsplit(sapply(strsplit(colnames(cormat),'_'),function(x)x[2]),'-'),function(x)x[1])

library(scales)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
pc <- prcomp(t(rmat[complete.cases(rmat),]), scale = T)

pcx=cbind(pc$x,data.frame(majtissue=majtis))
foo=pcx
pdf('./results/all_PCA.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
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
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
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
pcx=cbind(pc$x,data.frame(majtissue=majtis))
foo=pcx
pdf('./results/all_PCA_type_of_change.pdf',width = 10,height = 9)
ggplot(pcx,aes(x=PC1,y=PC2,color=majtissue))+
  geom_point(size=5,alpha=0.8)+
  theme_bw()+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
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
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))])+
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
annotcol=data.frame(TissueType=majtis,Dataset=dset)
rownames(annotcol) <- colnames(co)
TissueType=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(pcx$majtissue))]
names(TissueType) <-unique(annotcol$TissueType)
Dataset=brewer.pal(8,'Pastel2')
names(Dataset) <-unique(annotcol$Dataset)
anno_colors <- list(TissueType = TissueType,Dataset=Dataset)
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

###########################

rm(list=ls())

br_up=readRDS('./data/processed/humanBrainMicroarray/all_up_ensembl.rds')
br_down=readRDS('./data/processed/humanBrainMicroarray/all_down_ensembl.rds')
gtex_up=readRDS('./data/processed/GTEx/all_up.rds')
gtex_down=readRDS('./data/processed/GTEx/all_down.rds')

genelist=list(`Brain Microarray\nUp`=br_up,`Brain Microarray\nDown`=br_down,`GTEx Up`=gtex_up,`GTEx Down`=gtex_down)
library(VennDiagram)
library(RColorBrewer)
venn.diagram(genelist,
             filename = './results/numberofgenes.png',
             imagetype = 'png',fill=(brewer.pal(4,'Set1')),cat.pos=c(-8,+8,0,0))

gtex_cormat=readRDS(file='./data/processed/GTEx/cormat.rds')
brainmicro_cormat=readRDS(file='./data/processed/humanBrainMicroarray/cormat.rds')
genes=union(rownames(gtex_cormat),rownames(brainmicro_cormat))
cormat=matrix(NA,ncol=ncol(gtex_cormat)+ncol(brainmicro_cormat),nrow=length(genes),dimnames=list(genes,c(paste('GTEx',colnames(gtex_cormat),sep='_'),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep='')))))
cormat[rownames(gtex_cormat),paste('GTEx',colnames(gtex_cormat),sep='_')]=gtex_cormat
cormat[rownames(brainmicro_cormat),sapply(strsplit(colnames(brainmicro_cormat),'_'),function(x)paste(x[1],'_Brain-',x[2],sep=''))]=brainmicro_cormat

rmat=cormat

dset=sapply(strsplit(colnames(cormat),'_'),function(x)x[1])
majtis=sapply(strsplit(sapply(strsplit(colnames(cormat),'_'),function(x)x[2]),'-'),function(x)x[1])

xx=sapply(genelist,function(x)cormat[x,])
x=Reduce(rbind,xx)
annotcol=data.frame(TissueType=majtis,Dataset=dset)
rownames(annotcol) <- colnames(x)
TissueType=c(brewer.pal(9,'Set1'),brewer.pal(8,'Dark2'),brewer.pal(8,'Accent'))[1:length(unique(majtis))]
names(TissueType) <-unique(annotcol$TissueType)
Dataset=brewer.pal(8,'Pastel2')
names(Dataset) <-unique(annotcol$Dataset)
rownames(x)=1:nrow(x)
annotrow=data.frame(GeneSet=rep(gsub('\n',' ',names(xx)),sapply(xx,nrow)))
rownames(annotrow)=rownames(x)
GeneSet=muted(brewer.pal(4,'Set1'))
names(GeneSet)=unique(annotrow$GeneSet)
anno_colors <- list(TissueType = TissueType,Dataset=Dataset,GeneSet=GeneSet)
pdf('./results/microarray_vs_gtex_genesets-nobrain-_expression_heatmap.pdf',width=11,height=15)
pheatmap(x,
         annotation_col = annotcol,annotation_row = annotrow,
         annotation_colors = anno_colors,
         color = colorRampPalette(c(muted('lightblue'),'white',muted('pink')))(30),
         annotation_names_col = F,
         annotation_names_row = F,
         show_colnames = T,
         show_rownames = F, 
         breaks=seq(-1,1,length.out = 31),cluster_rows = F,gaps_row = cumsum(unname(sapply(xx,nrow)))[-4])
dev.off()

####################
rm(list=ls())

source('../shared/functions/functions.R')
br_up=readRDS('./data/processed/humanBrainMicroarray/all_up_ensembl.rds')
br_down=readRDS('./data/processed/humanBrainMicroarray/all_down_ensembl.rds')
br_bg=readRDS('./data/processed/humanBrainMicroarray/all_bg_ensembl.rds')
gtex_up=readRDS('./data/processed/GTEx/noBrain_up.rds')
gtex_down=readRDS('./data/processed/GTEx/noBrain_down.rds')
gtex_bg=readRDS('./data/processed/GTEx/noBrain_bg_ensembl.rds')

genelist=setNames(rep(0,length(br_bg)),br_bg)
genelist[br_up]=1
br_up_GO=lapply(c('BP','MF','CC'),function(ontx)go_enrich.test(genelist=genelist,selection = function(x)x==1,ontologyx = ontx,nodesize = 10,corrmethod = 'fdr'))
genelist=setNames(rep(0,length(br_bg)),br_bg)
genelist[br_down]=1
br_down_GO=lapply(c('BP','MF','CC'),function(ontx)go_enrich.test(genelist=genelist,selection = function(x)x==1,ontologyx = ontx,nodesize = 10,corrmethod = 'fdr'))

genelist=setNames(rep(0,length(gtex_bg)),gtex_bg)
genelist[gtex_up]=1
gtex_up_GO=lapply(c('BP','MF','CC'),function(ontx)go_enrich.test(genelist=genelist,selection = function(x)x==1,ontologyx = ontx,nodesize = 10,corrmethod = 'fdr'))
genelist=setNames(rep(0,length(gtex_bg)),gtex_bg)
genelist[gtex_down]=1
gtex_down_GO=lapply(c('BP','MF','CC'),function(ontx)go_enrich.test(genelist=genelist,selection = function(x)x==1,ontologyx = ontx,nodesize = 10,corrmethod = 'fdr'))

xx=list(br_up_GO,br_down_GO,gtex_up_GO,gtex_down_GO)
signifGOs=unique(unlist(lapply(xx,function(x)unlist(lapply(x,function(a)a[a$p.adjusted<0.05,2])))))
genelist=Reduce(union,list(br_up,br_down,gtex_up,gtex_down))
gomat=matrix(NA,ncol=2,nrow=0)

library(reshape2)

for(x in xx){
  for (a in x){
    aa=a$genelist[a$Term%in%signifGOs]
    names(aa)=a$Term[a$Term%in%signifGOs]
    gomat=rbind(gomat,melt(aa))
  }
}

gomat=unique(gomat)
library(igraph)

gomat=graph_from_data_frame(gomat)

library(GGally)

cx=as.factor(1+grepl('ENSG',rownames(gomat[])))

dx=rep(1,nrow(gomat[]))
dx[rownames(gomat[])%in%union(br_up,br_down)]=2
dx[rownames(gomat[])%in%union(gtex_up,gtex_down)]=3
dx[rownames(gomat[])%in%intersect(union(br_up,br_down),union(gtex_up,gtex_down))]=4

go_br_up=unname(unlist(sapply(br_up_GO,function(x)x$Term[x$p.adjusted<0.05])))
go_br_down=unname(unlist(sapply(br_down_GO,function(x)x$Term[x$p.adjusted<0.05])))
go_gtex_up=unname(unlist(sapply(gtex_up_GO,function(x)x$Term[x$p.adjusted<0.05])))
go_gtex_down=unname(unlist(sapply(gtex_down_GO,function(x)x$Term[x$p.adjusted<0.05])))

dx[rownames(gomat[])%in%union(go_br_up,go_br_down)]=2
dx[rownames(gomat[])%in%union(go_gtex_up,go_gtex_down)]=3
dx[rownames(gomat[])%in%intersect(union(go_br_up,go_br_down),union(go_gtex_down,go_gtex_up))]=3


nmx=rownames(gomat[])
nmx[grep('ENSG',nmx)]=''

dx=factor(c('Brain Datasets','GTEx','Both')[dx-1],levels=c('Brain Datasets','GTEx','Both'))
pdf('./results/GO-result.pdf',width = 35,height = 35)
ggnet2(gomat,label = nmx,label.alpha = 1,label.color = 'gray35',
       shape = cx,
       color = dx,alpha = 0.75,edge.color = 'gray90',size=c(50,16)[cx],color.palette = 'Set1')+guides(color=guide_legend(title=''),shape=F,size=F)
dev.off()

