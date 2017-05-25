## Libraries
source('../shared/functions/functions.R')
library(VennDiagram)
library(scales)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
###############

cormat=readRDS('./data/processed/allBrain/cormat.rds')
nms=colnames(cormat)
bregions=read.table('./docs/bregions.txt',header = T)
subregx=sapply(strsplit(nms,'-'),function(x)x[2])
subregx[subregx%in%bregions$Subregion]=as.character(bregions$Abb[which(bregions$Subregion%in%subregx)])
bregx=setNames(as.character(bregions$BrainRegion),as.character(bregions$Abb))[subregx]
bregx=gsub('Lobe_Cerebrum','Cortex',bregx)
bregx[bregx=='Cortex']='FrontalCortex'
bregx[bregx=='OccipitalCortex']='Other'
dsetname=sapply(strsplit(nms,'_'),function(x)x[1])
colnames(cormat)=paste(dsetname,subregx,sep='_')
DataSourceColors=setNames(brewer.pal(8,'Set2'),unique(dsetname))
BregionColors=setNames(brewer.pal(6,'Set3'),unique(bregx))
array_ages=sapply(grep('Colantuoni',grep('_age.rds',list.files('./data/processed/ages/',full.names = T),v=T),v=T,invert=T),function(x)readRDS(x)/365)
array_ages$Colantuoni_PFC_age=readRDS(grep('Colantuoni',grep('_age.rds',list.files('./data/processed/ages/',full.names = T),v=T),v=T))
summary(unlist(array_ages))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 20.00   30.50   40.61   48.65   64.91  106.00 
length(unique(subregx[dsetname=='GTEx']))
# 13
length(unique(subregx[dsetname!='GTEx']))
# 22


### correlation between datasets
corx_all=cor(cormat,method='s',use='pairwise')
corx_gtex=cor(cormat[,which(dsetname=='GTEx')],method='s',use='pairwise')
corx_array=cor(cormat[,which(dsetname!='GTEx')],method='s',use='pairwise')
annotx=data.frame(DataSource=dsetname,BrainRegion=bregx)
rownames(annotx)=colnames(corx_all)
pheatmap(corx_all,
         annotation_row = annotx,
         annotation_col = annotx,
         annotation_colors = list(DataSource=DataSourceColors,
                                  BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(-1,1,length.out = 31),
         show_rownames = F,
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = T,
         display_numbers = T,
         number_format = '%.1f',
         filename = './results/article/all_correlation.pdf',
         number_color = 'gray15',cellwidth = 15,cellheight = 15)

pheatmap(corx_array,
         annotation_row = annotx,
         annotation_col = annotx,
         annotation_colors = list(DataSource=DataSourceColors,
                                  BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(-1,1,length.out = 31),
         show_rownames = F,
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = T,
         display_numbers = T,
         number_format = '%.1f',
         filename = './results/article/array_correlation.pdf',
         number_color = 'gray15',cellwidth = 15,cellheight = 15)

annotx=data.frame(BrainRegion=bregx)
rownames(annotx)=colnames(corx_all)
pheatmap(corx_gtex,
         annotation_row = annotx,
         annotation_col = annotx,
         annotation_colors = list(DataSource=DataSourceColors,
                                  BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(-1,1,length.out = 31),
         show_rownames = F,
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = T,
         display_numbers = T,
         number_format = '%.1f',
         filename = './results/article/gtex_correlation.pdf',
         number_color = 'gray15',cellwidth = 15,cellheight = 15)

complete_cormat=cormat[complete.cases(cormat),]
dim(complete_cormat)
annotx=data.frame(DataSource=dsetname,BrainRegion=bregx)
rownames(annotx)=colnames(corx_all)
pheatmap(complete_cormat,
         annotation_col = annotx,
         annotation_colors = list(DataSource=DataSourceColors,
                                  BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(-1,1,length.out = 31),
         show_rownames = F,
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = T,
         kmeans_k = 100,
         # border_color = 'gray80',
         cellwidth = 5,cellheight = 5,
         filename='./results/article/all_heatmap.pdf')

complete_cormat_noLu=cormat[complete.cases(cormat[,which(dsetname!='Lu2004')]),which(dsetname!='Lu2004')]
dim(complete_cormat_noLu)
pheatmap(complete_cormat_noLu,
         annotation_col = annotx,
         annotation_colors = list(DataSource=DataSourceColors,
                                  BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(-1,1,length.out = 31),
         show_rownames = F,
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = T,
         kmeans_k = 100,
         # border_color = 'gray80',
         cellwidth = 5,cellheight = 5,
         filename='./results/article/all_heatmap_noLu.pdf')

complete_cormat_bin=cormat[complete.cases(cormat),]>0
annotx=data.frame(DataSource=dsetname,BrainRegion=bregx)
rownames(annotx)=colnames(corx_all)
pheatmap(complete_cormat_bin,
         annotation_col = annotx,
         annotation_colors = list(DataSource=DataSourceColors,
                                  BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(0,1,length.out = 31),
         show_rownames = F,
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = T,
         kmeans_k = 50,
         # border_color = 'gray80',
         cellwidth = 5,cellheight = 5,
         filename='./results/article/all_heatmap_direction.pdf')

complete_cormat_noLu_bin=cormat[complete.cases(cormat[,which(dsetname!='Lu2004')]),which(dsetname!='Lu2004')]>0
pheatmap(complete_cormat_noLu_bin,
         annotation_col = annotx,
         annotation_colors = list(DataSource=DataSourceColors,
                                  BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(0,1,length.out = 31),
         show_rownames = F,
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = T,
         kmeans_k = 50,
         # border_color = 'gray80',
         cellwidth = 5,cellheight = 5,
         filename='./results/article/all_heatmap_noLu_direction.pdf')

array_cormat=cormat[,which(dsetname!='GTEx')]
gtex_cormat=cormat[,which(dsetname=='GTEx')]

array_up=names(which(rowMeans(array_cormat>0)==1))
l_array_up=length(array_up)
#100
array_down=names(which(rowMeans(array_cormat<0)==1))
l_array_down=length(array_down)
#117
array_bg=names(which(rowMeans(is.na(array_cormat))==0))
l_array_bg=length(array_bg)
#5677

gtex_up=names(which(rowMeans(gtex_cormat>0)==1))
l_gtex_up=length(gtex_up)
#1189
gtex_down=names(which(rowMeans(gtex_cormat<0)==1))
l_gtex_down=length(gtex_down)
#1352
gtex_bg=names(which(rowMeans(is.na(gtex_cormat))==0))
l_gtex_bg=length(gtex_bg)
#23718

all_up=names(which(rowMeans(cormat>0)==1))
l_all_up=length(all_up)
#50
all_down=names(which(rowMeans(cormat<0)==1))
l_all_down=length(all_down)
#48
all_bg=names(which(rowMeans(is.na(cormat))==0))
l_all_bg=length(all_bg)
#5421


binom.test(p=l_array_up/l_array_bg * l_gtex_up/l_gtex_bg,x = l_all_up,n = l_all_bg)
# Exact binomial test
# 
# data:  l_all_up and l_all_bg
# number of successes = 50, number of trials = 5421, p-value < 2.2e-16
# alternative hypothesis: true probability of success is not equal to 0.0008830492
# 95 percent confidence interval:
#   0.006853309 0.012142015
# sample estimates:
#   probability of success 
# 0.009223391

binom.test(p=l_array_down/l_array_bg * l_gtex_down/l_gtex_bg,x = l_all_down,n = l_all_bg)
# Exact binomial test
# 
# data:  l_all_down and l_all_bg
# number of successes = 48, number of trials = 5421, p-value < 2.2e-16
# alternative hypothesis: true probability of success is not equal to 0.001174804
# 95 percent confidence interval:
#   0.006535587 0.011722767
# sample estimates:
#   probability of success 
# 0.008854455 

overlapplot=list()
overlapplot[['Microarray\nUp Genes']]=array_up
overlapplot[['GTEx\nUp Genes']]=gtex_up
overlapplot[['Microarray\nDown Genes']]=array_down
overlapplot[['GTEx\nDown Genes']]=gtex_down
venn.diagram(overlapplot,
             filename = './results/article/overlapplot.tiff',
             col=c('firebrick3','forestgreen','forestgreen','firebrick3'))

array_genelist=setNames(rep(0,l_array_bg),array_bg)
array_genelist[array_up]=1
array_genelist[array_down]=-1
array_go_up=go_enrich.test(genelist = array_genelist,selection = function(x)x==1,ontologyx = 'BP',nodesize = 10,corrmethod = 'fdr')
pdf('./results/article/array_go_up.pdf',width = 10,height = 10)
plot_goenrichment(array_go_up)
dev.off()
array_go_down=go_enrich.test(genelist = array_genelist,selection = function(x)x==-1,ontologyx = 'BP',nodesize = 10,corrmethod = 'fdr')
pdf('./results/article/array_go_down.pdf',width = 10,height = 10)
plot_goenrichment(array_go_down)
dev.off()

gtex_genelist=setNames(rep(0,l_gtex_bg),gtex_bg)
gtex_genelist[gtex_up]=1
gtex_genelist[gtex_down]=-1
gtex_go_up=go_enrich.test(genelist = gtex_genelist,selection = function(x)x==1,ontologyx = 'BP',nodesize = 10,corrmethod = 'fdr')
# pdf('./results/article/gtex_go_up.pdf',width = 10,height = 10)
plot_goenrichment(gtex_go_up,label = F)
# dev.off()
gtex_go_down=go_enrich.test(genelist = gtex_genelist,selection = function(x)x==-1,ontologyx = 'BP',nodesize = 10,corrmethod = 'fdr')
# pdf('./results/article/gtex_go_down.pdf',width = 10,height = 10)
plot_goenrichment(gtex_go_down,label = F)
# dev.off()

all_genelist=setNames(rep(0,l_all_bg),all_bg)
all_genelist[all_up]=1
all_genelist[all_down]=-1
all_go_up=go_enrich.test(genelist = all_genelist,selection = function(x)x==1,ontologyx = 'BP',nodesize = 10,corrmethod = 'fdr')
pdf('./results/article/all_go_up.pdf',width = 10,height = 10)
plot_goenrichment(all_go_up,label = F)
dev.off()
all_go_down=go_enrich.test(genelist = all_genelist,selection = function(x)x==-1,ontologyx = 'BP',nodesize = 10,corrmethod = 'fdr')
pdf('./results/article/all_go_down.pdf',width = 10,height = 10)
plot_goenrichment(all_go_down)
dev.off()

cmap_array=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
cmap_gtex=read.csv('./data/processed/brainGTEx/CMap_results/GTEx_brain.csv')
cmap_all=read.csv('./data/processed/allBrain/CMap_results/allbrain.csv')
cmap=list(cmap_array,cmap_gtex,cmap_all)
cmap=lapply(cmap,function(cm){
  ps=as.numeric(as.character(cm$p))
  cm$p.adj=p.adjust(ps,method='fdr')
  cm
})
cmap_array=cmap[[1]]
cmap_gtex=cmap[[2]]
cmap_all=cmap[[3]]

s_cmap_array=cmap_array[which(cmap_array$p.adj<0.05),]
s_cmap_gtex=cmap_gtex[which(cmap_gtex$p.adj<0.05),]
s_cmap_all=cmap_all[which(cmap_all$p.adj<0.05),]

s_cmap_array=s_cmap_array[order(-abs(s_cmap_array$mean)),]
s_cmap_gtex=s_cmap_gtex[order(-abs(s_cmap_gtex$mean)),]
s_cmap_all=s_cmap_all[order(-abs(s_cmap_all$mean)),]

s_array_cids=sapply(as.character(s_cmap_array$cmap.name),get_PubChemCID)
names(s_array_cids)=as.character(s_cmap_array$cmap.name)
s_array_cids=melt(s_array_cids[!sapply(s_array_cids,is.null)])
colnames(s_array_cids)=c('CID','Name')
s_array_chembl=sapply(s_array_cids$CID,getChembl_from_CID)
names(s_array_chembl)=as.character(s_array_cids$CID)
s_array_chembl=melt(s_array_chembl[!sapply(s_array_chembl,is.null)])
colnames(s_array_chembl)=c('CHEMBL','CID')
s_array_id=merge(s_array_chembl,s_array_cids,by='CID',all=T)
s_array_targets=lapply(as.character(unique(s_array_id$CHEMBL[complete.cases(s_array_id$CHEMBL)])),getTarget_from_ChemblID)
head(s_array_targets)

#######################
# save(list=ls(),file='./data/article.RData')
