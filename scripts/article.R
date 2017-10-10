## Libraries
source('../shared/functions/functions.R')
pdf.options(useDingbats=F)
library(VennDiagram)
library(scales)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(reshape2)
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
array_ages$Colantuoni2011_PFC_age=readRDS(grep('Colantuoni',grep('_age.rds',list.files('./data/processed/ages/',full.names = T),v=T),v=T))
summary(unlist(array_ages))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 20.00   30.50   40.61   48.65   64.91  106.00 
length(unique(subregx[dsetname=='GTEx']))
# 13
length(unique(subregx[dsetname!='GTEx']))
# 22

agedist_array=melt(array_ages)
colnames(agedist_array)=c('Age','fname')
agedist_array$dataset=gsub('_age|.rds','',sapply(strsplit(agedist_array$fname,'/'),function(x)x[length(x)]))
agedist_array$dataset=gsub('Somel2010_PFC','Somel2011_PFC',agedist_array$dataset)
head(agedist_array)

pdf('./results/article/array_agedist.pdf',width = 10)
ggplot(agedist_array,aes(x=Age))+
  geom_histogram(bins = 10,color='gray15')+
  facet_wrap(~dataset,scales = 'free_y')+
  theme_bw()+
  ylab('Frequency')
dev.off()

agedist_array$dname=as.factor(sapply(strsplit(agedist_array$dataset,'_'),function(x)x[1]))
pdf('./results/article/array_agedist_summary.pdf',height = 4.5)
ggplot(agedist_array,aes(x=Age,group=dataset,color=dname))+
  stat_density(geom='line',position = 'identity',size=1)+
  # geom_histogram(bins = 10,color='gray15')+
  # facet_wrap(~dataset,scales = 'free_y')+
  theme_bw()+
  ylab('Density')+
  scale_color_manual(values = DataSourceColors[levels(agedist_array$dname)],
                     guide = guide_legend('Dataset',override.aes = list(size=2)))
dev.off()

gtex_sample_att=readRDS('../GTEx/data/processed/all_sample_attributes.rds')
gtex_sample_age=setNames(gtex_sample_att$age,gtex_sample_att$sample_id)
gtex_sample_age=melt(sapply(grep('Brain',list.files('../GTEx/data/processed/exp2/',full.names = T),v=T),function(x)gtex_sample_age[colnames(readRDS(x))]))
gtex_sample_age$dataset=gsub('.rds','',sapply(strsplit(gtex_sample_age$L1,'/'),function(x)x[length(x)]))
gtex_sample_age$dataset=paste('GTEx',unname(setNames(bregions$Abb,bregions$Subregion)[sapply(strsplit(gtex_sample_age$dataset,'-'),function(x)x[2])]),sep='_')
head(gtex_sample_age)
colnames(gtex_sample_age)=c('Age','fname','dataset')
pdf('./results/article/gtex_agedist.pdf',width = 10)
ggplot(gtex_sample_age,aes(x=Age))+
  geom_bar(color='gray15')+
  facet_wrap(~dataset,scales = 'free_y')+
  theme_bw()+
  ylab('Frequency')+
  theme(axis.text.x=element_text(angle=330))
dev.off()

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
         annotation_row = annotx[grep('GTEx',rownames(annotx),invert=T),],
         annotation_col = annotx[grep('GTEx',rownames(annotx),invert=T),],
         annotation_colors = list(DataSource=DataSourceColors[grep('GTEx',names(DataSourceColors),invert=T)],
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
         filename='./results/article/all_heatmap.pdf',
         cellwidth = 5,cellheight = 5)

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

array_permutations=t(readRDS('./data/processed/allBrain/permutations/micro_shared.rds'))
colnames(array_permutations)=c('up','down')

hist(array_permutations[,1])
mean(array_permutations[,1]>=l_array_up)
mean(array_permutations[,2]>=l_array_down)

array_permutations=data.frame(Up=array_permutations[,1],Down=array_permutations[,2])
array_permutations$UpCol=array_permutations$Up>=l_array_up
array_permutations$DownCol=array_permutations$Down>=l_array_down
p1=ggplot(array_permutations,aes(x=Up))+
  geom_histogram(color='gray50',bins=200)+
  geom_segment(arrow = arrow(length = unit(0.1,'inches')),x=l_array_up,xend=l_array_up,y=100,yend=3,color='darkred',lwd=1)+
  geom_text(x=l_array_up,y=120,label='p<0.001')+ xlab('')+
  theme_bw()+
  ggtitle('Up-Regulated Genes')+
  xlim(-1,max(array_permutations$Up,l_array_up))+
  ylim(0,max(table(array_permutations$Up)))+
  ylab('Count')
p2=ggplot(array_permutations,aes(x=Down))+
  geom_histogram(color='gray50',bins=200)+
  geom_segment(arrow = arrow(length = unit(0.1,'inches')),x=l_array_down,xend=l_array_down,y=100,yend=3,color='darkred',lwd=1)+
  geom_text(x=l_array_down,y=120,label='p<0.001')+ xlab('')+
  theme_bw()+
  ggtitle('Down-Regulated Genes')+
  xlim(-1,max(array_permutations$Down,l_array_down))+
  ylim(0,max(table(array_permutations$Down)))+
  ylab('Count')
pdf('./results/article/array_shared_permutation.pdf')
multiplot(p1,p2)
dev.off()


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
pdf('./results/article/array_go_up.pdf',width = 10,height = 7)
plot_goenrichment(array_go_up,shapesx = c(16,18),alphax = c(1,1))
dev.off()
array_go_down=go_enrich.test(genelist = array_genelist,selection = function(x)x==-1,ontologyx = 'BP',nodesize = 10,corrmethod = 'fdr')
pdf('./results/article/array_go_down.pdf',width = 10,height = 10)
plot_goenrichment(array_go_down,shapesx = c(16,18),alphax = c(1,1))
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
save(list=ls(),file='./data/article.RData')
#######################
melike=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
ps=as.numeric(as.character(melike$p))
melike$padj=p.adjust(ps,method='fdr')
gtex=read.csv('./data/processed/brainGTEx/CMap_results/GTEx_brain.csv')
ps=as.numeric(as.character(gtex$p))
gtex$padj=p.adjust(ps,method='fdr')
matias=read.csv('./docs/matias_list/all.csv')
drosophila=read.csv('./docs/matthias_list/drosophila.csv')
celegans=read.csv('./docs/matthias_list/celegans.csv')
drugage=read.csv('./data/raw/DrugAge/DrugAge.csv')
idconvert2=readRDS('./data/processed/idconvert.rds')

resx=apply(idconvert2,1,function(x){
  list(microarray_brain=melike[melike$cmap.name%in%x[1],],
       gtex=gtex[gtex$cmap.name%in%x[1],],
       matias=matias[matias$drug%in%x[2],],
       matthias_drosophila=drosophila[drosophila$Drug..HetCode.%in%x[3],],
       matthias_celegans=celegans[celegans$Drug..HetCode.%in%x[4],],
       drugage=drugage[drugage$compound_name%in%x[5],])
})
resx=lapply(resx,function(res){
  sapply(res,function(x){
    if(nrow(x)==0){data.frame(t(setNames(rep(NA,ncol(x)),colnames(x))))}
    else(x)
  })
})
resx2=lapply(resx,function(res){
  data.frame(Array_Rank=as.numeric(res$microarray_brain$rank),
             CMap_name=as.character(res$microarray_brain$cmap.name),
             Array_mean=res$microarray_brain$mean,
             N=res$microarray_brain$n,
             Array_Enrichment=res$microarray_brain$enrichment,
             Array_Spesificity=as.numeric(res$microarray_brain$specificity),
             Array_NonNULL=res$microarray_brain$percent.non.null,
             Array_p=res$microarray_brain$padj,
             
             GTEx_Rank=as.numeric(res$gtex$rank),
             GTEx_mean=res$gtex$mean,
             GTEx_Enrichment=res$gtex$enrichment,
             GTEx_Spesificity=as.numeric(res$gtex$specificity),
             GTEx_NonNULL=res$gtex$percent.non.null,
             GTEx_p=res$gtex$padj,
             
             Target=as.character(res$matias$target),
             CHEMBL=as.character(res$matias$drug),
             Phase=res$matias$Max.Phase,
             Gene.Name=as.character(res$matias$Gene.Name),
             orthology=res$matias$WORMHOLE.Score,
             Lifespan=res$matias$Lifespan,
             Lifespan_score=res$matias$Score_target,
             Matias_score=res$matias$Score,
             
             HET=as.character(res$matthias_drosophila$Drug..HetCode.),
             Drosophila_Score=res$matthias_drosophila$Drosophila.Score,
             Drosophila_Domain=res$matthias_drosophila$Domain.conservation,
             Drosophila_BindingSite=res$matthias_drosophila$Binding.site.conservation,
             Drosophila_BindingAffinity=res$matthias_drosophila$Binding.affinity,
             Drosophila_Bioavailability=res$matthias_drosophila$Bio.availability,
             Drosophila_Lipinski=res$matthias_drosophila$Lipinski.Loss,
             Drosophila_Purchasability=res$matthias_drosophila$Purchasability.Bonus,
             Celegans_Score=res$matthias_celegans$C..elegans.Score,
             Celegans_Domain=res$matthias_celegans$Domain.conservation,
             Celegans_BindingSite=res$matthias_celegans$Binding.site.conservation,
             Celegans_BindingAffinity=res$matthias_celegans$Binding.affinity,
             Celegans_Bioavailability=res$matthias_celegans$Bio.availability,
             
             DrugAgeName=res$drugage$compound_name,
             TestedOrganism=res$drugage$species,
             Average_lifespanChange=res$drugage$avg_lifespan_change,
             Max_lifespanChange=res$drugage$max_lifespan_change,
             DrugAge_signif=res$drugage$significance)
})
library(reshape2)
res2=melt(resx2,id.vars=colnames(resx2$`4091`))

res2$Matias_rank=rank(1-res2$Matias_score)
res2$Matias_rank[is.na(res2$Matias_score)]=NA
res2$Drosophila_rank=rank(1-res2$Drosophila_Score)
res2$Drosophila_rank[is.na(res2$Drosophila_Score)]=NA
res2$Celegans_rank=rank(1-res2$Celegans_Score)
res2$Celegans_rank[is.na(res2$Celegans_Score)]=NA

res2$Array_names=as.character(res2$CMap_name)
res2$Array_names[-which(res2$Array_p<=0.05)]=NA
res2$Array_names[-which(res2$Array_p<=0.05)]=res2$DrugAgeName[-which(res2$Array_p<=0.05)]

res2$Array_signif='-'
res2$Array_signif[which(res2$Array_p<=0.05)]='FDR<=0.05'
res2$Array_signif[which(res2$Array_p>0.05)]='NS'
res2$Array_signif=factor(res2$Array_signif,levels=c('-','NS','FDR<=0.05'))

res2$GTEx_names=as.character(res2$CMap_name)
res2$GTEx_names[-which(res2$GTEx_p<=0.05)]=NA
res2$GTEx_names[-which(res2$GTEx_p<=0.05)]=res2$DrugAgeName[-which(res2$GTEx_p<=0.05)]

res2$GTEx_signif='-'
res2$GTEx_signif[which(res2$GTEx_p<=0.05)]='FDR<=0.05'
res2$GTEx_signif[which(res2$GTEx_p>0.05)]='NS'
res2$GTEx_signif=factor(res2$GTEx_signif,levels=c('-','NS','FDR<=0.05'))

res2$CMap_signif='-'
res2$CMap_signif[which(res2$GTEx_p<=0.05)]='GTEx FDR<=0.05'
res2$CMap_signif[which(res2$Array_p<=0.05)]='Array FDR<=0.05'
res2$CMap_signif[which(res2$Array_p<=0.05 & res2$GTEx_p<=0.05)]='Both FDR<=0.05'

res2$CMap_names=as.character(res2$CMap_name)
res2$CMap_names[-which(res2$Array_p<=0.05)]=NA
res2$CMap_names[-which(res2$Array_p<=0.05)]=res2$DrugAgeName[-which(res2$Array_p<=0.05)]
res2$CMap_names[-which(res2$Array_p<=0.05)]=res2$GTEx_names[-which(res2$Array_p<=0.05)]

compare_methods=res2
rm(list=setdiff(ls(),'compare_methods'))
load('./data/article.RData')
save(list=ls(),file='./data/article.RData')
############################
drugsindrugage=data.frame(compare_methods$CMap_name,compare_methods$DrugAgeName)
drugsindrugage=unique(drugsindrugage[complete.cases(drugsindrugage),])
drugsindrugage=as.character(drugsindrugage$compare_methods.CMap_name)
indrugageperm=sapply(1:10000,function(i)sum(as.character(sample(cmap_array$cmap.name,nrow(s_cmap_array)))%in%drugsindrugage))
mean(indrugageperm>=4)
hist(indrugageperm,br=5)
summary(indrugageperm)
save(list=ls(),file='./data/article.RData')
#########################
res3=unique(compare_methods[,c('Array_mean','GTEx_mean','CMap_names','Array_p','DrugAgeName','CMap_signif','Array_signif','GTEx_signif')])
res3$DrugNames=res3$CMap_names
res3$DrugNames[which(!(res3$Array_signif=='FDR<=0.05' | res3$GTEx_signif=='FDR<=0.05'))]=NA
summary(lm(GTEx_mean~Array_mean,data=res3))
pdf('./results/article/array_vs_gtex_drugres.pdf',width = 12)
ggplot(res3,aes(x=Array_mean,y=GTEx_mean))+
  geom_line(stat='smooth',method='lm',color='gray')+
  geom_point(alpha=0.5,aes(color=is.na(DrugAgeName),shape=CMap_signif,size=CMap_signif))+
  geom_text(nudge_y = 0.025,aes(color=is.na(DrugAgeName),label=DrugNames,size=CMap_signif))+
  theme_bw()+
  scale_size_manual(values=c(1,3,3,3))+
  scale_shape_manual(values=c(16,17,18,15))+
  # ggtitle('Microarray vs. GTEx')+
  guides(shape=guide_legend('Significant in Array'),
         size=guide_legend('Significant in Array'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))+
  geom_text(x=0.75,y=0.35,label='y=0.5124x-0.0067',fontface='italic')+
  xlab('Similarity Score - Based on Microarray Datasets')+
  ylab('Similarity Score - Based on GTEx')
dev.off()
################################

br_CMapRes_m=setNames(br_CMapRes$mean,br_CMapRes$cmap.name)
gtex_brain_CMapRes_m=setNames(gtex_brain_CMapRes$mean,gtex_brain_CMapRes$cmap.name)

##############

psx=as.numeric(as.character(br_CMapRes$p))
br_CMapRes$padj=p.adjust(psx,method='fdr')
br_cmap=as.character(br_CMapRes[which(br_CMapRes$padj<0.05),2])
psx=as.numeric(as.character(gtex_brain_CMapRes$p))
gtex_brain_CMapRes$padj=p.adjust(psx,method='fdr')
gtex_brain_cmap=as.character(gtex_brain_CMapRes[which(gtex_brain_CMapRes$padj<0.05),2])
drlist=union(br_cmap,gtex_brain_cmap)
array_brain=data.frame(drname=drlist,mean=br_CMapRes_m[drlist],p=drlist%in%br_cmap,type='Array\nBrain')
gtex_brain=data.frame(drname=drlist,mean=gtex_brain_CMapRes_m[drlist],p=drlist%in%gtex_brain_cmap,type='GTEx\nBrain')
mymat=rbind(array_brain,gtex_brain)
mymat$drname=factor(mymat$drname,levels=drlist[hclust(dist(cbind(br_CMapRes_m[drlist],gtex_brain_CMapRes_m[drlist])))$order])
idconvert=as.data.frame(readRDS('./data/processed/idconvert.rds'))
interlist=setdiff(unique(idconvert$CMap[!is.na(idconvert$DrugAge)]),NA)
# pdf('./results/article/array_vs_gtex_significant.pdf',width = 3,useDingbats = F,height = 3)
# ggplot(mymat,aes(x=type,y=drname,label=round(mean,2),fill=mean,color=p))+
#   geom_tile()+
#   geom_text(aes(size=p),fontface='bold')+
#   scale_color_manual(values = c('gray40','black'),guide=F)+
#   scale_fill_gradient2(midpoint = 0,guide = guide_colorbar('Mean\nSimilarity'))+
#   scale_size_manual(values=c(1,1.5),guide=F)+
#   theme_bw()+
#   theme(axis.title = element_blank(),
#         axis.text = element_text(size=6,face='bold'),
#         panel.border = element_blank(),
#         panel.grid = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_text(colour = c('gray25','darkred')[1+levels(mymat$drname)%in%interlist]),
#         legend.title = element_text(size=5,face='bold'),
#         legend.text = element_text(size=4),
#         legend.key.size = unit(5,units = 'pt'))
# dev.off()

save(list=ls(),file='./data/article.RData')

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
pdf('./results/article/gtex_shared_permutation.pdf')
multiplot(p1,p2)
dev.off()
gtex_permutations=null_sharedness
rm(list=setdiff(ls(),'gtex_permutations'))
load('./data/article.RData')
save(list=ls(),file='./data/article.RData')

indrugageperm_gtex=sapply(1:10000,function(i)sum(as.character(sample(cmap_gtex$cmap.name,nrow(s_cmap_gtex)))%in%drugsindrugage))
mean(indrugageperm_gtex>=sum(s_cmap_gtex$cmap.name%in%drugsindrugage))

library(ggplot2)
pdf('./results/article/drug_intersect.pdf',useDingbats = F,width = 10,height = 10)
plot(x=0:10,y=0:10,type='n',axes = 'n')
draw.circle(x=3,y=5,radius = 3,col = rgb(1,0,0,0.2))
draw.circle(x=6,y=5,radius = 3,col = rgb(0,0,1,0.2))
text(x=4.5,y=seq(3,7,length.out = length(intersect(s_cmap_array$cmap.name,s_cmap_gtex$cmap.name))),
     labels = intersect(s_cmap_array$cmap.name,s_cmap_gtex$cmap.name),
     col=c('gray35','darkred')[1+intersect(s_cmap_array$cmap.name,s_cmap_gtex$cmap.name)%in%drugsindrugage],font=2)
text(x=2,y=seq(3,7,length.out = length(setdiff(s_cmap_array$cmap.name,s_cmap_gtex$cmap.name))),
     labels = setdiff(s_cmap_array$cmap.name,s_cmap_gtex$cmap.name),
     col=c('gray35','darkred')[1+setdiff(s_cmap_array$cmap.name,s_cmap_gtex$cmap.name)%in%drugsindrugage],font=2)
text(x=7.3,y=seq(3,7,length.out = length(setdiff(s_cmap_gtex$cmap.name,s_cmap_array$cmap.name))),
     labels = setdiff(s_cmap_gtex$cmap.name,s_cmap_array$cmap.name),
     col=c('gray35','darkred')[1+setdiff(s_cmap_gtex$cmap.name,s_cmap_array$cmap.name)%in%drugsindrugage],font=2)
text(x=1,y=8.5,labels='Microarray',font=2,cex=2)
text(x=8,y=8.5,labels='GTEx',font=2,cex=2)
dev.off()
# ######################################
# # 06.06.2017
# load('./data/article.RData')
# library(pheatmap)
# library(RColorBrewer)
# library(data.table)
# 
# cmap_detailed_array=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_detailed.csv')
# cmap_detailed_array=setNames(cmap_detailed_array$score,cmap_detailed_array$instance_id)
# drugmat=fread('./data/raw/CMap/amplitudeMatrix.txt')
# inst=read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
# plist=drugmat$V1
# drugmat$V1=NULL
# colnames(drugmat)=strsplit(readLines('./data/raw/CMap/amplitudeMatrix.txt',n=1),'\t')[[1]][-1]
# rownames(drugmat)=plist
# instnames=setNames(inst$cmap_name,inst$instance_id)
# mymat=as.matrix(drugmat[,which(colnames(drugmat)%in%names(instnames)[which(instnames%in%s_cmap_array$cmap.name)]),with=F])
# colannot=data.frame(drugname=inst$cmap_name[which(inst$cmap_name%in%s_cmap_array$cmap.name)],cline=inst$cell2[which(inst$cmap_name%in%s_cmap_array$cmap.name)])
# rownames(colannot)=inst$instance_id[which(inst$cmap_name%in%s_cmap_array$cmap.name)]
# colannot$meanscore=cmap_detailed_array[rownames(colannot)]
# annot_colors=list(drugname=setNames(c(brewer.pal(9,'Set1'),brewer.pal(4,'Set2')),unique(colannot$drugname)),
#                   cline=setNames(brewer.pal(5,'Set3'),unique(colannot$cline))
#                   # meanscore=setNames(colorRampPalette(c('midnightblue','white','darkred'))(2001),as.character(round(seq(-1,1,by = 0.001),3)))[as.character(colannot$meanscore)]
#                   )
# pheatmap(mymat,
#          kmeans_k = 100,
#          annotation_col = colannot,
#          annotation_colors = annot_colors,
#          show_colnames = F,
#          color = colorRampPalette(c('blue','white','red'))(100),
#          breaks = seq(-2,2,length.out = 101),filename = '~/Desktop/deneme.pdf')
# mymat=mymat[,names(sort(setNames(colannot$meanscore,rownames(colannot)),dec=T))]
# ph=pheatmap(mymat,
#          kmeans_k = 10,
#          annotation_col = colannot,
#          annotation_colors = annot_colors,
#          show_colnames = F,
#          cluster_cols = F,
#          color = colorRampPalette(c('blue','white','red'))(100),
#          breaks = seq(-1,1,length.out = 101))
# 
# mygenelist=setNames(rep(0,nrow(drugmat)),rownames(drugmat))
# mygenelist[rownames(drugmat)[which(ph$kmeans$cluster==5)]]=1
# library(biomaRt)
# ids=getBM(attributes = c('ensembl_gene_id','affy_hg_u133a'),filters = 'affy_hg_u133a',values = names(mygenelist),mart = useMart('ensembl','hsapiens_gene_ensembl'))
# ids=setNames(ids$affy_hg_u133a,ids$ensembl_gene_id)
# head(ids)
# mygenelist2=as.numeric(sapply(names(ids),function(ens)all(mygenelist[unname(ids[names(ids)%in%ens])]==1)))
# names(mygenelist2)=names(ids)
# head(mygenelist2)
# gores_negativedrugs=go_enrich.test(genelist = mygenelist2,selection = function(x)x==1,nodesize = 5)
# xx=data.frame(term=gores_negativedrugs$Term,pval=gores_negativedrugs$p.adjusted)
# xx=xx[xx$pval<0.05,]
# br_gene_ids_array_down=readRDS('./data/processed/humanBrainMicroarray/all_down_affy.rds')
# 
# ph=pheatmap(mymat[rownames(drugmat)%in%br_gene_ids_array_down,],
#             kmeans_k = 4,
#             annotation_col = colannot,
#             annotation_colors = annot_colors,
#             show_colnames = F,
#             cluster_cols = F,
#             color = colorRampPalette(c('blue','white','red'))(100),
#             breaks = seq(-1,1,length.out = 101))
# xx=mymat[rownames(drugmat)%in%br_gene_ids_array_down,]
# apply(xx,1,function(x){
#   wilcox.test(x[ptx],x[ngx])$p.val
# })
##########
signifdruglist=unique(c(as.character(s_cmap_array$cmap.name),as.character(s_cmap_gtex$cmap.name)))
matias=sample(signifdruglist,length(signifdruglist)/2)
melike=setdiff(signifdruglist,matias)
matias
write.table(melike,file='./data/processed/drugstosearch.csv',quote = F,row.names = F,col.names = F)
save(list=ls(),file='./data/article.RData')
##############################
x1=setNames(cmap_array$mean,cmap_array$cmap.name)
x2=setNames(cmap_gtex$mean,cmap_gtex$cmap.name)
cmap_other=read.csv('./data/processed/GTEx/CMap_results/nobrain_consistent_byname.csv')
s_cmap_other=cmap_other[which(p.adjust(as.numeric(as.character(cmap_other$p)),method='fdr')<0.05),]
x3=setNames(cmap_other$mean,cmap_other$cmap.name)
# mydruglist=unique(c(as.character(s_cmap_array$cmap.name),as.character(s_cmap_other$cmap.name),as.character(s_cmap_gtex$cmap.name)))
mydruglist=unique(c(as.character(s_cmap_other$cmap.name)))
# mydruglist=unique(c(as.character(s_cmap_array$cmap.name),as.character(s_cmap_other$cmap.name)))

par(mfrow=c(1,2))
cor.test(x3[mydruglist],x1[mydruglist],method='s')
plot(x3[mydruglist],x1[mydruglist],type='n',xlim=c(-1,1),ylim=c(-1,1),main='Other Tissues')
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')
abline(0,1,lty=2,col='gray')
abline(lm(x1[mydruglist]~x3[mydruglist]),col='pink')
text(x3[mydruglist],x1[mydruglist],labels=mydruglist,cex=0.5)

mydruglist=unique(c(as.character(s_cmap_array$cmap.name)))

cor.test(x1[mydruglist],x3[mydruglist],method='s')
plot(x1[mydruglist],x3[mydruglist],type='n',xlim=c(-1,1),ylim=c(-1,1),main='Array')
abline(h=0,lty=2,col='gray')
abline(v=0,lty=2,col='gray')
abline(0,1,lty=2,col='gray')
abline(lm(x3[mydruglist]~x1[mydruglist]),col='pink')
text(x1[mydruglist],x3[mydruglist],labels=mydruglist,cex=0.5)

###################

library(ggplot2)
psx=as.numeric(as.character(br_CMapRes$p))
br_CMapRes$padj=p.adjust(psx,method='fdr')
br_cmap=as.character(br_CMapRes[which(br_CMapRes$padj<0.05),2])
psx=as.numeric(as.character(gtex_brain_CMapRes$p))
gtex_brain_CMapRes$padj=p.adjust(psx,method='fdr')
gtex_brain_cmap=as.character(gtex_brain_CMapRes[which(gtex_brain_CMapRes$padj<0.05),2])
drlist=union(br_cmap,gtex_brain_cmap)
array_brain=data.frame(drname=drlist,mean=br_CMapRes_m[drlist],p=drlist%in%br_cmap,type='Array\nBrain')
gtex_brain=data.frame(drname=drlist,mean=gtex_brain_CMapRes_m[drlist],p=drlist%in%gtex_brain_cmap,type='GTEx\nBrain')
mymat=cbind(array_brain$mean,gtex_brain$mean)
colnames(mymat)=c('Array','GTEx')
rownames(mymat)=array_brain$drname
annotx=data.frame(DrugAge=c('-','in DrugAge')[1+(rownames(mymat)%in%interlist)])
rownames(annotx)=rownames(mymat)
annotxcol=list(DrugAge=setNames(brewer.pal(3,'Set1')[1:2],c('in DrugAge','-')))

pheatmap(mymat,cutree_rows = 2,display_numbers = T,annotation_row = annotx,annotation_names_row = F,annotation_colors = annotxcol,
         fontsize_number = c(c(5,8)[1+array_brain$p[hclust(dist(mymat))$order]],c(5,8)[1+gtex_brain$p[hclust(dist(mymat))$order]]),
         number_color = c(c('gray50','black')[1+array_brain$p[hclust(dist(mymat))$order]],c('gray50','black')[1+gtex_brain$p[hclust(dist(mymat))$order]]),
         filename = './results/article/array_vs_gtex_significant2.pdf',
         cellwidth = 35,cellheight = 15,
         breaks = seq(-0.87,0.87,length.out = 101))
