library(clustermq)
fx = function(i){
  xx=readRDS(paste('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/cormat/perm',i,'.rds',sep=''))
  co=cor(xx,method='s',use='pairwise')
  co
  }

fx2=function(i){readRDS(paste('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/cormat/perm1.rds',sep=''))}
exprs_perm1=Q(fx2,i=1,n_jobs = 1)[[1]]
dim(exprs_perm1)

corsx=Q(fx, i=1:1000, n_jobs=50)

load('./data/article.RData')
rm(list=setdiff(ls(),c('corsx','corx_array','exprs_perm1')))

xx=t(sapply(corsx,function(x){
  co=x[upper.tri(x)]
  c(median(co),mean(co<0))
}))

head(xx)

par(mfrow=c(1,2))
hist(xx[,2],xlim=c(-0.1,0.55),br=100,main='Percent negative')
abline(v=mean(corx_array<0),col='red')

hist(xx[,1],xlim=c(-0.1,0.35),br=100,main='Median Correlation Coefficient')
abline(v=median(corx_array),col='red')


library(pheatmap)
library(scales)

pheatmap(corx_array,
         # annotation_row = annotx[grep('GTEx',rownames(annotx),invert=T),],
         # annotation_col = annotx[grep('GTEx',rownames(annotx),invert=T),],
         # annotation_colors = list(DataSource=DataSourceColors[grep('GTEx',names(DataSourceColors),invert=T)],
                                  # BrainRegion=BregionColors),
         color = colorRampPalette(c(muted('lightblue'),'white','darkred'))(30),
         breaks=seq(-1,1,length.out = 31),
         # show_rownames = F,
         # show_colnames = F,
         # annotation_names_row = F,
         # annotation_names_col = T,
         display_numbers = T,
         number_format = '%.1f',
         # filename = './results/article/array_correlation.pdf',
         number_color = 'gray15',cellwidth = 15,cellheight = 15)


plot(cmdscale(1-corx_array),col=as.factor(dsource),pch=19,cex=2)


sx=rnorm(2000,0,0.3)
x1=c(rnorm(n = 16000,mean = 0,sd = 0.2),sx)
x2=c(rnorm(n = 16000,mean = 0,sd = 0.2),sx)





cor.test(x1,x2,method='s')
plot(rank(x1),rank(x2))


dataset_genelists=apply(!is.na(exprs_perm1),2,function(aa)rownames(exprs_perm1)[aa])

totlen=15000
ageinggene=totlen*15/100

randomvector=sapply(1:1000,function(i){
  sx=rnorm(ageinggene,0,0.3)
  x1=c(rnorm(n = totlen-ageinggene,mean = 0,sd = 0.2),sx)
  x2=c(rnorm(n = totlen-ageinggene,mean = 0,sd = 0.2),sx)
  cor(x1,x2,method='s')
})
hist(randomvector,br=50)

dsource=sapply(strsplit(colnames(corx_array),'_'),function(x)x[[1]])

inds=sapply(dsource,function(x)which(x!=dsource))
median(unlist(sapply(1:nrow(corx_array),function(i){
  corx_array[inds[[i]],i]
})))


500/15000

exprs_cors=readRDS('./data/processed/humanBrainMicroarray/cormat.rds')
head(exprs_cors)

cor(exprs_cors[,1],exprs_cors[,2],method='s',use = 'pairwise')

mycor=exprs_cors[complete.cases(exprs_cors[,1:2]),1:2]

plot((mycor[,1]),(mycor[,2]))
plot(rank(mycor[,1]),rank(mycor[,2]))

mycor=data.frame(mycor)
library(ggplot2)
ggplot(mycor,aes(x=Barnes2011_STC,y=Berchtold2008_EC))+geom_hex()+geom_smooth(method='lm',color='darkred')+theme_bw()


genage_human=readRDS('../shared/ageingGenes/genAge/human/genage_human_ensemblID_20170313.rds')
load('./data/article.RData')
sum(array_up%in%genage)
sum(array_down%in%genage)

down_magal=read.csv('~/Desktop/signatures_supplement/underexp_genage2.csv',skip = 9,head=T)[,1]
up_magal=read.csv('~/Desktop/signatures_supplement/overexp_genage2.csv',skip = 14,head=T)[,1]


library(biomaRt)
down_magal_ens=getBM(attributes = c('ensembl_gene_id'),filters = 'entrezgene',values = down_magal,mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1]
up_magal_ens=getBM(attributes = c('ensembl_gene_id'),filters = 'entrezgene',values = up_magal,mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1]
genage=c(down_magal_ens,up_magal_ens)


down_magal2=read.csv('~/Desktop/signatures_supplement/underexp_genage.csv',skip = 7,head=T)[,1]
up_magal2=read.csv('~/Desktop/signatures_supplement/overexp_genage.csv',skip = 7,head=T)[,1]


library(biomaRt)
down_magal_ens2=getBM(attributes = c('ensembl_gene_id'),filters = 'entrezgene',values = down_magal2,mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1]
up_magal_ens2=getBM(attributes = c('ensembl_gene_id'),filters = 'entrezgene',values = up_magal2,mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1]
genage2=c(down_magal_ens,up_magal_ens2)

library(VennDiagram)
library(tidyverse)

venn.diagram(list(Up=array_up,Down=array_down,GenAge_Down=down_magal_ens,GenAge_Up=up_magal_ens),filename = '~/Desktop/venn2.png',imagetype = 'png')

modelgenage=read.csv('~/Desktop/models_genes/genage_models.csv')
unique(modelgenage$organism)

celegans=(filter(modelgenage,organism==unique(modelgenage$organism)[1]) %>% select(entrez.gene.id))[,1]
mouse=(filter(modelgenage,organism==unique(modelgenage$organism)[2]) %>% select(entrez.gene.id))[,1]
yeast=(filter(modelgenage,organism==unique(modelgenage$organism)[3]) %>% select(entrez.gene.id))[,1]
fly=(filter(modelgenage,organism==unique(modelgenage$organism)[4]) %>% select(entrez.gene.id))[,1]

worm=unique(getBM(attributes = c('hsapiens_homolog_ensembl_gene'),filters = 'entrezgene',values = celegans,mart = useMart('ensembl','celegans_gene_ensembl'))[,1])
mouse=unique(getBM(attributes = c('hsapiens_homolog_ensembl_gene'),filters = 'entrezgene',values = mouse,mart = useMart('ensembl','mmusculus_gene_ensembl'))[,1])
yeast=unique(getBM(attributes = c('hsapiens_homolog_ensembl_gene'),filters = 'entrezgene',values = yeast,mart = useMart('ensembl','scerevisiae_gene_ensembl'))[,1])
fly=unique(getBM(attributes = c('hsapiens_homolog_ensembl_gene'),filters = 'entrezgene',values = fly,mart = useMart('ensembl','dmelanogaster_gene_ensembl'))[,1])

AgeingLiterature=unique(getBM(attributes='ensembl_gene_id',filters = 'entrezgene',values = read.csv('~/Desktop/plosone.csv')[,2],mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1])

# AgeingLiterature2=unique(getBM(attributes='ensembl_gene_id',filters = 'entrezgene',values = read.csv('~/Desktop/plosone2.csv')[,2],mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1])

AgeingLiterature3=unique(getBM(attributes='ensembl_gene_id',filters = 'entrezgene',values = read.csv('~/Desktop/literature3.csv')[,2],mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1])

venn.diagram(list(Array=c(array_up,array_down),GTEx=c(gtex_up,gtex_down),worm=worm,genage=genage),filename = '~/Desktop/venn3.png',imagetype = 'png')

AgeFactDB=read.table('~/Desktop/AgeFactDB/gene-expmAnalysis.csv',head=T,sep=',')
AgeFactDB=as.character(AgeFactDB[AgeFactDB$Name..Species...Homology.Group.Gene.=='Homo sapiens',1])
AgeFactDB=unique(getBM(attributes='ensembl_gene_id',filters = 'hgnc_symbol',values = AgeFactDB,mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1])


atlas=unique(getBM(attributes='ensembl_gene_id',filters = 'hgnc_symbol',values = unique(sapply(strsplit(as.character(read.table('~/Desktop/digital_ageing_atlas_data/atlas_brain.txt',head=T,sep='\t')[,1]),'[ (]'),function(x)x[1])),mart = useMart('ensembl','hsapiens_gene_ensembl'))[,1])


mylist=list(
  # AgeingLiterature=AgeingLiterature,
            # AgeingLiterature2=AgeingLiterature2,
            # AgeingLiterature3=AgeingLiterature3,
  AgeFactDB=AgeFactDB,
            Array=c(array_up,array_down),
            # GTEx=c(gtex_up,gtex_down),
  atlas=atlas,
            GenAge_Exp1=genage,
            GenAge_Exp2=genage2,
            GenAge_Human=genage_human,
            Mouse=mouse,
            Fly=fly,
            Worm=worm,
            Yeast=yeast
)
names(mylist)=paste(names(mylist),unname(sapply(mylist,length)),sep='_')  
inters=sapply(mylist,function(x)sapply(mylist,function(y)mean(x%in%y)))
pheatmap(inters,
         # annotation_row = annotx[grep('GTEx',rownames(annotx),invert=T),],
         # annotation_col = annotx[grep('GTEx',rownames(annotx),invert=T),],
         # annotation_colors = list(DataSource=DataSourceColors[grep('GTEx',names(DataSourceColors),invert=T)],
         # BrainRegion=BregionColors),
         color = colorRampPalette(c('white','darkred'))(30),
         breaks=seq(0,1,length.out = 31),
         # show_rownames = F,
         # show_colnames = F,
         # annotation_names_row = F,
         # annotation_names_col = T,
         display_numbers = T,
         number_format = '%.2f',
         # filename = './results/article/array_correlation.pdf',
         number_color = 'gray15',cellwidth = 45,cellheight = 45,cluster_cols = F,cluster_rows = F)

venn.diagram(mylist,file='~/Desktop/human2.png',imagetype = 'png',col = "transparent",
             fill = c("cornflowerblue", "green", "yellow", "darkorchid1"))
genlist=Reduce(f = 'union',mylist)
mymat=matrix(0,nrow=length(genlist),ncol=length(mylist),dimnames = list(genlist,names(mylist)))
for(nm in names(mylist)){
  mymat[mylist[[nm]],nm]=1
}

# pheatmap(mymat,kmeans_k = 25)