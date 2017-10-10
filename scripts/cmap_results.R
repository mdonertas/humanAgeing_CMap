library(scales)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
rm(list=ls())
cmap=list(microbrain=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv'),
microbrain_noLu=read.csv('./data/processed/humanBrainMicroarray/CMap_results/noLu_consistent_byname.csv'),
gtex=read.csv('./data/processed/GTEx/CMap_results/all_consistent_byname.csv'),
gtex_nobrain=read.csv('./data/processed/GTEx/CMap_results/nobrain_consistent_byname.csv'),
allbrain=read.csv('./data/processed/allBrain/CMap_results/allbrain.csv'),
allbrain_noLu=read.csv('./data/processed/allBrain/CMap_results/allbrain_noLu.csv'))

cmap_padj=lapply(cmap,function(x){
  ps=as.numeric(as.character(x$p))
  padj=p.adjust(ps,method='fdr')
  cbind(x,padj)
})
druglist=unique(unname(unlist(lapply(cmap_padj,function(x)as.character(x$cmap.name[which(x$padj<0.05)])))))

mx=sapply(cmap_padj,function(x){
  setNames(x$mean,x$cmap.name)[druglist]
})
px=sapply(cmap_padj,function(x){
  setNames(x$padj,x$cmap.name)[druglist]
})<0.05


i=hclust(dist((mx)))$order
j=hclust(dist(t(mx)))$order
mx=mx[i,j]
px=px[i,j]
mydat=melt(mx)
colnames(mydat)=c('drug','dataset','mean')
mydat$comb=paste(mydat$drug,mydat$dataset,sep='_')
foo=melt(px)
colnames(foo)=c('drug','dataset','mean')
foo$comb=paste(foo$drug,foo$dataset,sep='_')
rownames(foo)=foo$comb
mydat$signif=foo[mydat$comb,3]

ggplot(mydat,aes(x=dataset,y=drug,fill=mean))+
  geom_tile(lwd=0.2,color='black')+
  geom_text(fontface=c(1,2)[mydat$signif+1],aes(label=round(mean,2)))+
  # scale_fill_gradient2()+
  scale_fill_gradient2(low=muted('blue'),high=muted('red'),midpoint = 0)+
  scale_color_manual(values = c(NA,'gray20',NA))+
  guides(color=F)+
  theme_bw()+
  theme(axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),axis.ticks = element_blank(),axis.text.x = element_text(size=10,face='bold'))+
  scale_x_discrete(position = "top",labels=c('GTEx\nNo Brain','GTEx\nAll','Array\nAll','Array\nNo Lu','GTEx+Array\nAll','GTEx+Array\nNo Lu')) 

pc=prcomp((mx))
biplot(pc)




library(pheatmap)
# rowannot=data.frame(numSig=rowSums(px,na.rm=T))
# rownames(rowannot)=rownames(px)
rowannot=data.frame(AgeingDrug=c('-','+')[1+grepl('sirolimus|wortmannin|geldanamycin|LY-294002|thioridazine|resveratrol',rownames(px),ignore.case = T)])
rownames(rowannot)=rownames(px)
px2=px
px2[is.na(px)]=F
px2[px2==T]=2
px2[px2==F]=1

colannot=data.frame(GTEx=c('No GTEx','GTEx')[grepl('gtex|all',colnames(mx))+1],
           Microarray=c('No Microarray','Microarray')[1+grepl('micro|all',colnames(mx))],
           Brain=c('NoBrain','Brain')[1+!grepl('nobrain',colnames(mx))],
           otherTissues=c('OtherTissues','OtherTissues',rep('Only Brain',4)))
annocols=list(AgeingDrug=setNames(c('white','gray25'),c('-','+')),
              GTEx=setNames(c('white',brewer.pal(3,'Set1')[2]),c('No GTEx','GTEx')),
              Microarray=setNames(c('white',brewer.pal(3,'Set1')[2]),c('No Microarray','Microarray')),
              Brain=setNames(c('white',brewer.pal(3,'Set1')[1]),c('NoBrain','Brain')),
              otherTissues=setNames(c('white',brewer.pal(3,'Set1')[1]),c('Only Brain','OtherTissues')),noLu=setNames(c('white','gray15'),c('Lu','noLu')))

rownames(colannot)=colnames(mx)
pheatmap(mx,
         annotation_row = rowannot,
         annotation_col = colannot,
         breaks=seq(-1,1,length.out = 51),
         color=colorRampPalette(c(muted('blue'),'white',muted('red')))(50),
         display_numbers = T,
         number_format = '%.2f',
         # number_color='gray15',
         fontsize_number =c(6,9)[melt(px2)[,3]],
         number_color =c('gray30','black')[melt(px2)[,3]],
         annotation_colors = annocols,
         show_colnames = F,
         cellwidth = 45,
         cellheight = 11,
         cutree_rows = 7,filename = '~/Desktop/all.pdf')

mx3=mx[names(which(rowSums(px2[,-c(3,5)]==2)>0)),-c(3,5)]
px3=px2[rownames(mx3),-c(3,5)]
px3=px3[hclust(dist(mx3))$order,hclust(dist(t(mx3)))$order]
mx3=mx3[hclust(dist(mx3))$order,hclust(dist(t(mx3)))$order]
pheatmap(mx3,
         annotation_row = rowannot,
         annotation_col = colannot,
         breaks=seq(-1,1,length.out = 51),
         color=colorRampPalette(c(muted('blue'),'white',muted('red')))(50),
         display_numbers = T,
         number_format = '%.2f',
         # number_color='gray15',
         fontsize_number =c(6,9)[melt(px3)[,3]],
         number_color =c('gray30','black')[melt(px3)[,3]],
         annotation_colors = annocols,
         show_colnames = F,
         cellwidth = 45,
         cellheight = 11,
         cutree_rows = 7,filename = '~/Desktop/noLu.pdf')

mx3=mx[names(which(rowSums(px2[,-c(4,6)]==2)>0)),-c(4,6)]
px3=px2[rownames(mx3),-c(4,6)]
trow=hclust(dist(mx3))$order
crow=hclust(dist(t(mx3)))$order
px3=px3[trow,crow]
mx3=mx3[trow,crow]
pheatmap(mx3,
         annotation_row = rowannot,
         annotation_col = colannot,
         breaks=seq(-1,1,length.out = 51),
         color=colorRampPalette(c(muted('blue'),'white',muted('red')))(50),
         display_numbers = T,
         number_format = '%.2f',
         # number_color='gray15',
         fontsize_number =c(6,9)[melt(px3)[,3]],
         number_color =c('gray30','black')[melt(px3)[,3]],
         annotation_colors = annocols,
         show_colnames = F,
         cellwidth = 45,
         cellheight = 11,
         cutree_rows = 7,filename = '~/Desktop/withLu.pdf')
