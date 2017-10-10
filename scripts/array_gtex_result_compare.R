rm(list=ls())
#################
library(reshape2)
library(ggplot2)
#################
br_CMapRes=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
gtex_brain_CMapRes=read.csv('./data/processed/brainGTEx/CMap_results/GTEx_brain.csv')
gtex_nobrain_CMapRes=read.csv('./data/processed/GTEx/CMap_results/nobrain_consistent_byname.csv')
#################

br_CMapRes_m=setNames(br_CMapRes$mean,br_CMapRes$cmap.name)
gtex_brain_CMapRes_m=setNames(gtex_brain_CMapRes$mean,gtex_brain_CMapRes$cmap.name)
gtex_nobrain_CMapRes_m=setNames(gtex_nobrain_CMapRes$mean,gtex_nobrain_CMapRes$cmap.name)

##############

psx=as.numeric(as.character(br_CMapRes$p))
br_CMapRes$padj=p.adjust(psx,method='fdr')
br_cmap=as.character(br_CMapRes[which(br_CMapRes$padj<0.05),2])
psx=as.numeric(as.character(gtex_brain_CMapRes$p))
gtex_brain_CMapRes$padj=p.adjust(psx,method='fdr')
gtex_brain_cmap=as.character(gtex_brain_CMapRes[which(gtex_brain_CMapRes$padj<0.05),2])
psx=as.numeric(as.character(gtex_nobrain_CMapRes$p))
gtex_nobrain_CMapRes$padj=p.adjust(psx,method='fdr')
gtex_nobrain_cmap=as.character(gtex_nobrain_CMapRes[which(gtex_nobrain_CMapRes$padj<0.05),2])

drlist=union(union(br_cmap,gtex_brain_cmap),gtex_nobrain_cmap)

array_brain=data.frame(drname=drlist,mean=br_CMapRes_m[drlist],p=drlist%in%br_cmap,type='Array\nBrain')
gtex_brain=data.frame(drname=drlist,mean=gtex_brain_CMapRes_m[drlist],p=drlist%in%gtex_brain_cmap,type='GTEx\nBrain')
gtex_other=data.frame(drname=drlist,mean=gtex_nobrain_CMapRes_m[drlist],p=drlist%in%gtex_nobrain_cmap,type='GTEx\nOther')

mymat=rbind(array_brain,gtex_brain,gtex_other)
mymat$drname=factor(mymat$drname,levels=drlist[hclust(dist(cbind(br_CMapRes_m[drlist],gtex_brain_CMapRes_m[drlist],gtex_nobrain_CMapRes_m[drlist])))$order])

idconvert=as.data.frame(readRDS('./data/processed/idconvert.rds'))
interlist=setdiff(unique(idconvert$CMap[!is.na(idconvert$DrugAge)]),NA)

pdf('./results/result_compare_array_gtex.pdf',width = 3.5,useDingbats = F,height = 6)
ggplot(mymat,aes(x=type,y=drname,label=round(mean,2),fill=mean,color=p))+
  geom_tile()+
  geom_text(aes(size=p),fontface='bold')+
  scale_color_manual(values = c('gray40','black'),guide=F)+
  scale_fill_gradient2(midpoint = 0,guide = guide_colorbar('Mean\nSimilarity'))+
  scale_size_manual(values=c(1,1.5),guide=F)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text = element_text(size=6,face='bold'),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(colour = c('gray25','darkred')[1+levels(mymat$drname)%in%interlist]),
        legend.title = element_text(size=5,face='bold'),
        legend.text = element_text(size=4),
        legend.key.size = unit(5,units = 'pt'))
dev.off()
