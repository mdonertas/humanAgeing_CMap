# rm(list=ls())
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

pdf('./results/compare_methods.pdf',useDingbats = F,width = 10)
library(ggplot2)
res3=unique(res2[,c('Array_mean','Average_lifespanChange','Array_names','Array_p','DrugAgeName','Array_signif')])
res3=res3[which(!is.na(res3$Array_mean) & !is.na(res3$DrugAgeName)),]
res3$Array_names=factor(res3$Array_names,levels=unique(names(sort(setNames(res3$Array_mean,res3$Array_names),dec=T))))
ggplot(res3,aes(y=Array_mean,x=Average_lifespanChange,color=Array_signif))+
  geom_vline(xintercept =0,lty=2,color='gray35')+
  geom_point(alpha=0.5)+
  theme_bw()+
  scale_size_manual(values=c(1,1,3))+
  ggtitle('Microarray vs. DrugAge')+
  guides(color=guide_legend('Significant in Array'))+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))+xlim(-120,120)


res3=unique(res2[,c('GTEx_mean','Average_lifespanChange','GTEx_names','GTEx_p','DrugAgeName','GTEx_signif')])
res3=res3[which(!is.na(res3$GTEx_mean) & !is.na(res3$DrugAgeName)),]
res3$GTEx_names=factor(res3$GTEx_names,levels=unique(names(sort(setNames(res3$GTEx_mean,res3$GTEx_names),dec=T))))
ggplot(res3,aes(y=GTEx_mean,x=Average_lifespanChange,color=GTEx_signif))+
  geom_vline(xintercept =0,lty=2,color='gray35')+
  geom_point(alpha=0.5)+
  theme_bw()+
  scale_size_manual(values=c(1,1,3))+
  ggtitle('GTEx vs. DrugAge')+
  guides(color=guide_legend('Significant in GTEx'))+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))+xlim(-120,120)

res3=unique(res2[,c('Matias_rank','Average_lifespanChange','DrugAgeName')])
n=sum(!is.na(res2$Matias_rank))
res3=res3[complete.cases(res3),]
res3$Matias_names=factor(paste(res3$Matias_rank,res3$DrugAgeName,sep='_'),levels=unique(names(sort(setNames(res3$Matias_rank,paste(res3$Matias_rank,res3$DrugAgeName,sep='_')),dec=F))))
ggplot(res3,aes(y=Matias_rank,x=Average_lifespanChange,color=Matias_rank<(n/20)))+
  geom_vline(xintercept =0,lty=2,color='gray35')+
  geom_point(alpha=0.5)+
  theme_bw()+
  ggtitle('Matias vs. DrugAge')+
  guides(color=guide_legend('in top 5%'))+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))+xlim(-120,120)

res3=unique(res2[,c('Celegans_rank','Average_lifespanChange','DrugAgeName')])
n=sum(!is.na(res2$Celegans_rank))
res3=res3[complete.cases(res3),]
res3$Celegans_names=factor(paste(res3$Celegans_rank,res3$DrugAgeName,sep='_'),levels=unique(names(sort(setNames(res3$Celegans_rank,paste(res3$Celegans_rank,res3$DrugAgeName,sep='_')),dec=F))))
ggplot(res3,aes(y=Celegans_rank,x=Average_lifespanChange,color=Celegans_rank<n/20))+
  geom_vline(xintercept =0,lty=2,color='gray35')+
  geom_point(alpha=0.5)+
  # geom_boxplot(outlier.shape = '.',fill=NA)+
  theme_bw()+
  ggtitle('Matthias - C.elegans vs. DrugAge')+
  guides(color=guide_legend('in top 5%'))+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))+xlim(-120,120)

res3=unique(res2[,c('Drosophila_rank','Average_lifespanChange','DrugAgeName')])
n=sum(!is.na(res2$Drosophila_rank))
res3=res3[complete.cases(res3),]
res3$Drosophila_names=factor(paste(res3$Drosophila_rank,res3$DrugAgeName,sep='_'),levels=unique(names(sort(setNames(res3$Drosophila_rank,paste(res3$Drosophila_rank,res3$DrugAgeName,sep='_')),dec=F))))
ggplot(res3,aes(y=Drosophila_rank,x=Average_lifespanChange,color=Drosophila_rank<n/20))+
  geom_vline(xintercept =0,lty=2,color='gray35')+
  geom_point(alpha=0.5)+
  # geom_boxplot(outlier.shape = '.',fill=NA)+
  theme_bw()+
  ggtitle('Matthias - Drosophila vs. DrugAge')+
  guides(color=guide_legend('in top 5%'))+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))+xlim(-120,120)


res3=unique(res2[,c('Array_mean','GTEx_mean','CMap_names','Array_p','DrugAgeName','CMap_signif')])
ggplot(res3,aes(x=Array_mean,y=GTEx_mean,color=is.na(DrugAgeName),shape=CMap_signif,label=CMap_names,size=CMap_signif))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(1,3,3,3))+
  scale_shape_manual(values=c(16,17,18,15))+
  ggtitle('Microarray vs. GTEx')+
  guides(shape=guide_legend('Significant in Array'),
         size=guide_legend('Significant in Array'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('Array_mean','Matias_rank','Array_names','Array_p','DrugAgeName','Array_signif')])
ggplot(res3,aes(x=Matias_rank,y=Array_mean,color=is.na(DrugAgeName),shape=Array_signif,label=Array_names,size=Array_signif))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('Microarray vs. Matias')+
  guides(shape=guide_legend('Significant in Array'),
         size=guide_legend('Significant in Array'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('GTEx_mean','Matias_rank','GTEx_names','GTEx_p','DrugAgeName','GTEx_signif')])
ggplot(res3,aes(x=Matias_rank,y=GTEx_mean,color=is.na(DrugAgeName),shape=GTEx_signif,label=GTEx_names,size=GTEx_signif))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('GTEx vs. Matias')+
  guides(shape=guide_legend('Significant in GTEx'),
         size=guide_legend('Significant in GTEx'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('GTEx_mean','Drosophila_rank','GTEx_names','GTEx_p','DrugAgeName','GTEx_signif')])
ggplot(res3,aes(x=Drosophila_rank,y=GTEx_mean,color=is.na(DrugAgeName),shape=GTEx_signif,label=GTEx_names,size=GTEx_signif))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('GTEx vs. Matthias - Drosophila')+
  guides(shape=guide_legend('Significant in GTEx'),
         size=guide_legend('Significant in GTEx'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('GTEx_mean','Celegans_rank','GTEx_names','GTEx_p','DrugAgeName','GTEx_signif')])
ggplot(res3,aes(x=Celegans_rank,y=GTEx_mean,color=is.na(DrugAgeName),shape=GTEx_signif,label=GTEx_names,size=GTEx_signif))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('GTEx vs. Matthias - C.elegans')+
  guides(shape=guide_legend('Significant in GTEx'),
         size=guide_legend('Significant in GTEx'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('Array_mean','Drosophila_rank','Array_names','Array_p','DrugAgeName','Array_signif')])
ggplot(res3,aes(x=Drosophila_rank,y=Array_mean,color=is.na(DrugAgeName),shape=Array_signif,label=Array_names,size=Array_signif))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('Microarray vs. Matthias - Drosophila')+
  guides(shape=guide_legend('Significant in Array'),
         size=guide_legend('Significant in Array'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('Array_mean','Celegans_rank','Array_names','Array_p','DrugAgeName','Array_signif')])
ggplot(res3,aes(x=Celegans_rank,y=Array_mean,color=is.na(DrugAgeName),shape=Array_signif,label=Array_names,size=Array_signif))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('Microarray vs. Matthias - Celegans')+
  guides(shape=guide_legend('Significant in Array'),
         size=guide_legend('Significant in Array'),
         color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('Drosophila_rank','Celegans_rank','DrugAgeName')])
ggplot(res3,aes(x=Celegans_rank,y=Drosophila_rank,color=is.na(DrugAgeName),label=DrugAgeName))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('Matthias - Drosophila vs. C.elegans')+
  guides(color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('Drosophila_rank','Matias_rank','DrugAgeName')])
ggplot(res3,aes(x=Matias_rank,y=Drosophila_rank,color=is.na(DrugAgeName),label=DrugAgeName))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('Matias vs. Matthias - Drosophila')+
  guides(color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

res3=unique(res2[,c('Celegans_rank','Matias_rank','DrugAgeName')])
ggplot(res3,aes(x=Matias_rank,y=Celegans_rank,color=is.na(DrugAgeName),label=DrugAgeName))+
  geom_point(alpha=0.3)+
  geom_text()+
  theme_bw()+
  scale_size_manual(values=c(3,3,5))+
  ggtitle('Matias vs. Matthias - C.elegans')+
  guides(color=guide_legend('DrugAge'))+
  scale_color_brewer(type = 'qual',palette='Set1',labels=c('in DrugAge','-'))

dev.off()
