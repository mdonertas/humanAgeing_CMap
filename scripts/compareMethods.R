library(reshape2)
idconvert=readRDS('../shared/chembl/CMap.rds')
melike=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
ps=as.numeric(as.character(melike$p))
melike$padj=p.adjust(ps,method='fdr')
matias=read.csv('./docs/matias_list/all.csv')

ranksinmatias=sapply(melike$cmap.name,function(drugname){
  x=which(matias$drug%in%setdiff(unique(idconvert$chembl[idconvert$Drug%in%drugname]),NA))
  if(length(x)<1){NA}
  else(x)
  })
names(ranksinmatias)=melike$cmap.name
xx=melt(ranksinmatias)
matias=setNames(xx$value,xx$L1)

drosophila=read.csv('./docs/matthias_list/drosophila.csv')

drosophilarank=sapply(as.character(melike$cmap.name),function(drugname){
  x=which(drosophila$Drug..HetCode.%in%setdiff(unique(unlist(idconvert$het[idconvert$Drug%in%drugname])),NA))
  if(length(x)<1){NA}
  else(x)
})
names(drosophilarank)=melike$cmap.name
xx=melt(drosophilarank)
drosophilarank=setNames(xx$value,rownames(xx))

celegans=read.csv('./docs/matthias_list/celegans.csv')

celegansrank=sapply(as.character(melike$cmap.name),function(drugname){
  x=which(celegans$Drug..HetCode.%in%setdiff(unique(unlist(idconvert$het[idconvert$Drug%in%drugname])),NA))
  if(length(x)<1){NA}
  else(x)
})
names(celegansrank)=melike$cmap.name
xx=melt(celegansrank)
celegansrank=setNames(xx$value,rownames(xx))

melike_mean=setNames(melike$mean,melike$cmap.name)
melike_p=setNames(melike$padj,melike$cmap.name)

drugage=readRDS('../shared/chembl/DrugAge.rds')





xx=lapply(melike$cmap.name,function(nm){
  data.frame(mean=subset(melike,cmap.name==nm,select = mean),
             p=as.character(c('NS','FDR<=0.05')[1+(subset(melike,cmap.name==nm,select = padj)<=0.05)]),
             matias=matias[names(matias)%in%nm],
             drosophila=drosophilarank[names(drosophilarank)%in%nm],
             celegans=celegansrank[names(celegansrank)%in%nm],
             drugage=any(setdiff(unique(unname(unlist(idconvert[idconvert$Drug%in%nm,]))),NA)%in%setdiff(unique(unname(unlist(drugage))),NA)),
             melike_rank=subset(melike,cmap.name==nm,select = rank))
})
names(xx)=as.character(melike$cmap.name)
xx=melt(xx,id.vars=c('p','mean','matias','drosophila','celegans','drugage','rank'))
head(xx)
colnames(xx)=c('Significant','Melike','Matias','Matthias-Drosophila','Matthias-Celegans','DrugAge','Rank-Melike','Drugname')

mydat=xx
mydat$Significant=as.character(mydat$Significant)
mydat$Significant[is.na(mydat$Significant)]='-'
mydat$Significant=factor(mydat$Significant,levels=c('-','NS','FDR<=0.05'))
mydat$DrugAge=factor(c('-','in DrugAge')[1+mydat$DrugAge],levels=c('in DrugAge','-'))
library(ggplot2)
ggplot(mydat,aes(x=Matias,y=Melike,color=DrugAge))+
  geom_point(aes(size=Significant))+
  scale_size_manual(values = c(1,2,4))+
  theme_bw()+
  scale_color_brewer(type='qual',palette='Set1')
head(mydat)

head(mydat)

mydat$names=mydat$Drugname
mydat$names[mydat$DrugAge!='in DrugAge']=''


ggplot(mydat,aes(x=`Rank-Melike`, y=Melike, color=DrugAge))+
  geom_point(aes(size=Significant),alpha=0.5)+
  geom_text(aes(label=names),color='gray25',alpha=0.75,nudge_y = sample(c(0.01,-0.01),rep=T,nrow(mydat)))+
  scale_size_manual(values = c(1,1,4))+
  theme_bw()+
  scale_color_brewer(type='qual',palette='Set1')


