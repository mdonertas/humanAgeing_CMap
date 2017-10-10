filesx=grep('age.rds',list.files('./data/processed/expression/log2qn/'),v=T)
# filesx=grep('Berchtold2008_SFG|Colantuoni2011_PFC|Kang2011_DFC|Lu2004_FC|Maycox2009_APFC|Somel2010_PFC',filesx,v=T)
pathx='./data/processed/expression/log2qn/'
expdata=list()
for(fx in filesx){
  nm=gsub('_age.rds','',fx)
  expdata[[nm]]=readRDS(paste(pathx,fx,sep=''))
}
gnames=Reduce('intersect',sapply(expdata,rownames))
dnames=Reduce('union',sapply(expdata,colnames))
mymat=matrix(NA,nrow=length(gnames),ncol=length(dnames),dimnames=list(gnames,dnames))
for(expr in expdata){
  mymat[gnames,colnames(expr)]=expr[gnames,]
}
source('../shared/functions/functions.R')
mymat_qn=quantileNorm(mymat)
filesx=grep('age.rds',list.files('./data/processed/ages/'),v=T)
# filesx=grep('Berchtold2008_SFG|Colantuoni2011_PFC|Kang2011_DFC|Lu2004_FC|Maycox2009_APFC|Somel2010_PFC',filesx,v=T)
pathx='./data/processed/ages/'
agedat=list()
for(fx in filesx){
  nm=gsub('_age.rds','',fx)
  agedat[[nm]]=readRDS(paste(pathx,fx,sep=''))
}
agedat$Colantuoni2011_PFC=agedat$Colantuoni2011_PFC*365

xx=sapply(agedat,function(x)rank(x,ties.method = 'first'))

# dnames=c(sapply(names(xx),function(nm)tapply(xx[[nm]],cut(xx[[nm]],11),function(x)sample(names(x),1))))

agex=setNames(rep(NA,ncol(mymat_qn)),colnames(mymat_qn))
for(agx in agedat){
  agex[names(agx)]=agx
}
mymat=mymat[,dnames]
agex=agex[dnames]
mymat_qn=quantileNorm(mymat)

library(reshape2)
xx=melt(sapply(expdata,colnames))
datasets=setNames(xx[,2],xx[,1])
pcx=pca(mymat,labelx = round(agex/365))
pcx=data.frame(pcx$x)
pcx$age=agex[rownames(pcx)]
pcx$dataset=datasets[rownames(pcx)]
pcx$datasource=sapply(strsplit(datasets[rownames(pcx)],'_'),function(x)x[1])
library(ggplot2)
ggplot(pcx,aes(x=PC1,y=PC2,color=datasource))+geom_point()+theme_bw()
library(sva)
age_025=agex^0.25
mod1=model.matrix(~age_025)
mod0=cbind(mod1[,1])
svs=sva(dat=mymat_qn,mod = mod1,mod0=mod0)
exprx_qn_sv=t(apply(mymat_qn,1,function(x)lm(x~svs$sv)$residuals))
exprx_qn_sv_qn=quantileNorm(exprx_qn_sv)
pcx2=pca(exprx_qn_sv_qn)
pcx2=data.frame(pcx2$x)
pcx2$age=agex[rownames(pcx2)]
pcx2$dataset=datasets[rownames(pcx2)]
pcx2$datasource=sapply(strsplit(datasets[rownames(pcx2)],'_'),function(x)x[1])
library(ggplot2)
ggplot(pcx2,aes(x=PC1,y=PC2,color=datasource))+geom_point()+theme_bw()
ggplot(pcx2,aes(x=PC1,y=PC2,color=age))+geom_point()+theme_bw()

mynewdat=exprx_qn_sv_qn

resx=data.frame(t(apply(mynewdat,1,function(x){
  lmx=summary(lm(x~age_025))
  cox=cor.test(x,age_025,method='s')
  setNames(c(lmx$coefficients['age_025','Estimate'],lmx$coefficients['age_025',4],cox$estimate,cox$p.value),c('beta','lm_p','rho','spearman_p'))
})))
resx$lm_p.adj=p.adjust(resx$lm_p,method='BY')
resx$spear_p.adj=p.adjust(resx$spearman_p,method='BY')
sum(resx$lm_p.adj<0.05)
ggplot(resx,aes(x=beta,y=-log10(lm_p.adj)))+geom_point()+theme_bw()

xx=data.frame(exprx=mynewdat[setdiff(oldgenes_down,names(mynewgenes_down))[3],],age=agex)
xx$age025=xx$age^0.25
xx$dsource=pcx$datasource
ggplot(xx,aes(x=age025,y=exprx))+geom_smooth(method='lm')+geom_point(size=0.3,aes(color=dsource))+theme_bw()+facet_wrap(~dsource)
ggplot(xx,aes(x=age025,y=exprx))+geom_smooth(method='lm')+geom_point(size=1,aes(color=dsource))+theme_bw()

oldgenes_up=read.table('./data/processed/humanBrainMicroarray_LuIncluded/ensembl_consistent_up.txt')[,1]
oldgenes_down=read.table('./data/processed/humanBrainMicroarray_LuIncluded/ensembl_consistent_down.txt')[,1]

resx$inold_up=rownames(resx)%in%oldgenes_up
resx$inold_down=rownames(resx)%in%oldgenes_down
resx$inold=factor(resx$inold_down|resx$inold_up,levels=c(T,F))
ggplot(resx,aes(x=beta,y=-log10(lm_p.adj),color=inold))+geom_point()+theme_bw()

newgenes_up=names(sort(setNames((resx$beta),rownames(resx)),dec=T)[1:500])
newgenes_down=names(sort(setNames((resx$beta),rownames(resx)),dec=F)[1:500])
sum(newgenes_up%in%oldgenes_up)
sum(newgenes_down%in%oldgenes_down)

glist=setNames(rep(0,nrow(mynewdat)),rownames(mynewdat))
glist[newgenes_up]=1
go_up=go_enrich.test(genelist = glist,selection = function(x)x==1,nodesize = 10)
glist=setNames(rep(0,nrow(mynewdat)),rownames(mynewdat))
glist[newgenes_down]=1
go_down=go_enrich.test(genelist = glist,selection = function(x)x==1,nodesize = 10)

go_up[go_up$p.adjusted<0.05,2]
go_down[go_down$p.adjusted<0.05,2]

xx=setNames(resx[resx$spear_p.adj<0.05,'rho'],rownames(resx)[resx$spear_p.adj<0.05])
mynewgenes_up=sort(xx[xx>0],dec=T)
mynewgenes_down=sort(xx[xx<0],dec=F)

library(biomaRt)
affy_up=getBM(attributes = c('affy_hg_u133a','ensembl_gene_id'),filters = 'ensembl_gene_id',values = names(mynewgenes_up),mart = useMart('ensembl','hsapiens_gene_ensembl'))
affy_down=getBM(attributes = c('affy_hg_u133a','ensembl_gene_id'),filters = 'ensembl_gene_id',values = names(mynewgenes_down),mart = useMart('ensembl','hsapiens_gene_ensembl'))

affy_up=setNames(affy_up[,1],affy_up[,2])
affy_up=affy_up[!affy_up=='']

affy_up=unname(unlist(sapply(newgenes_up,function(g)affy_up[names(affy_up)%in%g]))[1:500])

affy_down=setNames(affy_down[,1],affy_down[,2])
affy_down=affy_down[!affy_down=='']
affy_down=unname(unlist(sapply(newgenes_down,function(g)affy_down[names(affy_down)%in%g]))[1:500])

oldaffy_up=as.character(read.table('./data/processed/humanBrainMicroarray_LuIncluded/consistent_up.grp')[,1])
oldaffy_down=as.character(read.table('./data/processed/humanBrainMicroarray_LuIncluded/consistent_down.grp')[,1])

sum(affy_down%in%oldaffy_down)
sum(affy_up%in%oldaffy_up)

write.table(affy_up,file='~/Desktop/combined_up.grp',quote = F,row.names = F,col.names = F)
write.table(affy_down,file='~/Desktop/combined_down.grp',quote = F,row.names = F,col.names = F)

combined_res=read.csv('~/Desktop/combined_all.csv')
head(combined_res)
old_res=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')

newres=setNames(combined_res$mean,combined_res$cmap.name)
oldres=setNames(old_res$mean,old_res$cmap.name)
# oldres=oldres[names(newres)]
newres=newres[names(oldres)]
cor.test(oldres,newres,method='s')
xx=data.frame(old=oldres,new=newres)
head(xx)
head(xx)
xx$signifold=rownames(xx)%in%names(newres[as.character(old_res$cmap.name[which(old_res$padj<0.05)])])
ggplot(xx,aes(x=old,y=new,color=signifold))+geom_point()+theme_bw()

xx=list(olddown=oldgenes_down,oldup=oldgenes_up,newdown=names(mynewgenes_down),newup=names(mynewgenes_up))
library(VennDiagram)
venn.diagram(xx,filename='~/Desktop/venn.tiff',imagetype = 'tiff')

combined_res=combined_res[which(combined_res$mean<0),]

