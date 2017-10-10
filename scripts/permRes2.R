rm(list=ls())
path1='/nfs/research2/thornton/melike/GTEx/data/permutations/permres/'
path2='/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/cormat/'
filesx=list.files(path1)
permress=list()
for(fx in filesx){
  print(fx)
  f1=paste(path1,fx,sep='')
  f2=gsub('.RData','.rds',paste(path2,fx,sep=''))
  load(f1)
  micro=readRDS(f2)
  nm=c(colnames(cormat),colnames(micro))
  genelist=union(rownames(cormat),rownames(micro))
  rmat=matrix(NA,ncol=length(nm),nrow=length(genelist),dimnames=list(genelist,nm))
  rmat[rownames(cormat),colnames(cormat)]=cormat
  rmat[rownames(micro),colnames(micro)]=micro
  saveRDS(rmat,file=paste('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/gtex_micro/',gsub('.RData','.rds',fx),sep=''))
}

filesx=list.files('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/gtex_micro/',full.names = T)
permres=sapply(filesx,function(fx){
  cormat=readRDS(fx)
  cormat=cormat[,c(grep('Brain',colnames(cormat)[1:48],v=T),colnames(cormat)[49:ncol(cormat)])]
  c(sum(rowMeans(cormat>0)==1,na.rm=T),sum(rowMeans(cormat<0)==1,na.rm=T))
})
saveRDS(permres,file='/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/gtex_micro_shared.rds')
permres_noLu=sapply(filesx,function(fx){
  cormat=readRDS(fx)
  cormat=cormat[,c(grep('Brain',colnames(cormat)[1:48],v=T),colnames(cormat)[49:ncol(cormat)])]
  cormat=cormat[,grep('Lu2004',colnames(cormat),invert=T)]
  c(sum(rowMeans(cormat>0)==1,na.rm=T),sum(rowMeans(cormat<0)==1,na.rm=T))
})
saveRDS(permres_noLu,file='/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/gtex_micro_shared_noLu.rds')
