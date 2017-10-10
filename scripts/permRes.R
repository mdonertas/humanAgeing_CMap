path='/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation'
setwd(path)
fx=setdiff(list.files(path),'cormat')
for(permnum in 1:1000){
  print(permnum)
  filesx=paste(fx,'/perm',permnum,'.rds',sep='')
  corss=lapply(filesx,readRDS)
  names(corss)=fx
  genelist=unique(unname(unlist(sapply(corss,rownames))))
  rmat=matrix(NA,nrow=length(genelist),ncol=length(fx),dimnames=list(genelist,fx))
  for(nm in names(corss)){
    corx=corss[[nm]]
    genes=rownames(corx)
    rmat[genes,nm]=corx$rho
  }
  saveRDS(rmat,file=paste('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/permutation/cormat/perm',permnum,'.rds',sep=''))
  rm(list=setdiff(ls(),c('path','fx','permnum')))
}
