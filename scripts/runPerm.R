args=commandArgs(trailingOnly = T)
exps=readRDS(args[1])
age=readRDS(args[2])[[as.numeric(args[3])]]
cortab=as.data.frame(t(apply(exps,1,function(x){
  co=cor.test(x,age[names(x)],method='s')
  c(co$est,co$p.val)
})))
colnames(cortab)=c('rho','p')
cortab$padj=p.adjust(cortab$p,method='fdr')
saveRDS(cortab,file=paste(args[4],'perm',args[3],'.rds',sep=''))