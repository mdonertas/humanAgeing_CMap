allres=lapply(inst$instance_id,function(ins)k.calc(drugmat,ins,list(up=br_up,down=br_down)))
res=as.data.frame(t(sapply(allres,function(res)res$result)))
rownames(res)=inst$instance_id
maxs=max(res$up)
mins=min(res$up)

res$finup=sapply(res$up,function(x){
  if(x>=0){x/maxs}
  else{-x/mins}
})

maxs=max(res$down)
mins=min(res$down)

res$findown=sapply(res$down,function(x){
  if(x>=0){x/maxs}
  else{-x/mins}
})

maxs=max(res$score,na.rm=T)
mins=min(res$score,na.rm=T)

res$finscore=sapply(res$score,function(x){
  if(is.na(x)==F){if(x>=0){x/maxs}
  else{-x/mins}}else{NA}
})

head(res)

finres=as.data.frame(t(sapply(unique(as.character(inst$cmap_name)),function(x){
  mymat=res[as.character(inst$instance_id[as.character(inst$cmap_name)==x]),]
  colMeans(mymat)
})))

finres$calcp=sapply(unique(as.character(inst$cmap_name)),function(x){
  mymat=res[as.character(inst$instance_id[as.character(inst$cmap_name)==x]),]
  mean(mymat[,1]*mymat[,2]<0)
})>0

doperm=function(drugmat,inst,br_up,br_down){
  mup=sample(rownames(drugmat),length(br_up))
  mdown=sample(setdiff(rownames(drugmat),mup),length(br_down))
  allres=lapply(inst$instance_id,function(ins)k.calc(drugmat,ins,list(up=mup,down=mdown)))
  res=as.data.frame(t(sapply(allres,function(res)res$result)))
  rownames(res)=inst$instance_id
  maxs=max(res$up)
  mins=min(res$up)
  res$finup=sapply(res$up,function(x){
    if(x>=0){x/maxs}
    else{-x/mins}
  })
  maxs=max(res$down)
  mins=min(res$down)
  res$findown=sapply(res$down,function(x){
    if(x>=0){x/maxs}
    else{-x/mins}
  })
  maxs=max(res$score)
  mins=min(res$score)
  res$finscore=sapply(res$score,function(x){
    if(x>=0){x/maxs}
    else{-x/mins}
  })
  finres=t(sapply(unique(as.character(inst$cmap_name)),function(x){
    colMeans(res[as.character(inst$instance_id[as.character(inst$cmap_name)==x]),])
  }))
  return(finres)
}

library(parallel)
permres=mclapply(1:1000,function(i)doperm(drugmat,inst,br_up,br_down),mc.cores = 8)
saveRDS(permres,file='./data/processed/cmap_permutations2.rds')

head(finres)

pvals=as.data.frame(t(sapply(rownames(finres),function(nm){
  c(mean(abs(finres[nm,'finscore'])<=abs(sapply(permres,function(x)x[nm,'finscore']))),
    mean(abs(finres[nm,'finup'])<=sapply(permres,function(x)x[nm,'finup'])),
    mean(abs(finres[nm,'findown'])<=abs(sapply(permres,function(x)x[nm,'findown']))))
})))
colnames(pvals)=c('total','up','down')
# pvals[finres$calcp==F,]=NA
# pvals$totFDR=p.adjust(pvals$total,method='fdr')
# pvals$upFDR=p.adjust(pvals$up,method='fdr')
# pvals$downFDR=p.adjust(pvals$down,method='fdr')

pvals[pvals[,2]<=0.05 & pvals[,3]>0.05,]

finres[pvals[,3]<=0.05 & finres[,'findown']>0,]

finres[pvals[,2]<=0.05 & finres[,'finup']>0,]

##################

perms=readRDS('./data/processed/cmap_permutations.rds')
res=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_detailed.csv')
maxs=max(res$up)
mins=min(res$up)
res$finup=sapply(res$up,function(x){
  if(x>=0){x/maxs}
  else{-x/mins}
})
maxs=max(res$down)
mins=min(res$down)
res$findown=sapply(res$down,function(x){
  if(x>=0){x/maxs}
  else{-x/mins}
})
realres=t(sapply(unique(res$cmap.name),function(nm){
  xx=res[res$cmap.name==nm,c('finup','findown')]
  aa=names(which(colMeans(xx>0)>=0.9 | colMeans(xx>0)<=0.1))
  resx=colMeans(xx)[aa]
  resx[setdiff(c('finup','findown'),aa)]=NA
  resx
}))

rownames(realres)=unique(res$cmap.name)
realres=realres[rownames(perms[[1]]),]
upp=setNames(p.adjust(sapply(1:nrow(realres),function(i)mean(na.rm=T,(realres[i,'finup'])<=sapply(perms,function(x)(x[i,'finup'])))),method='fdr'),rownames(realres))
downp=setNames(p.adjust(sapply(1:nrow(realres),function(i)mean(na.rm=T,(realres[i,'findown'])<=sapply(perms,function(x)(x[i,'findown'])))),method='fdr'),rownames(realres))
which(upp<0.05)
which(downp<0.05)

plot(realres)

d=realres[,'findown']
u=realres[,'finup']
xxx=d>0 & u>0
d=d[xxx]
u=u[xxx]
plot(d,u)
sort(d*u,dec=T)[1:20]


sort(d,dec=T)[1:10]
