setwd('/nfs/research2/thornton/melike/humanAgeing_CMap')
rm(list=ls())

myfunc=function(i){
  print(i)
  library(data.table)
  #################
  load('./data/processed/brainGTEx/all.RData')
  br_up=up_affy3
  br_down=down_affy3
  rm(list=setdiff(ls(),c('br_up','br_down')))
  drugmat=fread('./data/raw/CMap/rankMatrix.txt')
  inst=read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
  plist=drugmat$V1
  drugmat$V1=NULL
  colnames(drugmat)=as.character(inst$instance_id)
  rownames(drugmat)=plist
  significantres=read.csv('./data/processed/brainGTEx/CMap_results/significantres.csv')
  #################
  k.calc=function(drugmat,i,myprobelist){
    ks=c()
    for(nm in names(myprobelist)){
      myprobes=myprobelist[[nm]]
      plist=rownames(drugmat)
      v=sort(setNames(drugmat[plist%in%myprobes,colnames(drugmat)%in%i,with=F][[1]],plist[plist%in%myprobes]))
      n=length(plist)
      t=length(v)
      j=1:t
      x1=j/t
      x2=v/n
      x3=(j-1)/t
      a=max(x1-x2)
      b=max(x2-x3)
      if(a>b){k=a
      } else {k= -b}
      ks[nm]=k
    }
    ks['score']=ks['up']-ks['down']
    return(ks)
  }
  res=as.data.frame(t(sapply(inst$instance_id,function(ins)k.calc(drugmat,ins,list(up=sample(unique(names(br_up)),500),down=sample(unique(names(br_down)),500))))))
  maxs=max(res$score)
  mins=min(res$score)
  res$finscore=sapply(res$score,function(x){
    if(x>=0){x/maxs}
    else{-x/mins}
  })
  rownames(res)=inst$instance_id
  setNames(sapply(1:nrow(significantres),function(i)mean(res[as.character(inst$instance_id[inst$cmap_name%in%(significantres$cmap.name[i])]),4])),significantres$cmap.name)
}
allres=Q(myfunc,i=1:1000, n_jobs=1000)
saveRDS(allres,file='./data/processed/brainGTEx/500of1000.rds')
