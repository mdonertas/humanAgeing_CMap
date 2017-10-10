myfunc=function(i){
  print(i)
  library(data.table)
  #################
  load('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/brainGTEx/all.RData')
  br_up=up_affy3
  br_down=down_affy3
  rm(list=setdiff(ls(),c('br_up','br_down')))
  drugmat=fread('/nfs/research2/thornton/melike/humanAgeing_CMap/data/raw/CMap/rankMatrix.txt')
  inst=read.csv('/nfs/research2/thornton/melike/humanAgeing_CMap/data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
  plist=drugmat$V1
  drugmat$V1=NULL
  colnames(drugmat)=as.character(inst$instance_id)
  rownames(drugmat)=plist
  significantres=read.csv('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/brainGTEx/CMap_results/significantres.csv')
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
  res$instance_id=inst$instance_id
  res=setNames(res$finscore,res$instance_id)
  return(res)
}
# allres=Q(myfunc,i=1:1000, n_jobs=1000)
# saveRDS(allres,file='./data/processed/brainGTEx/1000times500.rds')
allres=readRDS('./data/processed/brainGTEx/1000times500.rds')
myfunc2=function(i){
  # print(i)
  library(data.table)
  #################
  load('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/brainGTEx/all.RData')
  br_up=up_affy3
  br_down=down_affy3
  rm(list=setdiff(ls(),c('br_up','br_down')))
  drugmat=fread('/nfs/research2/thornton/melike/humanAgeing_CMap/data/raw/CMap/rankMatrix.txt')
  inst=read.csv('/nfs/research2/thornton/melike/humanAgeing_CMap/data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
  plist=drugmat$V1
  drugmat$V1=NULL
  colnames(drugmat)=as.character(inst$instance_id)
  rownames(drugmat)=plist
  significantres=read.csv('/nfs/research2/thornton/melike/humanAgeing_CMap/data/processed/brainGTEx/CMap_results/significantres.csv')
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
  res=as.data.frame(t(sapply(inst$instance_id,function(ins)k.calc(drugmat,ins,list(up=setdiff(unique(names(sort(br_up,dec=T))),names(br_down))[1:500],down=setdiff(unique(names(sort(br_down,dec=F))),names(br_up))[1:500])))))
  maxs=max(res$score)
  mins=min(res$score)
  res$finscore=sapply(res$score,function(x){
    if(x>=0){x/maxs}
    else{-x/mins}
  })
  res$instance_id=inst$instance_id
  res=setNames(res$finscore,res$instance_id)
  return(res)
}
realgtex=Q(myfunc2,i=1,n_jobs = 1)
realgtex=realgtex[[1]]

resx=sapply(allres,function(x)x)
corx=cor(resx,method='s')
gtex=read.csv('./data/processed/brainGTEx/CMap_results/gtex_brain_detailed.csv')

gtex=setNames(gtex$score,gtex$instance_id)
gtex=gtex[names(allres[[1]])]

realgtexcorx=apply(resx,2,function(x)cor(x,gtex,method='s'))

par(mfrow=c(2,1))
hist(corx[upper.tri(corx)],br=25,xlim=c(0.75,0.89))
hist(realgtexcorx,br=25,xlim=c(0.75,0.89))

summary(realgtexcorx)
