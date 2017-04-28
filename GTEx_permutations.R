permfilex=list.files('../GTEx/data/processed/permutations/permres/',full.names = T)
library(data.table)
library(biomaRt)
library(parallel)
k.calc=function(drugmat,i,myprobelist){
  ks=c()
  for( nm in names(myprobelist)){
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

drugmat=fread('./data/raw/CMap/rankMatrix.txt')
inst=read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
plist=drugmat$V1
drugmat$V1=NULL
colnames(drugmat)=as.character(inst$instance_id)
rownames(drugmat)=plist

real_up=readRDS('./data/processed/GTEx/all_up.rds')
real_down=readRDS('./data/processed/GTEx/all_down.rds')
real_mat=readRDS('./data/processed/GTEx/cormat.rds')
humanmart=useMart('ensembl','hsapiens_gene_ensembl')
idconv=getBM(attributes = c('affy_hg_u133a','ensembl_gene_id'),filters = 'ensembl_gene_id',values = rownames(real_mat),mart = humanmart)
idconv=setNames(idconv$affy_hg_u133a,idconv$ensembl_gene_id)

allres=mclapply(permfilex,function(permfile){
  load(permfile)
  cormat=cormat[rownames(real_mat),intersect(colnames(cormat),colnames(real_mat))]
  perm_up=names(sort(rowSums(cormat>0),dec=T)[1:length(real_up)])
  perm_down=names(sort(rowSums(cormat<0),dec=T)[1:length(real_down)])
  perm_up_affy=idconv[names(idconv)%in%perm_up]
  perm_down_affy=idconv[names(idconv)%in%perm_down]
  int=union(intersect(perm_down_affy,perm_up_affy),'')
  perm_down_affy=setdiff(perm_down_affy,int)
  perm_up_affy=setdiff(perm_up_affy,int)
  res=as.data.frame(t(sapply(inst$instance_id,function(ins)k.calc(drugmat,ins,list(up=perm_up_affy,down=perm_down_affy)))))
  rownames(res)=inst$instance_id
  maxs=max(res$score)
  mins=min(res$score)
  res$finscore=sapply(res$score,function(x){
    if(x>=0){x/maxs}
    else{-x/mins}
  })
  res
},mc.cores = 8)
saveRDS(allres,file='./data/processed/GTEx_all_CMap_permres.rds')

scoredist=sapply(allres,function(x)x$finscore)

rownames(scoredist)=rownames(allres[[1]])

real_up_affy=as.character(read.table('./data/processed/GTEx/all_up.grp')[,1])
real_down_affy=as.character(read.table('./data/processed/GTEx/all_down.grp')[,1])
realres=as.data.frame(t(sapply(inst$instance_id,function(ins)k.calc(drugmat,ins,list(up=real_up_affy,down=real_down_affy)))))
rownames(realres)=inst$instance_id
maxs=max(realres$score)
mins=min(realres$score)
realres$finscore=sapply(realres$score,function(x){
  if(x>=0){x/maxs}
  else{-x/mins}
})
head(realres)

pvals=sapply(rownames(realres),function(inst){
  mean(abs(scoredist[inst,])>=abs(realres[inst,4]))
})

padjvals=p.adjust(pvals,method='fdr')

head(inst)
rownames(inst)=inst$instance_id
unique(inst[names(which(padjvals<0.05)),3])
# nitrendipine 
# coralyne     
# butamben     
# deptropine   
# alprostadil  
# triprolidine