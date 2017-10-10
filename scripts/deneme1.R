rm(list=ls())
#################
library(RColorBrewer)
library(data.table)
library(pheatmap)
library(biomaRt)
#################
br_up=readRDS('./data/processed/humanBrainMicroarray/all_up_affy.rds')
br_down=readRDS('./data/processed/humanBrainMicroarray/all_down_affy.rds')
drugmat=fread('./data/raw/CMap/rankMatrix.txt')
inst=read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
plist=drugmat$V1
drugmat$V1=NULL
colnames(drugmat)=as.character(inst$instance_id)
rownames(drugmat)=plist
detailedRes=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_detailed.csv')
#################
head(detailedRes)

k.calc=function(drugmat,i,myprobelist){
  ks=c()
  glist=melt(myprobelist)
  glist$res=0
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
    if(a>b){
      k=a
      mylim=which.max(x1-x2)
    } else {
        k= -b
        mylim=which.max(x2-x3)
        }
    ks[nm]=k
    if(nm == 'up'){
      genes=v[mylim:length(v)]
    }
    if(nm == 'down'){
      genes=v[1:mylim]
    }
    glist$res[glist$value%in%names(genes)]=1
  }
  ks['score']=ks['up']-ks['down']
  # ks=data.frame(ks,up=glist[['up']],down=glist[['down']])
  res=list(ks,setNames(glist$res,glist$value))
  names(res)=c('result','genelist')
  return(res)
}
# k.calc(drugmat,1,list(up=br_up,down=br_down))


allres=lapply(inst$instance_id,function(ins)k.calc(drugmat,ins,list(up=br_up,down=br_down)))
res=as.data.frame(t(sapply(allres,function(res)res$result)))
rownames(res)=inst$instance_id
maxs=max(res$score)
mins=min(res$score)

res$finscore=sapply(res$score,function(x){
  if(x>=0){x/maxs}
  else{-x/mins}
})

genesthatmatter=sapply(allres,function(res)res$genelist[c(br_up,br_down)])

ups=sapply(allres,function(res)res$genelist$up[br_up])
colnames(ups)=inst$instance_id

downs=sapply(allres,function(res)res$genelist$down[br_down])
colnames(downs)=inst$instance_id

res=read.csv('data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
psx=as.numeric(as.character(res$p))
# psx[is.na(psx)]=1
res$padj=p.adjust(psx,method='fdr')
br_cmap=res[which(res$padj<0.05),]
# br_cmap_res=inst$instance_id[inst$cmap_name%in%br_cmap$cmap.name]
br_cmap_res=inst$instance_id[inst$cmap_name%in%c('sirolimus','geldanamycin')]

myhits=drugmat[,as.character(br_cmap_res),with=F]
# myhits=downs
instids=colnames(drugmat)[which(colnames(drugmat)%in%unique(as.character(br_cmap_res)))]
annotcol=data.frame(drug=c('-','p<0.05')[1+colnames(myhits)%in%instids],
                    score=sign(unname(sapply(colnames(myhits),function(id){
                      detailedRes$score[detailedRes$instance_id==as.numeric(id)]
                    }))))
rownames(annotcol)=colnames(myhits)
drug=setNames(c(brewer.pal(8,'Accent'),brewer.pal(8,'Set3'))[1:length(unique(annotcol$drug))],unique(annotcol$drug))
updownlist=c(br_down)
annotrow=data.frame(ExpChange=c('down','up')[1+plist[plist%in%updownlist]%in%br_up])
rownames(annotrow)=rownames(myhits)
ExpChange=setNames((c('lightblue','pink')),c('down','up'))
anno_cols=list(ExpChange=ExpChange,drug=drug)
myhits=myhits[,names(sort(sapply(instids,function(id){
  detailedRes$score[detailedRes$instance_id==id]
})))]


myhits=myhits[,names(sort(setNames(annotcol[colnames(myhits),2],colnames(myhits))))]
# pdf('~/Desktop/deneme.pdf',width = 50)
pheatmap(myhits,cluster_cols = T,
         show_rownames = F,cluster_rows=T,
         annotation_row = annotrow,
         annotation_colors = anno_cols,
         annotation_col = annotcol,
         show_colnames = F,
         annotation_names_row = F,
         color=rev(brewer.pal(5,'RdBu')),
         breaks = seq(0,1,length.out = 6),cutree_rows = 10,cutree_cols = 5)
# dev.off()

mygrp=as.factor(unname(sort(setNames(as.character(unname(sapply(colnames(myhits),function(x)inst$cmap_name[inst$instance_id==x]))),colnames(myhits)))))

tres=apply(myhits,1,function(x)t.test(x~mygrp)$p.val)




att=getBM(attributes = c('ensembl_gene_id','hgnc_symbol','description'),filters = 'affy_hg_u133a',values = names(which(p.adjust(tres,method='fdr')<0.05)),mart = useMart('ensembl','hsapiens_gene_ensembl'))
