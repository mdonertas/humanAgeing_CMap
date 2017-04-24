rm(list=ls())
library(data.table)
library(pheatmap)
res=read.csv('data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
psx=as.numeric(as.character(res$p))
# psx[is.na(psx)]=1
res$padj=p.adjust(psx,method='fdr')
br_cmap=res[which(res$padj<0.05),]
br_up=readRDS('./data/processed/humanBrainMicroarray/all_up_affy.rds')
br_down=readRDS('./data/processed/humanBrainMicroarray/all_down_affy.rds')
drugmat=fread('./data/raw/CMap/rankMatrix.txt')
inst=read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
plist=drugmat$V1
drugmat$V1=NULL
colnames(drugmat)=as.character(inst$instance_id)
rownames(drugmat)=plist
detailedRes=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_detailed.csv')
head(detailedRes)

updownlist=sample(union(br_down,br_up))
br_cmap_res=inst$instance_id[inst$cmap_name%in%br_cmap$cmap.name]
# br_cmap_res=inst$instance_id[inst$cmap_name%in%c('sirolimus','geldanamycin')]
myhits=as.matrix(drugmat[plist%in%updownlist,which(colnames(drugmat)%in%unique(as.character(br_cmap_res))),with=F])
rownames(myhits)=plist[plist%in%updownlist]



instids=colnames(drugmat)[which(colnames(drugmat)%in%unique(as.character(br_cmap_res)))]
annotcol=data.frame(drug=as.character(unname(sapply(instids,function(id){
  inst$cmap_name[inst$instance_id==id]
}))),
celltype=as.character(unname(sapply(instids,function(id){
  inst$cell2[inst$instance_id==id]
}))),score=sign(unname(sapply(instids,function(id){
  detailedRes$score[detailedRes$instance_id==id]
}))))
rownames(annotcol)=colnames(myhits)
drug=setNames(c(brewer.pal(8,'Accent'),brewer.pal(8,'Set3'))[1:length(unique(annotcol$drug))],unique(annotcol$drug))
celltype=setNames(c(brewer.pal(8,'Pastel2'))[1:length(unique(annotcol$celltype))],unique(annotcol$celltype))
annotrow=data.frame(ExpChange=c('down','up')[1+plist[plist%in%updownlist]%in%br_up])
rownames(annotrow)=rownames(myhits)
ExpChange=setNames((c('lightblue','pink')),c('down','up'))
anno_cols=list(ExpChange=ExpChange,drug=drug,celltype=celltype)
myhits=myhits[,names(sort(sapply(instids,function(id){
  detailedRes$score[detailedRes$instance_id==id]
})))]

pheatmap((myhits),cluster_cols = T,border_color = NA,
         show_rownames = F,cluster_rows=T,
         annotation_row = annotrow,
         annotation_colors = anno_cols,
         annotation_col = annotcol,
         show_colnames = F,
         annotation_names_row = F,
         color=rev(brewer.pal(5,'RdBu')),
         breaks = seq(0,nrow(drugmat),length.out = 6))

pc=prcomp((as.matrix(drugmat[,which(colnames(drugmat)%in%unique(as.character(br_cmap_res))),with=F])),scale=T)

colx=setNames(rep(1,length(plist)),plist)
colx[br_up]=2
colx[br_down]=3
plot(x=pc$x[,1],y=pc$x[,8],pch=19,col=c(rgb(1,1,1,0),rgb(1,0,0,1),rgb(0,0,1,1))[colx])

km=kmeans(pc$x,2)
mean(br_up%in%plist[km$cluster==1])
mean(br_down%in%plist[km$cluster==1])
