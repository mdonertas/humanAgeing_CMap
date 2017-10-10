br_up=as.character(read.table('./data/processed/GTEX/all_up.grp')[,1])
br_down=as.character(read.table('./data/processed/GTEX/all_down.grp')[,1])
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
colnames(genesthatmatter)=colnames(drugmat)

res=read.csv('data/processed/GTEx//CMap_results/nobrain_consistent_byname.csv')
psx=as.numeric(as.character(res$p))
# psx[is.na(psx)]=1
res$padj=p.adjust(psx,method='fdr')
br_cmap=res[which(res$padj<0.05),]
br_cmap_res=inst$instance_id[inst$cmap_name%in%br_cmap$cmap.name]
# br_cmap_res=inst$instance_id[inst$cmap_name%in%c('sirolimus','geldanamycin')]

myhits=genesthatmatter[,as.character(br_cmap_res)]

myhits=myhits[,as.character(detailedRes$instance_id[detailedRes$instance_id%in%colnames(myhits)])]
instids=colnames(myhits)
annotcol=data.frame(drug=as.character(unname(sapply(instids,function(id){
  inst$cmap_name[inst$instance_id==id]
}))),score=c('opposite','-','similar')[sign(unname(sapply(instids,function(id){
  detailedRes$score[detailedRes$instance_id==id]
})))+2])
rownames(annotcol)=colnames(myhits)
drug=setNames(colorRampPalette(c(brewer.pal(8,'Accent'),brewer.pal(8,'Set3')))(length(unique(annotcol$drug))),unique(annotcol$drug))
celltype=setNames(c(brewer.pal(8,'Pastel2'))[1:length(unique(annotcol$celltype))],unique(annotcol$celltype))
score=setNames(c('firebrick','white','midnightblue'),unique(annotcol$score))
updownlist=c(br_up,br_down)
annotrow=data.frame(ExpChange=c('down','up')[1+rownames(myhits)%in%br_up])
rownames(annotrow)=rownames(myhits)
ExpChange=setNames((c('lightblue','pink')),c('down','up'))
anno_cols=list(ExpChange=ExpChange,drug=drug,celltype=celltype,score=score)

ph=pheatmap((myhits),cluster_cols = T,
         show_rownames = F,cluster_rows=T,
         annotation_row = annotrow,
         annotation_colors = anno_cols,
         annotation_col = annotcol,
         show_colnames = F,
         annotation_names_row = F,
         color=rev(brewer.pal(3,'RdBu')),cutree_rows = 10,cutree_cols = 9)
hc=ph$tree_row
i=9
pheatmap((myhits[cutree(hc,10)==i,]),cluster_cols = T,
         show_rownames = F,cluster_rows=T,
         annotation_row = annotrow,
         annotation_colors = anno_cols,
         annotation_col = annotcol,
         show_colnames = F,
         annotation_names_row = F,
         color=rev(brewer.pal(3,'RdBu')),cutree_cols = 3)


att2=getBM(attributes = c('ensembl_gene_id','hgnc_symbol','description'),filters = 'affy_hg_u133a',values = names(which(cutree(hc,10)==i)),mart = useMart('ensembl','hsapiens_gene_ensembl'))
