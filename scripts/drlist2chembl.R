library(reshape2)
library(jsonlite)
library(RCurl)

# drlist=as.character(read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)[,3])
# drlist=as.character(read.csv('./data/processed/brainGTEx/CMap_results/significantres.csv')$cmap.name)
# drlist=as.character(read.csv('./data/processed/humanBrainMicroarray/CMap_results/')$cmap.name)

cmap=list(microbrain=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv'),
          microbrain_noLu=read.csv('./data/processed/humanBrainMicroarray/CMap_results/noLu_consistent_byname.csv'),
          gtex=read.csv('./data/processed/GTEx/CMap_results/all_consistent_byname.csv'),
          gtex_nobrain=read.csv('./data/processed/GTEx/CMap_results/nobrain_consistent_byname.csv'),
          allbrain=read.csv('./data/processed/allBrain/CMap_results/allbrain.csv'),
          allbrain_noLu=read.csv('./data/processed/allBrain/CMap_results/allbrain_noLu.csv'),
          gtex_brain=read.csv('./data/processed/brainGTEx/CMap_results/GTEx_brain.csv'))

cmap=lapply(c('microbrain','gtex_nobrain','gtex_brain','allbrain'),function(nm)cmap[[nm]])
names(cmap)=c('microbrain','gtex_nobrain','gtex_brain','allbrain')

cmap_padj=lapply(cmap,function(x){
  ps=as.numeric(as.character(x$p))
  padj=p.adjust(ps,method='fdr')
  cbind(x,padj)
})
druglist=unique(unname(unlist(lapply(cmap_padj,function(x)as.character(x$cmap.name[which(x$padj<0.05)])))))


drlist=(unique(druglist))

drres=sapply(drlist,function(nm){
  xx=getURL(paste('https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact=',nm,sep=''),async = T)
  xx=fromJSON(xx)
  xx$molecules$molecule_chembl_id
})
dr_rest=names(which(sapply(drres,is.null)))
drres2=sapply(dr_rest,function(nm){
  xx=getURL(paste('https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=',nm,sep=''),async = T)
  xx=fromJSON(xx)
  xx$molecules$molecule_chembl_id  
})

for(nm in dr_rest){
  drres[[nm]]=drres2[[nm]]
}

idmap=sapply(unique(names(which(sapply(drres,length)==1))),function(nm)drres[[nm]])



mx=as.data.frame(sapply(cmap_padj,function(x){
  setNames(x$mean,x$cmap.name)[druglist]
}))

druglist=names(sort(rowMeans(mx),dec=T))
mx=mx[druglist,]

px=as.data.frame(sapply(cmap_padj,function(x){
  setNames(x$padj,x$cmap.name)[druglist]
})<0.05)

finres=data.frame(drugname=rownames(mx),
           chembl_id=unname(idmap[rownames(mx)]),
           similarity_gtex_others=mx$gtex_nobrain,
           significance_gtex_others=px$gtex_nobrain,
           similarity_gtex_brain=mx$gtex_brain,
           significance_gtex_brain=px$gtex_brain,
           similarity_microarray_brain=mx$microbrain,
           significance_microarray_brain=px$microbrain,
           similarity_combined_brain=mx$allbrain,
           significance_combined_brain=px$allbrain)

head(finres)
write.table(finres,file='~/Desktop/finalres.tsv',quote = F,row.names = F,col.names = T,sep='\t')


