library(clustermq)
drlist=unique(as.character(read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)[,3]))
get_cid <- function(nm){
  library(RCurl)
  library(jsonlite)
  nm=gsub(' ','%20',nm)
  if(!grepl('#',nm)){
    name2cid=getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", nm, "/cids/JSON",sep = ""))
    name2cid=fromJSON(name2cid)
    return(name2cid$IdentifierList$CID)
  }else{return(NA)}
}
cmap=Q(get_cid,nm=drlist,n_jobs = 40)
names(cmap)=drlist

drlist=unique(as.character(read.csv('./docs/matias_list/all.csv')[,2]))
getCID_fromCHEMBL=function(cid){
  library(RCurl)
  library(jsonlite)
  cid2chembl=getURL(paste('https://www.ebi.ac.uk/unichem/rest/src_compound_id/',cid,'/1/22',sep=''))
  cid2chembl=fromJSON(cid2chembl)
  return(cid2chembl$src_compound_id)
}
matias=Q(getCID_fromCHEMBL,cid=drlist,n_jobs = 40)
names(matias)=drlist

drlist=unique(as.character(read.csv('./docs/matthias_list/drosophila.csv')[,1]))
getPubChem_from_HET=function(cid){
  library(RCurl)
  library(jsonlite)
  cid2chembl=getURL(paste('https://www.ebi.ac.uk/unichem/rest/src_compound_id/',cid,'/3/22',sep=''))
  cid2chembl=fromJSON(cid2chembl)
  return(cid2chembl$src_compound_id)
}
drosophila=Q(getPubChem_from_HET,cid=drlist,n_jobs = 40)
names(drosophila)=drlist

drlist=unique(as.character(read.csv('./docs/matthias_list/celegans.csv')[,1]))
celegans=Q(getPubChem_from_HET,cid=drlist,n_jobs = 40)
names(celegans)=drlist

drlist=as.character(unique(read.csv('./data/raw/DrugAge/DrugAge.csv')$compound_name))
drugage=Q(get_cid,nm=drlist,n_jobs=40)
names(drugage)=drlist

pubchem=unique(c(unique(unname(unlist(cmap))),
      unique(unname(unlist(matias))),
      unique(unname(unlist(drosophila))),
      unique(unname(unlist(celegans))),
      unique(unname(unlist(drugage)))))

idconvert=lapply(pubchem,function(cid){
  list(CMap=names(cmap)[which(sapply(cmap,function(x)any(x%in%cid)))],
  Matias=names(matias)[which(sapply(matias,function(x)any(x%in%cid)))],
  Drosophila=names(drosophila)[which(sapply(drosophila,function(x)any(x%in%cid)))],
  Celegans=names(celegans)[which(sapply(celegans,function(x)any(x%in%cid)))],
  DrugAge=names(drugage)[which(sapply(drugage,function(x)any(x%in%cid)))])
})
names(idconvert)=pubchem
head(idconvert)
library(reshape2)
idconvert2=t(sapply(idconvert,function(x)sapply(x,function(y)ifelse(length(y)==0,NA,y))))
saveRDS(idconvert2,'./data/processed/idconvert.rds')
