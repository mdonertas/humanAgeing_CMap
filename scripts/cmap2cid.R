drlist=unique(as.character(read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)[,3]))

get_cid <- function(nm){
  library(RCurl)
  library(jsonlite)
  nm=gsub(' ','%20',nm)
  name2cid=getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", nm, "/cids/JSON",sep = ""))
  name2cid=fromJSON(name2cid)
  return(name2cid$IdentifierList$CID)
}
cids=Q(get_cid,nm=drlist,n_jobs = 40)
names(cids)=drlist

xx=melt(cids[which(!sapply(cids,is.null))])

head(xx)
colnames(xx)=c('CID','Drug')

length(unique(xx$Drug))

saveRDS(xx,file='../shared/chembl/CMap2CID.rds')

cids=readRDS('../shared/chembl/CMap2CID.rds')

cidvec=setNames(cids$CID,cids$Drug)

getChembl_from_CID=function(cid){
  library(RCurl)
  library(jsonlite)
  cid2chembl=getURL(paste('https://www.ebi.ac.uk/unichem/rest/src_compound_id/',cid,'/22/1',sep=''))
  cid2chembl=fromJSON(cid2chembl)
  return(cid2chembl$src_compound_id)
}

library(clustermq)
chembl=Q(getChembl_from_CID,cid=cidvec,n_jobs = 50)
names(chembl)=names(cidvec)
cids$chembl=sapply(1:nrow(cids),function(i){
  ch=chembl[[i]]
  if(is.null(ch)){NA}
  else(ch)
})

getHET_from_CID=function(cid){
  library(RCurl)
  library(jsonlite)
  cid2het=getURL(paste('https://www.ebi.ac.uk/unichem/rest/src_compound_id/',cid,'/22/3',sep=''))
  cid2het=fromJSON(cid2het)
  return(cid2het$src_compound_id)
}

library(clustermq)
het=Q(getHET_from_CID,cid=cidvec,n_jobs = 50)
names(het)=names(cidvec)
cids$het=sapply(1:nrow(cids),function(i){
  ch=het[[i]]
  if(is.null(ch)){NA}
  else(ch)
})

saveRDS(cids,file='../shared/chembl/CMap.rds')

#################################

drugage=read.csv('./data/raw/DrugAge/DrugAge.csv')
library(clustermq)
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
drlist=unique(as.character(drugage$compound_name))
cids=Q(get_cid,nm=drlist,n_jobs = 50)
names(cids)=drlist

cids=reshape2::melt(cids[!sapply(cids,is.null)])
cidvec=setNames(cids$value,cids$L1)

getChembl_from_CID=function(cid){
  library(RCurl)
  library(jsonlite)
  cid2chembl=getURL(paste('https://www.ebi.ac.uk/unichem/rest/src_compound_id/',cid,'/22/1',sep=''))
  cid2chembl=fromJSON(cid2chembl)
  return(cid2chembl$src_compound_id)
}

library(clustermq)
chembl=Q(getChembl_from_CID,cid=cidvec,n_jobs = 50)
names(chembl)=names(cidvec)
cids$chembl=sapply(1:nrow(cids),function(i){
  ch=chembl[[i]]
  if(is.null(ch)){NA}
  else(ch)
})

getHET_from_CID=function(cid){
  library(RCurl)
  library(jsonlite)
  cid2het=getURL(paste('https://www.ebi.ac.uk/unichem/rest/src_compound_id/',cid,'/22/3',sep=''))
  cid2het=fromJSON(cid2het)
  return(cid2het$src_compound_id)
}

het=Q(getHET_from_CID,cid=cidvec,n_jobs = 50)
names(het)=names(cidvec)
cids$het=sapply(1:nrow(cids),function(i){
  ch=het[[i]]
  if(is.null(ch)){NA}
  else(ch)
})
head(cids)
saveRDS(cids,file='../shared/chembl/DrugAge.rds')
