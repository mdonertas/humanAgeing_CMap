ages=readRDS('./data/processed/ages/Maycox2009_APFC_age.rds')
permage=lapply(1:1000,function(i)setNames(sample(ages),names(ages)))
saveRDS(permage,file='./data/processed/ages/permutations/Maycox2009_APFC_age.rds')

rm(list=ls())

ages=readRDS('./data/processed/ages/Barnes2011_STC_age.rds')
permage=lapply(1:1000,function(i)setNames(sample(ages),names(ages)))
saveRDS(permage,file='./data/processed/ages/permutations/Barnes2011_STC_age.rds')

rm(list=ls())

ages=readRDS('./data/processed/ages/Lu2004_FC_age.rds')
permage=lapply(1:1000,function(i)setNames(sample(ages),names(ages)))
saveRDS(permage,file='./data/processed/ages/permutations/Lu2004_FC_age.rds')

rm(list=ls())

##########Berchtold

library(RCurl)
library(GEOquery)
g=getGEO('GSE11882')
pd=pData(g[[1]])
rm(g)
xx=data.frame(individual=sapply(strsplit(as.character(pd$title),'_'),function(x)x[4]),age=as.numeric(gsub('[age (yrs):]','',pd$characteristics_ch1.3)),id=rownames(pd))
ages=unique(data.frame(ind=xx$individual,age=xx$age))
gsm=setNames(as.character(xx$individual),as.character(xx$id))
permage=lapply(1:1000,function(i){
  ages=setNames(sample(ages$age),ages$ind)
  aa=setNames(ages[gsm],names(gsm))
  aa
})
saveRDS(permage,file='./data/processed/ages/permutations/Berchtold2008_EC_age.rds')
saveRDS(permage,file='./data/processed/ages/permutations/Berchtold2008_HIP_age.rds')
saveRDS(permage,file='./data/processed/ages/permutations/Berchtold2008_PCG_age.rds')
saveRDS(permage,file='./data/processed/ages/permutations/Berchtold2008_SFG_age.rds')

#########others
rm(list=ls())
permage=readRDS('~/Desktop/permage.rds')
nms=names(permage[[1]])
ages=lapply(nms,function(nm){
  lapply(permage,function(x){
    x[[nm]]
  })
})
names(ages)=nms
sapply(nms,function(nm){
  age=ages[[nm]]
  saveRDS(age,file=paste('./data/processed/ages/permutations/',nm,'.rds',sep=''))
})

############
rm(list=ls())

fx=list.files('./data/processed/ages/permutations/')
length(fx)
