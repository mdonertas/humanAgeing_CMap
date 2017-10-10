rm(list=ls())
#################
# library(RColorBrewer)
library(data.table)
# library(pheatmap)
# library(biomaRt)
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
k.calc=function(drugmat,i,myprobelist,inst,filename){
  pdf(filename,useDingbats = F,width=10)
  par(mfrow=c(2,3))
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
      hist(x1,br=15,col='gray',main=paste('expectation - no effect'),xlab='Rank')
      hist(x2,br=15,col='gray35',main=paste(nm,'genes in',inst$cmap_name[inst$instance_id==i]),xlab='Rank')
      plot(x1,type='o',pch=19,cex=0.3,main='ks-test',col='gray',ylab='Rank')
      points(x2,type='o',pch=10,cex=0.3,col='gray35')
      arrows(x0=mylim,x1=mylim,y0=x1[mylim],y1=x2[mylim],code = 3,lwd=2,length = 0.1,col='red')
      text(x=mylim+10,y=x2[mylim]-0.01,labels = paste('k=',round(k,2),sep=''),col='red',font=2)
    } else {
      k= -b
      mylim=which.max(x2-x3)
      hist(x3,br=15,col='gray',main=paste('expectation - no effect'),xlab='Rank')
      hist(x2,br=15,col='gray35',main=paste(nm,'genes in',inst$cmap_name[inst$instance_id==i]),xlab='Rank')
      plot(x3,type='o',pch=19,cex=0.3,main='ks-test',col='gray',ylab='Rank')
      points(x2,type='o',pch=10,cex=0.3,col='gray35')
      arrows(x0=mylim,x1=mylim,y0=x1[mylim],y1=x2[mylim],code = 3,lwd=2,length = 0.1,col='red')
      text(x=mylim+10,y=x3[mylim]-0.01,labels = paste('k=',round(k,2),sep=''),col='red',font=2)
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
  dev.off()
  return(res)
}
i=7338
myprobelist=list(up=br_up,down=br_down)
filex='~/GD_ebi/TAC/1st/presentation/figures/extra/kstest_7338.pdf'
k.calc(drugmat = drugmat,i = i,myprobelist = myprobelist,inst = inst,filename=filex)
