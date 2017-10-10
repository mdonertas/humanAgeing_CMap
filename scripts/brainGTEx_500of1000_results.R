xx=readRDS('./data/processed/brainGTEx/500of1000.rds')
xx2=sapply(1:length(xx[[1]]),function(i)sapply(xx,function(x)x[i]))
colnames(xx2)=names(xx[[1]])
gtex=read.csv('./data/processed/brainGTEx/CMap_results/GTEx_brain.csv')
pdf('./results/gtex500of1000.pdf',width = 10,height = 8)
par(mfrow=c(4,5))
resx=sapply(colnames(xx2),function(nm){
  distx=xx2[,nm]
  realx=gtex[gtex$cmap.name==nm,'mean']
  p=mean(realx>=distx)
  p=2*ifelse(p>0.5,1-p,p)
  hist(distx,br=25,main=paste(nm,'\np=',round(p,3),'\nmedian=',round(median(distx),3),sep=''),xlim=c(min(c(realx,distx)),max(c(realx,distx))),col='gray80')
  abline(v=realx,col='red',lwd=2)
  p
})
dev.off()
which(p.adjust(resx,method='fdr')<0.05)

gtex[gtex$cmap.name%in%colnames(xx2),]
