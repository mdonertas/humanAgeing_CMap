library(pheatmap)
library(RColorBrewer)
library(data.table)
br_up=readRDS('./data/processed/humanBrainMicroarray/all_up_affy.rds')
br_down=readRDS('./data/processed/humanBrainMicroarray/all_down_affy.rds')
drugmat=fread('./data/raw/CMap/amplitudeMatrix.txt')
inst=read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
plist=drugmat$V1
drugmat$V1=NULL


colnames(drugmat)=strsplit(readLines('./data/raw/CMap/amplitudeMatrix.txt',n = 1),'\t')[[1]][-1]
rownames(drugmat)=plist
br_detailedRes=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_detailed.csv')
br_CMapRes=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')

br_CMapRes$p.adj=p.adjust(as.numeric(as.character(br_CMapRes$p)),method='fdr')
inst2=inst[inst$cmap_name%in%as.character(br_CMapRes$cmap.name[which(br_CMapRes$p.adj<0.05)]),]


drugmat_all_sc=apply(drugmat[,colnames(drugmat)%in%as.character(inst2$instance_id),with=F],2,scale)


drugmat_up=drugmat_all_sc[rownames(drugmat)%in%br_up,]
annotcols=data.frame(score=br_detailedRes$score[br_detailedRes$instance_id%in%colnames(drugmat_up)],
name=br_detailedRes$cmap.name[br_detailedRes$instance_id%in%colnames(drugmat_up)])
rownames(annotcols)=br_detailedRes$instance_id[br_detailedRes$instance_id%in%colnames(drugmat_up)]
annocols=list(
  name=setNames(c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(8,'Set3'),brewer.pal(8,'Dark2'),brewer.pal(9,'Pastel1'))[1:length(unique(annotcols$name))],unique(annotcols$name)),
  score=setNames(colorRampPalette(c('blue','white','red'))(length(seq(-1,1,0.001))),round(seq(-1,1,0.001),4)))
drugmat2=drugmat_up[,names(sort(setNames(annotcols$score,rownames(annotcols)),dec=T))]
pheatmap(drugmat2,
         annotation_col = annotcols,annotation_colors = annocols,
         color=rev(colorRampPalette(brewer.pal(10,'RdYlBu'))(25)),
         breaks=seq(-8.26,8.26,length.out = 26),
         # cluster_cols = F,
         # main = 'Drug Perturbed Expression Profil for Genes Up-regulated in Ageing',
         show_colnames = F,
         filename = '~/Desktop/drugPerturbation_ups.pdf')

drugmat_down=drugmat_all_sc[rownames(drugmat)%in%br_down,]
dim(drugmat_down)
annotcols=data.frame(score=br_detailedRes$score[br_detailedRes$instance_id%in%colnames(drugmat_down)],
                     name=br_detailedRes$cmap.name[br_detailedRes$instance_id%in%colnames(drugmat_down)])
rownames(annotcols)=br_detailedRes$instance_id[br_detailedRes$instance_id%in%colnames(drugmat_down)]
drugmat2=drugmat_down[,names(sort(setNames(annotcols$score,rownames(annotcols)),dec=T))]
pheatmap(drugmat2,
         annotation_col = annotcols,annotation_colors = annocols,
         # cluster_cols = F,
         color=rev(colorRampPalette(brewer.pal(10,'RdYlBu'))(25)),
         breaks=seq(-7,7,length.out = 26),
         # main = 'Drug Perturbed Expression Profil for Genes Down-regulated in Ageing',
         show_colnames = F,
         filename = '~/Desktop/drugPerturbation_downs.pdf')
