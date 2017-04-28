rm(list=ls())
#################
library(scales)
library(RColorBrewer)
library(data.table)
library(pheatmap)
library(biomaRt)
#################
br_up=readRDS('./data/processed/humanBrainMicroarray/all_up_affy.rds')
br_down=readRDS('./data/processed/humanBrainMicroarray/all_down_affy.rds')
gtex_up=as.character(read.table('./data/processed/GTEx/nobrain_up.grp')[,1])
gtex_down=as.character(read.table('./data/processed/GTEx/nobrain_down.grp')[,1])
drugmat=fread('./data/raw/CMap/rankMatrix.txt')
inst=read.csv('./data/raw/CMap/cmap_instances_02.csv',nrows = 6100)
plist=drugmat$V1
drugmat$V1=NULL
colnames(drugmat)=as.character(inst$instance_id)
rownames(drugmat)=plist
br_detailedRes=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_detailed.csv')
gtex_detailedRes=read.csv('./data/processed/GTEX/CMap_results/nobrain_consistent_detailed.csv')
br_CMapRes=read.csv('./data/processed/humanBrainMicroarray/CMap_results/all_consistent_byname.csv')
gtex_CMapRes=read.csv('./data/processed/GTEX/CMap_results/nobrain_consistent_byname.csv')
#################

psx=as.numeric(as.character(br_CMapRes$p))
br_CMapRes$padj=p.adjust(psx,method='fdr')
br_cmap=as.character(br_CMapRes[which(br_CMapRes$padj<0.05),2])
psx=as.numeric(as.character(gtex_CMapRes$p))
gtex_CMapRes$padj=p.adjust(psx,method='fdr')
gtex_cmap=as.character(gtex_CMapRes[which(gtex_CMapRes$padj<0.05),2])


signif=setNames(rep(1,length(inst$instance_id)),inst$cmap_name)
signif[names(signif)%in%br_cmap]=2
signif[names(signif)%in%gtex_cmap]=3
signif[names(signif)%in%intersect(gtex_cmap,br_cmap)]=4

col_annot=data.frame(p=c('-','Brain','GTEx','Both')[signif],
                     BrainScore=setNames(br_detailedRes$score,br_detailedRes$instance_id)[as.character(inst$instance_id)],
                     GTExScore=setNames(gtex_detailedRes$score,gtex_detailedRes$instance_id)[as.character(inst$instance_id)],
                     DrugName=inst$cmap_name,CellType=inst$cell2)
rownames(col_annot)=inst$instance_id

row_annot=data.frame(BrainReg=c('-','Down','Up')[plist%in%c(br_up,br_down) + plist%in%br_up +1],
                     GTExReg=c('-','Down','Up')[plist%in%c(gtex_up,gtex_down) + plist%in%gtex_up +1]
                     )
rownames(row_annot)=plist

mymat=as.matrix(drugmat[plist%in%c(br_up,br_down,gtex_up,gtex_down),c(rownames(col_annot)[(which(col_annot$p!='-'))]),with=F])
rownames(mymat)=plist[plist%in%c(br_up,br_down,gtex_up,gtex_down)]

col_annot2=col_annot[colnames(mymat),]
row_annot2=row_annot[rownames(mymat),]

annocols=list(
  # DrugName=setNames(c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(8,'Set3'))[1:length(unique(col_annot2$DrugName))],unique(col_annot2$DrugName)),
              BrainScore=setNames(colorRampPalette(c('midnightblue','white','firebrick4'))(length(seq(-1,1,0.001))),round(seq(-1,1,0.001),4)),
              GTExScore=setNames(colorRampPalette(c('midnightblue','white','firebrick4'))(length(seq(-1,1,0.001))),round(seq(-1,1,0.001),4)),
              p=setNames(c('gray90',brewer.pal(3,'Set1')),unique(col_annot2$p)),
              CellType=setNames(brewer.pal(5,'Pastel1'),unique(col_annot2$CellType)),
              BrainReg=setNames(c('white','lightblue','pink'),c('-','Down','Up')),
              GTExReg=setNames(c('white','lightblue','pink'),c('-','Down','Up')))

# col_annot2$DrugName=NULL
# col_annot2$GTExScore=NULL
# row_annot2$GTExReg=NULL
perc=10
mymat2=t(apply(mymat,1,scale))

pdf('./results/all_heatmap.pdf',width = 20,height = 30)
ph=pheatmap(mymat2,
            annotation_col = col_annot2,
            annotation_colors = annocols,
            show_rownames = F,
            show_colnames = F,
            annotation_row = row_annot2,
            color=colorRampPalette(brewer.pal(11,'RdYlBu'))(100)[10:90],
            # color=c(muted('pink'),'white',muted('lightblue')),
            # kmeans=20,
            # breaks = c(0,perc,100-perc,100)*nrow(drugmat)/100,
            cutree_cols = 10,
            cutree_rows = 15,
            border_color = 'gray95')
dev.off()