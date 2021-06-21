library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)

setwd("/Volumes/easystore/Human_DCM/Seurat/sctransform")

## GATHERING DATA TOGETHER

HDCM1.data<-Read10X(data.dir="/Volumes/easystore/Human_DCM/HDCM1/outs/filtered_feature_bc_matrix")
HDCM3.data<-Read10X(data.dir="/Volumes/easystore/Human_DCM/HDCM3/outs/filtered_feature_bc_matrix")
HDCM4.data<-Read10X(data.dir="/Volumes/easystore/Human_DCM/HDCM4/outs/filtered_feature_bc_matrix")
HDCM5.data<-Read10X(data.dir="/Volumes/easystore/Human_DCM/HDCM5/outs/filtered_feature_bc_matrix")
HDCM6.data<-Read10X(data.dir="/Volumes/easystore/Human_DCM/HDCM6_1/outs/filtered_feature_bc_matrix")
HDCM7.data<-Read10X(data.dir="/Volumes/easystore/Human_DCM/HDCM7/outs/filtered_feature_bc_matrix")
HDCM8.data<-Read10X(data.dir="/Volumes/easystore/Human_DCM/HDCM8/outs/filtered_feature_bc_matrix")


obj.HDCM1 <- CreateSeuratObject(counts=HDCM1.data, project = "HDCM1")
obj.HDCM3 <- CreateSeuratObject(counts=HDCM3.data, project = "HDCM3")
obj.HDCM4 <- CreateSeuratObject(counts=HDCM4.data, project = "HDCM4")
obj.HDCM5 <- CreateSeuratObject(counts=HDCM5.data, project = "HDCM5")
obj.HDCM6 <- CreateSeuratObject(counts=HDCM6.data, project = "HDCM6")
obj.HDCM7 <- CreateSeuratObject(counts=HDCM7.data, project = "HDCM7")
obj.HDCM8 <- CreateSeuratObject(counts=HDCM8.data, project = "HDCM8")

obj.HDCM1$orig.ident<-'HDCM1'
obj.HDCM3$orig.ident<-'HDCM3'
obj.HDCM4$orig.ident<-'HDCM4'
obj.HDCM5$orig.ident<-'HDCM5'
obj.HDCM6$orig.ident<-'HDCM6'
obj.HDCM7$orig.ident<-'HDCM7'
obj.HDCM8$orig.ident<-'HDCM8'

obj.HDCM1<-RenameCells(obj.HDCM1,add.cell.id = "HDCM1")
obj.HDCM3<-RenameCells(obj.HDCM3,add.cell.id = "HDCM3")
obj.HDCM4<-RenameCells(obj.HDCM4,add.cell.id = "HDCM4")
obj.HDCM5<-RenameCells(obj.HDCM5,add.cell.id = "HDCM5")
obj.HDCM6<-RenameCells(obj.HDCM6,add.cell.id = "HDCM6")
obj.HDCM7<-RenameCells(obj.HDCM7,add.cell.id = "HDCM7")
obj.HDCM8<-RenameCells(obj.HDCM8,add.cell.id = "HDCM8")

obj.HDCM1<-subset(obj.HDCM1, nCount_RNA>2000 & nCount_RNA < 10000)
obj.HDCM3<-subset(obj.HDCM3, nCount_RNA>2000 & nCount_RNA < 10000)
obj.HDCM4<-subset(obj.HDCM4, nCount_RNA>2000 & nCount_RNA < 10000)
obj.HDCM5<-subset(obj.HDCM5, nCount_RNA>2000 & nCount_RNA < 10000)
obj.HDCM6<-subset(obj.HDCM6, nCount_RNA>2000 & nCount_RNA < 10000)
obj.HDCM7<-subset(obj.HDCM7, nCount_RNA>2000 & nCount_RNA < 10000)
obj.HDCM8<-subset(obj.HDCM8, nCount_RNA>2000 & nCount_RNA < 10000)


M1<-merge(obj.HDCM1,obj.HDCM3)
M2<-merge(obj.HDCM5,obj.HDCM4)
M3<-merge(obj.HDCM7,obj.HDCM6)
M4<-merge(M1,M2)
M5<-merge(M3,obj.HDCM8)

HDCM<-merge(M4,M5)

HDCM<-SubsetData(HDCM,cells=cells)#List of cells from final clustered object after each round of doublet removal, repeated until it appears no doublets or contamination remain
rm(list=setdiff(ls(), "HDCM"))

save(HDCM, "HDCM_raw.Robj")

## NORMALIZATION
HDCM <- PercentageFeatureSet(HDCM, pattern = "^MT-", col.name = "percent.mito")

HDCM <- subset(x = HDCM, subset = percent.mito < 10)

HDCM <- SCTransform(HDCM, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)
gc()

## PCA
HDCM <- RunPCA(HDCM, npcs=60, verbose = TRUE)
ElbowPlot(object = HDCM, ndims=60)

## TSNE
HDCM <- RunTSNE(object = HDCM, dims = 1:40)

## UMAP
HDCM <- RunUMAP(object = HDCM, dims = 1:40, verbose=FALSE)

## SAVING
save(HDCM, file="HDCM.Robj")

## CLUSTERING
HDCM <- FindNeighbors(object = HDCM, dims = 1:40)
HDCM <- FindClusters(object = HDCM, reduction.type = "pca",
                     dims = 1:40,
                     resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
                     print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
Idents(HDCM)<-HDCM$SCT_snn_res.0.6
DimPlot(object = HDCM, reduction = 'umap',pt.size=0.3)
ggsave("UMAP_0.6.png", width=6, height=6)

## SAVING
save(HDCM, file="HDCM.Robj")

## FINDING ANS SAVING MARKERS
HDCM.markers <- FindAllMarkers(object = HDCM,
                               only.pos = TRUE,
                               min.pct = 0.10,
                               thresh.use = 0.10)
write.table(HDCM.markers, "markers.tsv", sep="\t", quote=F, row.names=F)

## SAVING
save(HDCM, file="HDCM.Robj")

#CLUSTER TABLE BY GROUP
clusters.by.group<-table(HDCM@active.ident, HDCM@meta.data$orig.ident)
write.table(clusters.by.group, "clusters_by_group_0.6.tsv", sep="\t", quote=F)

#HEATMAP
markers<-FindAllMarkers(object=HDCM, only.pos=TRUE, min.pct =0.25)
write.table(markers, "markers_0.6.tsv", sep="\t", quote=F, row.names=F)

top5 <- HDCM.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

my_levels<-c(0,1,10,11,12,13,14,15,16,17,18,19,2,20,3,4,5,6,7,8,9)
levels(HDCM)<-my_levels
levels(top5)<-my_levels

Heatmap<-DoHeatmap(object = HDCM, features = top5$gene, size=2, angle=90) + NoLegend()
Heatmap + theme(axis.text.y = element_text(size=5), plot.margin=unit(c(1,1,1,1),"cm"))
ggsave("Top5_Heatmap_0.6.png", width=16, height=12)

save(HDCM, file="HDCM.Robj")
