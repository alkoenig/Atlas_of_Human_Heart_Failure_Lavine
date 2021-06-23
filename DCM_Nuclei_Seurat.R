library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tibble)
library(sctransform)
library(purrr)

setwd("/scratch/alkoenig/Nuclei/Subsets")

## GATHERING DATA TOGETHER
dataFolders<-as.matrix(read.csv("/scratch/alkoenig/Nuclei/Subsets/folders_chpc.txt", header=FALSE))
dataFolders<-as.list(dataFolders[,1])
sampleNames<-as.matrix(read.csv("/scratch/alkoenig/Nuclei/Subsets/names.txt", header=FALSE))
sampleNames<-as.list(sampleNames[,1])

raw.data<-as.list(NULL)
for(i in 1:length(dataFolders)){
  raw.data[i]<-Read10X(data.dir=dataFolders[[i]]) 
}


objs<-as.list(NULL)
for(i in 1:length(raw.data)){
  objs[i] <- CreateSeuratObject(counts=raw.data[[i]], project = sampleNames[[i]])
  objs[[i]]$orig.ident<-sampleNames[[i]]                              
  objs[[i]]<-RenameCells(objs[[i]],add.cell.id = sampleNames[[i]])
}

rm(raw.data)

#Save Raw List of Objects
save(objs,file="obj_list.Robj")

#Merging
nuclei<-reduce(objs,merge)
rm(objs)

#Save Raw object
save(nuclei,file="raw_nuclei.Robj")

#QC
nuclei <- PercentageFeatureSet(nuclei, pattern = "^MT-", col.name = "percent.mito")
nuclei<-subset(nuclei, percent.mito<5)
nuclei<-subset(whole, nCount_RNA > 1000 & nCount_RNA < 10000)

#Normalize
nuclei <- SCTransform(nuclei, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory=TRUE)

## PCA
nuclei <- RunPCA(nuclei, npcs=80, verbose = TRUE)
ElbowPlot(object = nuclei, ndims=80) 

## TSNE
HDCM <- RunTSNE(object = HDCM, dims = 1:40)

## UMAP
HDCM <- RunUMAP(object = HDCM, dims = 1:40, verbose=FALSE)

## SAVING
save(nuclei, file="nuclei.Robj")

## CLUSTERING
nuclei <- FindNeighbors(object = nuclei, reduction='pca', dims = 1:80)
nuclei <- FindClusters(object = nuclei,
                       resolution = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
Idents(nuclei)<-nuclei$SCT_snn_res.0.5
DimPlot(object = nuclei, reduction = 'umap',pt.size=0.5)
ggsave("UMAP_0.5.png", width=12, height=6)

## FINDING ANS SAVING MARKERS
markers <- FindAllMarkers(object = nuclei,
                               only.pos = TRUE,
                               min.pct = 0.10,
                               thresh.use = 0.10)
write.table(markers, "nuclei.markers_0.5.tsv", sep="\t", quote=F, row.names=F)

## SAVING
save(nuclei, file="nuclei.Robj")

## Add Condition meta
nuclei$Condition<- ifelse(nuclei$orig.ident %in% c('H_ZC-LVAD',
                                                   'TWCM-11-3',
                                                   'TWCM-11-93',
                                                   'TWCM-13-17',
                                                   'TWCM-13-47',
                                                   'TWCM-13-84',
                                                   'TWCM-13-102',
                                                   'TWCM-13-181',
                                                   'TWCM-13-208',
                                                   'TWCM-13-285',
                                                   'TWCM-LVAD2',
                                                   'TWCM-LVAD3'), "DCM", "Donor")

## Add Age Group Meta
nuclei$Age_Group_Tertile<- ifelse(nuclei$orig.ident %in% c('TWCM-13-181',
                                                               'TWCM-11-93',
                                                               'TWCM-13-198',
                                                               'TWCM-13-235',
                                                               'TWCM-13-96',
                                                               'TWCM-13-192',
                                                               'TWCM-11-78',
                                                               'TWCM-13-285',
                                                               'TWCM-13-152',
                                                               'TWCM-11-256',
                                                               'TWCM-13-47',
                                                               'TWCM-13-102',
                                                               'TWCM-11-82'), "Young",ifelse(nuclei$orig.ident %in% c('TWCM-11-103',
                                                                                                                   'TWCM-13-80',
                                                                                                                   'TWCM-11-42',
                                                                                                                   'TWCM-13-104',
                                                                                                                   'TWCM-13-132',
                                                                                                                   'TWCM-13-208',
                                                                                                                   'TWCM-10-5',
                                                                                                                   'TWCM-13-84',
                                                                                                                   'H_ZC-LVAD',
                                                                                                                   'TWCM-11-41',
                                                                                                                   'TWCM-13-17',
                                                                                                                   'TWCM-13-101',
                                                                                                                   'TWCM-11-3',
                                                                                                                   'TWCM-13-1',
                                                                                                                   'TWCM-LVAD3'),"Middle","Old"))
## Add Sex Meta
nuclei$Sex<- ifelse(nuclei$orig.ident %in% c('H_ZC-11-292',
                                                 'H_ZC-LVAD',
                                                 'TWCM-11-41',
                                                 'TWCM-11-74',
                                                 'TWCM-11-78',
                                                 'TWCM-11-82',
                                                 'TWCM-11-93',
                                                 'TWCM-11-103',
                                                 'TWCM-11-192',
                                                 'TWCM-13-47',
                                                 'TWCM-13-80',
                                                 'TWCM-13-84',
                                                 'TWCM-13-102',
                                                 'TWCM-13-152',
                                                 'TWCM-13-168',
                                                 'TWCM-13-192',
                                                 'TWCM-13-198',
                                                 'TWCM-13-208',
                                                 'TWCM-LVAD3'), "Male", "Female")

#Add Names Meta
current.cluster.ids<-levels(nuclei)                  

new.cluster.ids<-c('Fibroblasts',
                   'Cardiomyocytes',
                   'Fibroblasts',
                   'Cardiomyocytes',
                   'Cardiomyocytes',
                   'Lymphatic',
                   'Cardiomyocytes',
                   'Neurons',
                   'Endothelium',
                   'Epicardium',
                   'Mast_Cells',
                   'Adipocytes',
                   'Macrophages',
                   'Macrophages',
                   'Monocytes',
                   'Macrophages',
                   'Cardiomyocytes',
                   'Fibroblasts',
                   'Pericytes',
                   'B_Cells',
                   'Cardiomyocytes',
                   'Endothelium',
                   'Pericytes',
                   'Endocardium',
                   'Endothelium',
                   'Smooth_Muscle',
                   'T/NK_Cells',
                   'Fibroblasts')

Idents(nuclei) <- plyr::mapvalues(x = nuclei$SCT_snn_res.0.5, from = current.cluster.ids, to = new.cluster.ids)

nuclei$Names<-Idents(nuclei)
