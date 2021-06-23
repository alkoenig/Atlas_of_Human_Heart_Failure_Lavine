library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(SeuratDisk)

## Gathering And Preparing Data
load("./DCM_Cells_Filtered.Robj") #Filtered cell object used as query (This object is the raw untransformed/unclustered object filtered to include only cells in the final clustered cell object)
load("./DCM_Nuclei.Robj") #Final clustered nuclei object used as reference

HDCM[["percent.mt"]] <- PercentageFeatureSet(HDCM, pattern = "^MT-")
HDCM <- SCTransform(HDCM, verbose = FALSE, conserve.memory=TRUE, vars.to.regress = "percent.mt")

##Find Anchors
anchors <- FindTransferAnchors(
  reference = nuclei,
  query = HDCM,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:80, recompute.residuals = FALSE
)
save(anchors,file="reference_anchors.Robj")

##Mapping
HDCM <- MapQuery(
  anchorset = anchors,
  query = HDCM,
  reference = nuclei,
  refdata = list(
    celltype.l1 = "seurat_clusters",
    celltype.l2 = "names"
    ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

save(HDCM,file="HDCM_RefMapped.Robj")

##Pre-Merge Plots
png("RefUMAP_clusters.png",
       width = 10*100, height = 8*100,
        res=100)
DimPlot(object = HDCM, reduction = 'ref.umap',pt.size=0.1, group.by="predicted.celltype.l1")
dev.off()

png("RefUMAP_names.png",
          width = 10*100, height = 8*100,
          res=100)
DimPlot(object = HDCM, reduction = 'ref.umap',pt.size=0.1, group.by="predicted.celltype.l2")
dev.off()

png("Myeloid_PredictionScore.png",
           width = 20*100, height = 8*100,
           res=100)
FeaturePlot(object = HDCM, features = c("Macrophages","Monocytes"),pt.size=0.1,reduction = "ref.umap",
                           cols =  c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
dev.off()

## Merging
nuclei$id <- 'reference'
HDCM$id <- 'query'
RefMerge <- merge(nuclei, HDCM)
RefMerge[["pca"]] <- merge(nuclei[["pca"]], HDCM[["ref.pca"]])
RefMerge <- RunUMAP(RefMerge, reduction = 'pca', dims = 1:80)

png("Merged_UMAP_IDs.png",
    width = 10*100, height = 8*100,
    res=100)
DimPlot(RefMerge, group.by = 'id',raster=FALSE)
dev.off()

save(RefMerge,file="RefMerged.Robj")

## Merged Plots
png("Merged_UMAP_IDs_Split.png",
    width = 20*100, height = 8*100,
    res=100)
DimPlot(RefMerge, group.by = 'id', split.by = "id",raster = FALSE)
dev.off()

png("Merge_RefUMAP_0.6_Nuclei.png",
    width = 10*100, height = 8*100,
    res=100)
DimPlot(object = RefMerge, reduction = 'umap',pt.size=0.1, group.by="SCT_snn_res.0.6",cells = WhichCells(RefMerge,idents="reference"))
dev.off()

png("Merge_RefUMAP_0.6_Cells.png",
    width = 10*100, height = 8*100,
    res=100)
DimPlot(object = RefMerge, reduction = 'umap',pt.size=0.1, group.by="predicted.celltype.l1",cells = WhichCells(RefMerge,idents="query"))
dev.off()

##Re-clustering
RefMerge <- FindNeighbors(object = RefMerge, dims = 1:80)
RefMerge <- FindClusters(object = RefMerge, reduction.type = "pca",
                                       dims = 1:80,
                                       resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
                                       print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
Idents(RefMerge)<-RefMerge$SCT_snn_res.0.6

##Re-clustered plots
png("UMAP_0.6.png",
    width = 10*100, height = 8*100,
    res=100)
DimPlot(object = RefMerge, reduction = 'umap',pt.size=0.1)
dev.off()

png("UMAP_0.6_SplitTech.png",
    width = 20*100, height = 8*100,
    res=100)
DimPlot(object = RefMerge, reduction = 'umap',pt.size=0.1,split.by="id")
dev.off()

save(RefMerge,file="RefMerged_Reclustered.Robj")


## Add Condition Meta
RefMerge$Condition<-ifelse(RefMerge@meta.data$orig.ident %in% c("HDCM1","HDCM3","HDCM4","HDCM6","HDCM8","H_ZC-LVAD",
                                                                "TWCM-11-93",
                                                                "TWCM-13-181",
                                                                "TWCM-13-201",
                                                                "TWCM-13-208",
                                                                "TWCM-LVAD2",
                                                                "TWCM-LVAD3",
                                                                "TWCM-11-3",
                                                                "TWCM-13-285",
                                                                "TWCM-13-280",
                                                                "TWCM-13-47",
                                                                "TWCM-13-102",
                                                                "TWCM-13-84",
                                                                "TWCM-13-17"), "DCM","Donor")

##Adding Age Tertile Meta
RefMerge$Age_Group_Tertile<- ifelse(RefMerge$orig.ident %in% c('TWCM-13-181',
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
                                                               'TWCM-11-82',
                                                               "HDCM3",
                                                               "HDCM8"), "Young",ifelse(RefMerge$orig.ident %in% c('TWCM-11-103',
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
                                                                                                                   'TWCM-LVAD3',
                                                                                                                   'HDCM1',
                                                                                                                   'HDCM6'),"Middle","Old"))

##Adding Sex Meta
RefMerge$Sex<- ifelse(RefMerge$orig.ident %in% c('H_ZC-11-292',
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
                                                 'TWCM-LVAD3',
                                                 'HDCM1',
                                                 'HDCM6',
                                                 'HDCM7',
                                                 'HDCM4'), "Male", "Female")

##Adding Additional Meta
RefMerge$tech<-ifelse(RefMerge$id =="query","SC","SN")

RefMerge$CondTech<-paste(RefMerge$condition,RefMerge$tech, sep="_")

##Find Marker Genes
markers <- FindAllMarkers(object = RefMerge,
                               only.pos = TRUE,
                               min.pct = 0.10)
write.table(markers, "markers_0.6.tsv", sep="\t", quote=F, row.names=F)

## Adding Cell Names Meta
current.cluster.ids<-levels(RefMerge)

new.cluster.ids<-c('Cardiomyocytes',
                   'Fibroblasts',
                   'NK/T-Cells',
                   'Macrophages',
                   'Smooth_Muscle',
                   'Fibroblasts',
                   'Fibroblasts',
                   'Fibroblasts',
                   'Fibroblasts',
                   'Cardiomyocytes',
                   'Lymphatic',
                   'Cardiomyocytes',
                   'Macrophages',
                   'Neurons',
                   'Fibroblasts',
                   'Cardiomyocytes',
                   'Endothelium',
                   'Epicardium',
                   'Mast',
                   'Adipocytes',
                   "Macrophages",
                   "Macrophages",
                   "Fibroblasts",
                   "Endothelium",
                   "B-Cells",
                   "Cardiomyocytes",
                   "Pericytes",
                   "Pericytes",
                   "Fibroblasts",
                   "Endocardium",
                   "Endothelium",
                   "Endothelium",
                   "Pericytes")

Idents(RefMerge) <- plyr::mapvalues(x = RefMerge$SCT_snn_res.0.6, from = current.cluster.ids, to = new.cluster.ids)

RefMerge$Names<-Idents(RefMerge)

Save(RefMerge,file="RefMerged_Reclustered_AddlMeta.Robj")