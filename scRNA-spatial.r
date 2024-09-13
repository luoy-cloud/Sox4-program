# Load packages ----
library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(paletteer)
library(patchwork)
setwd("E:/SpaDecon/")
Early <- readRDS('Early.rds')
# Use spatial locations to narrow E10-E13 dataset (Early.sub1)
DefaultAssay(Early) <- 'integrated'
Early <- FindClusters(Early, resolution = 0.6)
DefaultAssay(Early) <- 'RNA'
p1<-DimPlot(Early, label = T, pt.size = 0.3)  & NoLegend()#& NoAxes()
Early.sub1 <- subset(Early, idents = c(1,3,4,7,10,14,17))
p2<-DimPlot(Early, cells.highlight = WhichCells(Early, idents = c(1,3,4,7,10,14,17)), 
            sizes.highlight = NULL) #& NoLegend() #& NoAxes()
DefaultAssay(Early.sub1) <- "integrated"
Early.sub1 <- ScaleData(Early.sub1, vars.to.regress = c('S.Score','G2M.Score',
                                                        'percent.mt','percent.hist'))
Early.sub1 <- RunPCA(Early.sub1, verbose = FALSE, npcs = 50)
Early.sub1 <- RunUMAP(Early.sub1, reduction = "pca", dims = 1:30)
Early.sub1 <- FindNeighbors(Early.sub1, reduction = "pca", dims = 1:30)
Early.sub1 <- FindClusters(Early.sub1, resolution = 0.3)
DefaultAssay(Early.sub1) <- "RNA"
p3<-DimPlot(Early.sub1, label = T, pt.size = 0.4) & #ggsci::nrc_npg#cols= paletteer_d("ggsci::nrc_npg")
     NoLegend()#& NoAxes()
	 
#findmarkers
data_early.markers <- FindAllMarkers(Early.sub1 , only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
#cell annotation
cluster2celltype <- c("0"="Ebf2+/Cdh4+",
                      "1"="Sclerotome", 
                      "2"="SKM", 
                      "3"= "Dermis", 
                      "4"= "Sclerotome", 
                      "5"= "DM",
                      "6"= "SKM")
Early.sub1[['cell_type']] = unname(cluster2celltype[Early.sub1@meta.data$seurat_clusters])
#umap
DimPlot(Early.sub1, 
            label = T,
            label.size = 5,
            reduction = "umap", 
            group.by = "cell_type", 
            cols= paletteer_d("ggsci::category10_d3"),  
            pt.size = 0.5)
