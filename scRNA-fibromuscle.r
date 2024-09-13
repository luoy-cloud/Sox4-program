rm(list=ls())
# Load packages ----
library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(paletteer)
library(patchwork)
setwd("E:/GSE233955_RAW")
# Load data ----
E10.data <- ReadMtx(mtx = 'GSM6604425_E10_Ebf2GFP_matrix.mtx.gz',
                    cells = 'GSM6604425_E10_Ebf2GFP_barcodes.tsv.gz',
                    features = 'GSM6604425_E10_Ebf2GFP_features.tsv.gz')
E11.data <- ReadMtx(mtx = 'GSM6604426_E11_Pdgfra_matrix.mtx.gz',
                    cells = 'GSM6604426_E11_Pdgfra_barcodes.tsv.gz',
                    features = 'GSM6604426_E11_Pdgfra_features.tsv.gz')
E12.data <- ReadMtx(mtx = 'GSM6604427_E12_Pdgfra_matrix.mtx.gz',
                    cells = 'GSM6604427_E12_Pdgfra_barcodes.tsv.gz',
                    features = 'GSM6604427_E12_Pdgfra_features.tsv.gz')
E13.data <- ReadMtx(mtx = 'GSM6604428_E13_Pdgfra_matrix.mtx.gz',
                    cells = 'GSM6604428_E13_Pdgfra_barcodes.tsv.gz',
                    features = 'GSM6604428_E13_Pdgfra_features.tsv.gz')
# Create Seurat objects and consolidate into list ----
E10 <- CreateSeuratObject(E10.data, min.cells = 3, project = 'E10')
E11 <- CreateSeuratObject(E11.data, min.cells = 3, project = 'E11')
E12 <- CreateSeuratObject(E12.data, min.cells = 3, project = 'E12')
E13 <- CreateSeuratObject(E13.data, min.cells = 3, project = 'E13')
list3 <- c('E10' = E10,
           'E11' = E11,
           'E12' = E12,
           'E13' = E13
)

# Determine mitochondrial, ribosomal, and replication-dependent histone scores
list3 <- lapply(X = list3, FUN = function(x) {
  x[['percent.mt']] <- PercentageFeatureSet(x, pattern = "^mt-")
  x[['percent.ribo']] <- PercentageFeatureSet(x, pattern = "^Rp[sl]")
  x[['percent.hist']] <- PercentageFeatureSet(x, pattern = "^Hist")
  return(x)
})
# Choose cutoffs 
list3.minFeatures <- c('E10' = 2500,
                       'E11' = 2000,
                       'E12' = 2000,
                       'E13' = 1500
)
list3.maxMito <- c('E10' = 7.5,
                   'E11' = 7,
                   'E12' = 7.5,
                   'E13' = 7.5
)

# Import doublet scores from Scrublet and remove doublets
list3 <- lapply(X = list3, FUN = function(x) {
  doublet_scores <- read.table(paste0('E:/ws/code/scrublet/',x@project.name,'/doublet_scores.txt'),
                               quote="\"", comment.char="", na.strings="NM")
  x <- AddMetaData(x, doublet_scores$V1, col.name = "DS")
  x <- subset(x, subset = DS < 0.35)
  return(x)
})

# Subset objects based on cutoffs  (nFeature_RNA > minGene & nCount_RNA < maxUMI &#nFeature_RNA < maxGene & percent.mt < pctMT )
list3 <- lapply(X = list3, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > list3.minFeatures[x@project.name] &
                nFeature_RNA < quantile(x$nFeature_RNA, 0.99) & 
                percent.mt < list3.maxMito[x@project.name])
  return(x)
})
# Import gene lists for cell cycle scoring
m.s.genes <- readRDS("E:/cell-cycle_mouse/mouse.S.genes.RDS")
m.g2m.genes <- readRDS("E:/cell-cycle_mouse/mouse.G2M.genes.RDS")
# Process
list3 <- lapply(X = list3, FUN = function(x) {
  x <- NormalizeData(x, verbose = F)
  x <- FindVariableFeatures(x, verbose = F)
  x <- CellCycleScoring(x, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = F)
  x <- ScaleData(x, vars.to.regress = c('S.Score','G2M.Score','percent.mt','percent.hist'))
  x <- RunPCA(x, verbose = F)
  x <- RunUMAP(x, dims = 1:30, verbose = F)
  x <- FindNeighbors(x, dims = 1:30, verbose = F)
  x <- FindClusters(x, resolution = 0.3, verbose = F)
  return(x)
})

# Select muscle and fibroblast subsets for E10-E13 (Early) integration
list3 <- lapply(list3, FindClusters, resolution = 0.3)
lapply(list3[1:4], function(x) {
  p1 <- FeaturePlot(x, features = c('Pdgfra')) + theme(plot.title = element_text(face = 'italic')) & NoAxes() & NoLegend() 
  p2 <- FeaturePlot(x, features = c('Tnnt1')) & NoAxes() & NoLegend() + theme(plot.title = element_text(face = 'italic')) & NoAxes() & NoLegend() 
  p3 <- DimPlot(x, label = T) & NoAxes() & NoLegend()
  plot((p3 | p1 | p2) & plot_annotation(title = x@project.name, theme = theme(plot.title = element_text(hjust = 0.5, size = 24))))
  return(NULL)
})


list3.early.fibro.muscle.clusters <- list('E10' = c(0:4),
                                          'E11' = c(0,1,3:5,7),
                                          'E12' = c(0,1:3,5,12),
                                          'E13' = c(0:6))

## Subset selected clusters from datasets ####
list3.early <- list3[1:4]
list3.early <- lapply(X = list3.early, FUN = function(x) {
  x <- subset(x, idents = list3.early.fibro.muscle.clusters[[x@project.name]])
  return(x)
})


## Integrate ####
Early.features <- SelectIntegrationFeatures(list3.early, nfeatures = 3000)
list3.early <- lapply(X = list3.early, FUN = function(x) {
  x <- ScaleData(x, features = Early.features)
  x <- RunPCA(x, features = Early.features, verbose = FALSE)
  return(x)
})
#FindIntegrationAnchors
Early.anchors <- FindIntegrationAnchors(object.list = list3.early, 
                                        normalization.method = "LogNormalize",
                                        anchor.features = Early.features,
                                        reduction = 'rpca')

Early <- IntegrateData(Early.anchors, normalization.method = "LogNormalize")
DefaultAssay(Early) <- "integrated"
Early <- ScaleData(Early, vars.to.regress = c('S.Score','G2M.Score','percent.mt'))
Early <- RunPCA(Early, verbose = FALSE, npcs = 100)
Early <- FindNeighbors(Early, reduction = "pca", dims = 1:50)
Early <- FindClusters(Early, resolution = 0.6)
DefaultAssay(Early) <- "RNA"
Early <- RunUMAP(Early, reduction = "pca", dims = 1:50)
########umap
DimPlot(Early, 
        label = T,
        label.size = 5,
        reduction = "umap", 
        cols= paletteer_d("ggsci::nrc_npg"), 
        pt.size = 0.5) 