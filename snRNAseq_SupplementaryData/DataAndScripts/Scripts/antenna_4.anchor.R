suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
})

theme_set(theme_cowplot())
projectFolder='/Users/tzuchiao/ProjectAnalysis/Leslie_Vosshall'
setwd(projectFolder)



#####################################

### 4. anchor

for (normalizationMethod in c('LogNormalize')) {
  print(c('normalizationMethod: ', normalizationMethod))
  
  anchorFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/4.anchor')
  dir.create(anchorFolderPath, showWarnings = FALSE)
  picFolderPath=paste0(anchorFolderPath, '/pic')
  dir.create(picFolderPath, showWarnings = FALSE)
  picFolderPath=paste0(anchorFolderPath, '/pic')
  dir.create(picFolderPath, showWarnings = FALSE)
  
  # print(anchorFolderPath)
  # print(picFolderPath)
  
  
  seurat.merge.filter <- readRDS(file = file.path('Analysis/integratedData_decontX', 
                                          'antenna_3.merge.filter.rds') )
  print(anchorFolderPath)
  print(picFolderPath)
  
  
  seurat.list <- SplitObject(seurat.merge.filter, split.by = "dataset")
  
  ### 
  if (normalizationMethod == 'LogNormalize'){
    
    # normalize and identify variable features for each dataset independently
    seurat.list <- lapply(X = seurat.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    
    # select features that are repeatedly variable across datasets for integration
    # run PCA on each dataset using these features
    features <- SelectIntegrationFeatures(object.list = seurat.list)
    anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)
    # this command creates an 'integrated' data assay
    seurat.combined <- IntegrateData(anchorset = anchors)
    
  }
  if (normalizationMethod == 'SCTransform') {
    
    seurat.list <- lapply(X = seurat.list, FUN = SCTransform)
    features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
    seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT",
                                             anchor.features = features)
    seurat.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    
  }
  
  
  ### Perform an integrated analysis
  DefaultAssay(seurat.combined) <- "integrated"
  
  if (normalizationMethod == 'LogNormalize'){
    seurat.combined <- ScaleData(seurat.combined, vars.to.regress = c('nCount_RNA'))
  }

  # PCA, UMAP, tsne, clustering
  seurat.combined    <- RunPCA(seurat.combined, npcs = 50, verbose = F)
  seurat.combined    <- RunUMAP(seurat.combined, dims = 1:50, verbose = F)
  seurat.combined    <- RunTSNE(seurat.combined, dims = 1:50, verbose = F)
  seurat.combined    <- FindNeighbors(seurat.combined, dims = 1:50, verbose = F)
  seurat.combined    <- FindClusters(seurat.combined, verbose = T)
  
  
  pdf(file = file.path(picFolderPath, 
                       paste0('UMAP_tSNE_', 'merge.pdf')),
      width=15, height=6)
  print(DimPlot(seurat.combined, reduction = "umap", label = TRUE) + DimPlot(seurat.combined, reduction = "tsne", label = TRUE))
  dev.off()
  
  pdf(file = file.path(picFolderPath, 
                       paste0('UMAP_tSNE_', 'merge_dataset.pdf')),
      width=13, height=6)
  print(DimPlot(seurat.combined, reduction = "umap", group.by = 'dataset') + DimPlot(seurat.combined, reduction = "tsne", group.by = 'dataset'))
  dev.off()
  
  pdf(file = file.path(picFolderPath, 
                       paste0('UMAP_tSNE_', 'merge_dataset_split.pdf')),
      width=20, height=15)
  print(DimPlot(seurat.combined, reduction = "umap", split.by = 'dataset') / 
          DimPlot(seurat.combined, reduction = "tsne", split.by = 'dataset'))
  dev.off()
  
  
  
  saveRDS(seurat.combined, file = 
            file.path('Analysis/integratedData_decontX', 'antenna_4.seurat.combined.rds') )
}



