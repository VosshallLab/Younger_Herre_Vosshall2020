suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(DoubletFinder)
})

theme_set(theme_cowplot())
projectFolder='/Users/tzuchiao/ProjectAnalysis/Leslie_Vosshall'
setwd(projectFolder)


#####################################

### 3. Integrate data 

### Cell filtering

dir.create(file.path('Analysis', 'integratedData_decontX'), showWarnings = FALSE)
for (normalizationMethod in c('LogNormalize')) {
  print(c('normalizationMethod: ', normalizationMethod))
  
  integrateFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/3.integrate')
  dir.create(integrateFolderPath, showWarnings = FALSE)
  picFolderPath=paste0(integrateFolderPath, '/pic')
  dir.create(picFolderPath, showWarnings = FALSE)
  picFolderPath=paste0(integrateFolderPath, '/pic')
  dir.create(picFolderPath, showWarnings = FALSE)
  
  # print(integrateFolderPath)
  # print(picFolderPath)
  
  
  seurat.0819 <- readRDS(file = file.path('Analysis/integratedData_decontX', 
                                         'antenna_2.DoubletFinder.08192021.rds') )
  seurat.0819$dataset <- '0819'
  seurat.1006 <- readRDS(file = file.path('Analysis/integratedData_decontX', 
                                          'antenna_2.DoubletFinder.10062021.rds') )
  seurat.1006$dataset <- '1006'
  
  seurat.merge <- merge(seurat.0819 , y =seurat.1006,
                        add.cell.ids = c("0819", "1006"), project = "antenna")
  seurat.merge.filter <- subset(seurat.merge, subset = nFeature_RNA < 4000 & nFeature_RNA>= 400 & percent.mt < 5)
  
  # print(integrateFolderPath)
  # print(picFolderPath)
  
  
  p1 = VlnPlot(seurat.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'dataset', pt.size = 0.1)
  p2 = VlnPlot(seurat.merge.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'dataset', pt.size = 0.1)
  
  pdf(file = file.path(picFolderPath, 
                       paste0('integrate_nFeature_nCount_ptMt_', normalizationMethod, '_beforeFiler.pdf')),
      width=10, height=5)
  print(p1)
  dev.off()
  
  pdf(file = file.path(picFolderPath, 
                       paste0('integrate_nFeature_nCount_ptMt_', normalizationMethod, '_afterFiler.pdf')),
      width=10, height=5)
  print(p2)
  dev.off()
  
  
  
  ###
  if (normalizationMethod == 'LogNormalize'){
    seurat.merge.filter <- NormalizeData(seurat.merge.filter)
    seurat.merge.filter <- FindVariableFeatures(seurat.merge.filter, selection.method = "vst", nfeatures = 2000)
    seurat.merge.filter <- ScaleData(seurat.merge.filter, vars.to.regress = c('nCount_RNA'))
  }
  if (normalizationMethod == 'SCTransform') {
    seurat.merge.filter    <- SCTransform(seurat.merge.filter, verbose = F)
  }
  
  # Idents(seurat.merge.filter) <- seurat.merge.filter@meta.data$seurat_clusters
  seurat.merge.filter    <- RunPCA(seurat.merge.filter, npcs = 50, verbose = F)
  seurat.merge.filter    <- RunUMAP(seurat.merge.filter, dims = 1:50, verbose = F)
  seurat.merge.filter    <- RunTSNE(seurat.merge.filter, dims = 1:50, verbose = F)
  seurat.merge.filter    <- FindNeighbors(seurat.merge.filter, dims = 1:50, verbose = F)
  seurat.merge.filter    <- FindClusters(seurat.merge.filter, verbose = T)
  
  
  pdf(file = file.path(picFolderPath, 
                       paste0('UMAP_tSNE_', 'merge.pdf')),
      width=15, height=6)
  print(DimPlot(seurat.merge.filter, reduction = "umap", label = TRUE) + DimPlot(seurat.merge.filter, reduction = "tsne", label = TRUE))
  dev.off()
  
  pdf(file = file.path(picFolderPath, 
                       paste0('UMAP_tSNE_', 'merge_dataset.pdf')),
      width=13, height=6)
  print(DimPlot(seurat.merge.filter, reduction = "umap", group.by = 'dataset') + DimPlot(seurat.merge.filter, reduction = "tsne", group.by = 'dataset'))
  dev.off()
  
  
  pdf(file = file.path(picFolderPath, 
                       paste0('UMAP_tSNE_', 'merge_dataset_split.pdf')),
      width=20, height=15)
  print(DimPlot(seurat.merge.filter, reduction = "umap", split.by = 'dataset') / 
          DimPlot(seurat.merge.filter, reduction = "tsne", split.by = 'dataset'))
  dev.off()
  

  saveRDS(seurat.merge.filter, file = 
            file.path('Analysis/integratedData_decontX', 'antenna_3.merge.filter.rds') )
  
}

