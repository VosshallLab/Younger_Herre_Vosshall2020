suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(scater)
  library(celda)
  library(DropletUtils)
})

theme_set(theme_cowplot())
projectFolder='/Users/tzuchiao/ProjectAnalysis/Leslie_Vosshall'
setwd(projectFolder)


#####################################

### 1. Use SoupX for removing ambient RNAs

dir.create(file.path('Analysis', '10282021_MaxPalp'), showWarnings = FALSE)
for (normalizationMethod in c('LogNormalize')) {
  dir.create(file.path('Analysis/10282021_MaxPalp', normalizationMethod), showWarnings = FALSE)
  
  dataFolder_L = c('10282021_MaxPalp')
  for (dataFolder in dataFolder_L) {
    decontxFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/1.decontx')
    dir.create(decontxFolderPath, showWarnings = FALSE)
    dataFolderPath=file.path(decontxFolderPath, dataFolder)
    dir.create(dataFolderPath, showWarnings = FALSE)
    picFolderPath=paste0(dataFolderPath, '/pic')
    dir.create(picFolderPath, showWarnings = FALSE) # Folder for the output pic
    
    
    filt.matrix <- Read10X_h5(
      paste("Analysis", '10282021_MaxPalp', "outs/filtered_feature_bc_matrix.h5", sep = '/'), 
      use.names = T)
    raw.matrix  <- Read10X_h5(
      paste("Analysis", '10282021_MaxPalp', "outs/raw_feature_bc_matrix.h5", sep = '/'), 
      use.names = T)
    
    # str(raw.matrix)
    # str(filt.matrix)
    
    srat  <- CreateSeuratObject(counts = filt.matrix)
    # Idents(srat) <- srat@meta.data$orig.ident
    srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
    

    
    pdf(file = file.path(picFolderPath, 
                         paste0('nFeature_nCount_ptMt_', dataFolder, '.pdf')),
        width=10, height=5)
    print(VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    dev.off()
    
    
    if (normalizationMethod == 'LogNormalize'){
      srat <- NormalizeData(srat)
      srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
      srat <- ScaleData(srat, vars.to.regress = c('nCount_RNA'))
    }
    
    # Idents(srat) <- srat@meta.data$seurat_clusters
    srat    <- RunPCA(srat, npcs = 50, verbose = F)
    srat    <- RunUMAP(srat, dims = 1:50, verbose = F)
    srat    <- RunTSNE(srat, dims = 1:50, verbose = F)
    srat    <- FindNeighbors(srat, dims = 1:50, verbose = F)
    srat    <- FindClusters(srat, verbose = T)
    
    pdf(file = file.path(picFolderPath, 
                         paste0('UMAP_tSNE_', dataFolder, '.pdf')),
        width=15, height=6)
    print(DimPlot(srat, reduction = "umap", label = TRUE) + DimPlot(srat, reduction = "tsne", label = TRUE))
    dev.off()
    
    pdf(file = file.path(picFolderPath, 
                         paste0('coreceptor_VlnPlot_', dataFolder, '.pdf')),
        width=15, height=8)
    print( VlnPlot(srat, features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Gr3'), ncol = 3, pt.size = 0) )
    dev.off()
    
    pdf(file = file.path(picFolderPath, 
                         paste0('coreceptor_FeaturePlot_', dataFolder, '.pdf')),
        width=15, height=10)
    print(FeaturePlot(srat, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Gr3'), ncol = 3))
    dev.off()
    
    
    ### 
    sce <- read10xCounts(paste("Analysis", dataFolder, "outs/filtered_feature_bc_matrix/", sep = '/'))
    raw.sce <- read10xCounts(paste("Analysis", dataFolder, "outs/raw_feature_bc_matrix/", sep = '/'))
    
    sce.decontX <- decontX(x = sce, background = raw.sce)
    
    # Cluster labels on UMAP
    umap <- reducedDim(sce.decontX, "decontX_UMAP")
    
    pdf(file = file.path(picFolderPath, 
                         paste0('decontX_contamination_', dataFolder, '.pdf')),
        width=15, height=6)
    print(plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
            plotDecontXContamination(sce.decontX))
    dev.off()
    
    

    
    # # 5.3 Expression of markers on UMAP
    # sce.decontX <- logNormCounts(sce.decontX)
    
    decontX_outFolder=paste0("Analysis/", dataFolder,'/decontX_test')
    if (dir.exists(decontX_outFolder)) {
      unlink( decontX_outFolder, recursive = TRUE ) #remove folder
    }
    DropletUtils:::write10xCounts(decontX_outFolder, round(decontXcounts(sce.decontX)),
                                  barcodes = colData(sce)$Barcode)
    
    
    
    #########################################
    ### Check UMAP/tSNE for adjusted matrix
    
    
    # filt.matrix <- Read10X_h5(
    #   paste("Analysis", dataFolder, "outs/filtered_feature_bc_matrix.h5", sep = '/'), 
    #   use.names = T)
    # 
    # # str(raw.matrix)
    # # str(filt.matrix)
    
    
    #########################################
    ### decontX
    sratDecontx  <- 
      Read10X(decontX_outFolder) %>%
      CreateSeuratObject(project = strsplit(dataFolder, '_')[[1]][2], min.cells = 3, min.features = 200)
    
    sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-")
    
    if (normalizationMethod == 'LogNormalize'){
      sratDecontx <- NormalizeData(sratDecontx)
      sratDecontx <- FindVariableFeatures(sratDecontx, selection.method = "vst", nfeatures = 2000)
      sratDecontx <- ScaleData(sratDecontx, vars.to.regress = c('nCount_RNA'))
    }
    
    sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
    sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:50, verbose = F)
    sratDecontx    <- FindClusters(sratDecontx, verbose = T)
    
    pdf(file = file.path(picFolderPath, 
                         paste0('UMAP_tSNE_', dataFolder, '_Decontx.pdf')),
        width=15, height=6)
    print( DimPlot(sratDecontx, reduction = "umap", label = TRUE) + DimPlot(sratDecontx, reduction = "tsne", label = TRUE) )
    dev.off()
    
    pdf(file = file.path(picFolderPath, 
                         paste0('coreceptor_VlnPlot_', dataFolder, '_Decontx.pdf')),
        width=15, height=6)
    print( VlnPlot(sratDecontx, features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Gr3'), ncol = 3, pt.size = 0) )
    dev.off()
    
    pdf(file = file.path(picFolderPath, 
                         paste0('coreceptor_FeaturePlot_', dataFolder, '_Decontx.pdf')),
        width=15, height=10)
    print(FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Gr3'), ncol = 3))
    dev.off()
    
    
    rdataFile=file.path('Analysis/10282021_MaxPalp', 'maxpalp_1.Decontx.08192021.RData')
    save.image(file = rdataFile )
    
  }
}




