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

### 2. Remove doublets
for (normalizationMethod in c('LogNormalize')) {
  
  # normalizationMethod='SCTransform'
  
  dir.create(file.path('Analysis/integratedData_decontX', normalizationMethod), showWarnings = FALSE)
  
  dataFolder_L = c('Olivia_08192021', 'Olivia_10062021')
  for (dataFolder in dataFolder_L) {
    dfFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/2.doubletfinder')
    dir.create(dfFolderPath, showWarnings = FALSE)
    dataFolderPath=file.path(dfFolderPath, dataFolder)
    dir.create(dataFolderPath, showWarnings = FALSE)
    picFolderPath=paste0(dataFolderPath, '/pic')
    dir.create(picFolderPath, showWarnings = FALSE) # Folder for the output pic
    print(dfFolderPath)
    print(dataFolderPath)
    print(picFolderPath)
    
    rdataFile=file.path('Analysis/integratedData_decontX', 
                        paste('antenna_1.Decontx', dataFolder,'RData', sep = '.'))
    load(file = rdataFile )
    
    dfFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/2.doubletfinder')
    dataFolderPath=file.path(dfFolderPath, dataFolder)
    picFolderPath=paste0(dataFolderPath, '/pic')
    
    print(dfFolderPath)
    print(dataFolderPath)
    print(picFolderPath)
    
    
    ###
    cellN=nrow(sratDecontx@meta.data)
    expDoubletRate = cellN*0.000008
    
    ## pK Identification (no  ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_scData <-
      paramSweep_v3(sratDecontx, PCs = 1:50, sct = normalizationMethod == 'SCTransform', num.cores = 4) #, num.cores = 4
    # print(head(sweep.res.list_scData))
    
    sweep.stats_scData <- summarizeSweep(sweep.res.list_scData, GT = FALSE)
    # print(head(sweep.stats_scData))
    bcmvn_scData <- find.pK(sweep.stats_scData)
    bcmvn_scData$pK <- as.numeric(as.character(bcmvn_scData$pK))
    
    # print(head(bcmvn_scData))
    pK1=bcmvn_scData$pK[bcmvn_scData$BCmetric==max(bcmvn_scData$BCmetric)]
    # print(head(pK1))
    
    p1=ggplot(data=bcmvn_scData, aes(x=pK, y=BCmetric, group=2)) +
      geom_line(color="blue")+
      geom_point()+
      geom_vline(xintercept=pK1, linetype="dashed", color = "red")+
      labs(title="pK Selection",x="pK", y = "BCmvn")+
      theme_classic()
    
    pdf(file = file.path(picFolderPath, 
                         paste0('pK_Selection_', dataFolder, '.pdf')),
        width=3.5, height=2.5)
    print(p1)
    dev.off()
    
    pK1=as.numeric(as.character( pK1 ))
    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- sratDecontx@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)   
    
    nExp_poi <- round(expDoubletRate*nrow(sratDecontx@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    # data.combined.list[[i]] <- doubletFinder_v3( data.combined.list[[i]], PCs = data.combined.list[[i]]@commands$RunUMAP.RNA.pca$dims,
    #                                              pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    
    
    if (normalizationMethod == 'LogNormalize') {
      sratDecontx <- doubletFinder_v3( sratDecontx, PCs = sratDecontx@commands$RunUMAP.RNA.pca$dims,
                                     pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    }
    if (normalizationMethod == 'SCTransform') {
      sratDecontx <- doubletFinder_v3( sratDecontx, PCs = sratDecontx@commands$RunUMAP.SCT.pca$dims,
                                     pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    }
    

    
    sratDecontx@meta.data$DoubletFinder =  sratDecontx@meta.data[,grep('DF.classifications', colnames( sratDecontx@meta.data))]
    
    print('Finish DoubletFinder')
    
    
    plot1 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'DoubletFinder')
    plot2 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'DoubletFinder')
    
    pdf(file = file.path(picFolderPath, 
                         paste0('doubletFinder_nCount_pt.mt_', dataFolder, '.pdf')),
        width=10, height=5)
    print(plot1 + plot2)
    dev.off()
    
    pdf(file = file.path(picFolderPath, 
                         paste0('doubletFinder_nFeature_nCount_ptMt_', dataFolder, '.pdf')),
        width=10, height=5)
    print(VlnPlot(sratDecontx, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = 'DoubletFinder', pt.size = 0))
    dev.off()
    
    
    sratDF <- subset(sratDecontx, subset = DoubletFinder == "Singlet")
    rm(sratDecontx)
    
    saveRDS(sratDF, file = 
              file.path('Analysis/integratedData_decontX', 
                        paste('antenna_2.DoubletFinder', dataFolder, 'rds', sep = '.')) )
    
  }

  
}


