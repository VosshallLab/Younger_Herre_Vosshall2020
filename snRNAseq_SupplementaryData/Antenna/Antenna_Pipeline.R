
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(sctransform)
  library(ComplexHeatmap)
  library(future)
  library(scales)
  library(data.table)
  library(ggplot2)
  library(stats)
  library(gridExtra)
  library(scater)
  library(celda)
  library(DropletUtils)
  library(harmony)
  library(DoubletFinder)
  library(ggrepel)
})

set.seed(96)
theme_set(theme_cowplot())

cat("\014")
rm(list = ls())   #remove all objects
dev.off()

#setwd()
projectFolder = getwd()

######################################################################################
# 1_DECONTX:
# Import matrix files, apply decontX
######################################################################################

for (dataset.name in c('antenna1_08192021', 'antenna2_10062021')) {
  
  #Import .tsv (matrix) files
  sce <- read10xCounts(paste(projectFolder, dataset.name, "filtered_feature_bc_matrix/", sep = '/'), col.names = TRUE)
  raw.sce <- read10xCounts(paste(projectFolder, dataset.name, "raw_feature_bc_matrix/", sep = '/'), col.names = TRUE)
  
  # Apply decontx
  sce.decontX <- decontX(x = sce, background = raw.sce)
  
  # PLOT: Illustration of decontamination with decontX
  umap <- reducedDim(sce.decontX, "decontX_UMAP")
  
  p1 <- plotDimReduceCluster(x = sce.decontX$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2]) + 
    plotDecontXContamination(sce.decontX)
  filename <- paste0("1a_decontX_decontaminationUMAP_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000, p1)
  rm(p1, filename)
  
  # OUTPUT: Decontx files
  DropletUtils:::write10xCounts(paste(projectFolder,dataset.name,"decontX_out",sep='/'), round(decontXcounts(sce.decontX)),
                                gene.symbol = rowData(sce)$Symbol,
                                barcodes = colData(sce)$Barcode)
  
  #Save workspace
  #filename <- paste0("1b_decontX_",dataset.name, ".rdata")
  #save.image(filename)
  
  rm(raw.sce, sce, sce.decontX, umap, filename)
  
}

rm(list = ls())   #remove all objects

######################################################################################
# 2_FILTERING
# Evaluate min.cells argument, create seurat object, normalize, cluster
######################################################################################

#....................................................................................
# EVALUATE MIN.CELLS ARGUMENT (CELLS PER GENE)
#....................................................................................
projectFolder = getwd()

for (dataset.name in c('antenna1_08192021', 'antenna2_10062021')) {

  # Determine path of decontX output files
  decontX_outFolder <- paste(projectFolder, dataset.name,"decontX_out", sep="/")
  
  # Create seurat object with no lowerbound on min.cells argument
  sratDecontx  <- 
    Read10X(decontX_outFolder) %>%
    CreateSeuratObject(project = dataset.name, min.cells = 0, min.features = 200)
  
  #Create UMI count matrix and cell counts
  countMtx <- sratDecontx@assays[["RNA"]]@counts
  nonZeroCellNDf = data.frame(nonZeroCellN=rowSums(countMtx!=0))
  
  #PLOT: Cells per gene histogram (log)
  p1=nonZeroCellNDf %>%
    ggplot(aes(log(nonZeroCellN))) + geom_histogram(binwidth = 0.1)+
    xlab('Number of cells per gene (log)') + 
    ylab('Genes')+
    geom_vline(xintercept = c(log(12)), color='red')+
    ggtitle(paste(dataset.name,": Cells per Gene, vline = 12", sep=""))
  filename <- paste0("2a_filtering_CellsPerGene_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000, p1)
  
  #Save workspace
  #filename <- paste0("2b_filtering_",dataset.name, ".rdata")
  #save.image(filename)
  
  rm(sratDecontx, p1, countMtx, nonZeroCellNDf, filename, decontX_outFolder)
  
}

rm(list = ls())   #remove all objects

#....................................................................................
# CREATE SEURAT OBJECT & NORMALIZE
#....................................................................................

projectFolder = getwd()

for (dataset.name in c('antenna1_08192021', 'antenna2_10062021')) {
  
  # Determine path of decontX output files
  decontX_outFolder <- paste(projectFolder, dataset.name,"decontX_out", sep="/")
  
  sratDecontx  <- 
    Read10X(decontX_outFolder) %>%
    CreateSeuratObject(project = dataset.name, min.cells = 12, min.features = 200)
  
  # Find percent MT
  sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-")
  
  # Normalizing
  sratDecontx    <- SCTransform(sratDecontx, verbose = F)
  
  # Run cluster analysis
  sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
  sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindClusters(sratDecontx, verbose = T)
  
  # PLOT: umap/tsne
  p1 <- DimPlot(sratDecontx, reduction = "umap", label = TRUE)
  p2 <- DimPlot(sratDecontx, reduction = "tsne", label = TRUE)
  
  g = arrangeGrob(p1,p2, ncol = 2)
  filename <- paste0("2c_filtering_UMAP-tSNE_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(g)
  
  # PLOT: coreceptor violin
  p3 <- VlnPlot(sratDecontx, features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b'), ncol = 2, pt.size = 0)
  filename <- paste0("2d_filtering_CoReceptor-Violin_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p3)
  
  # PLOT: coreceptor tsne
  p4 <- FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b'), ncol = 2)
  filename <- paste0("2e_filtering_CoReceptor-tSNE_",dataset.name,".png")
  ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p4)
  
  # Save workspace
  #filename <- paste0("2f_filtering_",dataset.name, ".rdata")
  #save.image(filename)
  
  # Save seurat file
  filename <- paste0("2g_filtering_",dataset.name, ".rds")
  saveRDS(sratDecontx, filename)
  
  rm(decontX_outFolder, sratDecontx, p1, p2, p3, p4, filename)
  
}

rm(list = ls())   #remove all objects


######################################################################################
# 3_DOUBLETFINDER
######################################################################################

#Computational intensive, run garbage collection
gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)

projectFolder = getwd()

for (dataset.name in c('antenna1_08192021', 'antenna2_10062021')) {
  print(dataset.name)

  # Load Seurat Object
  filename <- paste0("2g_filtering_",dataset.name,".rds")
  sratDecontx <- readRDS(filename)
  
  # Compute expected doublet rate
  cellN=nrow(sratDecontx@meta.data)
  expDoubletRate = (cellN*0.0008 + 0.0527)*0.01
  
  normalizationMethod='SCTransform'
  
  sweep.res.list_scData <- paramSweep_v3(sratDecontx, 
                                         PCs = 1:50, 
                                         sct = normalizationMethod == 'SCTransform', 
                                         num.cores = 1) #num.cores = 4
  #head(sweep.res.list_scData, n=10)
  
  sweep.stats_scData <- summarizeSweep(sweep.res.list_scData, GT = FALSE)
  #print(head(sweep.stats_scData))
  
  bcmvn_scData <- find.pK(sweep.stats_scData)
  bcmvn_scData$pK <- as.numeric(as.character(bcmvn_scData$pK))
  #print(head(bcmvn_scData))
  
  pK1=bcmvn_scData$pK[bcmvn_scData$BCmetric==max(bcmvn_scData$BCmetric)]
  #print(head(pK1))
  
  # PLOT: pK selection
  p1=ggplot(data=bcmvn_scData, 
            aes(x=pK, y=BCmetric, group=2)) +
    geom_line(color="blue")+
    geom_point()+
    geom_vline(xintercept=pK1, linetype="dashed", color = "red")+
    labs(title="pK Selection",x="pK", y = "BCmvn")+
    theme_classic()
  filename <- paste0("3a_doubletfinder_pkselection_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)
  
  # More doublet finder
  pK1=as.numeric(as.character( pK1 ))
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- sratDecontx@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)   
  
  nExp_poi <- round(expDoubletRate*nrow(sratDecontx@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  # data.combined.list[[i]] <- doubletFinder_v3( data.combined.list[[i]], PCs = data.combined.list[[i]]@commands$RunUMAP.RNA.pca$dims,
  #                                              pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  sratDecontx <- doubletFinder_v3( sratDecontx, PCs = sratDecontx@commands$RunUMAP.SCT.pca$dims,
                                   pN = 0.25, pK = pK1, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  
  sratDecontx@meta.data$DoubletFinder =  sratDecontx@meta.data[,grep('DF.classifications', colnames( sratDecontx@meta.data))]
  
  # PLOT: Doublet Finder graphs
  p2 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'DoubletFinder')
  p3 <- FeatureScatter(sratDecontx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'DoubletFinder')
  g = arrangeGrob(p2,p3, ncol = 2)
  filename <- paste0("3b_doubletfinder_ScatterPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,g)
  rm(filename, g)
  
  # PLOT: Violin Plots
  p4 <- VlnPlot(sratDecontx, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = 'DoubletFinder', pt.size = 0)
  filename <- paste0("3c_doubletfinder_DoubletFinder-ViolinPlots_",dataset.name,".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p4)
  
  # Save workspace
  #filename <- paste0("3d_doubletfinder_",dataset.name, ".rdata")
  #save.image(filename)
  
  # Save seurat file
  filename <- paste0("3e_doubletfinder_",dataset.name, ".rds")
  saveRDS(sratDecontx, filename)

  rm(bcmvn_scData, p1, p2, p3, p4, homotypic.prop, sweep.res.list_scData, pK1, sratDecontx, filename,
     sweep.stats_scData, annotations, cellN, expDoubletRate, nExp_poi, nExp_poi.adj, normalizationMethod)
}

rm(list = ls())   #remove all objects

######################################################################################
# 4_FEATURE & PERCENT.MT FILTERING
######################################################################################

projectFolder = getwd()

#....................................................................................
# CALCULATE FILTERING THRESHOLDS
#....................................................................................

for (dataset.name in c('antenna1_08192021', 'antenna2_10062021')) {
  #dataset.name <- 'antenna1_08192021'
  print(dataset.name)
  
  # Load RDS object
  filename <- paste0("3e_doubletfinder_",dataset.name, ".rds")
  srat <- readRDS(filename)
  
  # Subset identified singlets
  srat <- subset(srat, subset = DoubletFinder == "Singlet")
  srat$dataset <- substr(dataset.name,10,13)
  
  # Calculate feature filtering parameters
  medianFeature<- median(srat$nFeature_RNA)
  madFeature <- mad(srat$nFeature_RNA)
  
  lowerBound<- srat@meta.data %>% arrange(nFeature_RNA) %>% select(nFeature_RNA) %>% pull %>% quantile(probs=0.03)
  upperBound <- medianFeature+madFeature*3
  
  # Filter seurat object
  srat.filter=subset(srat, subset = (nFeature_RNA < medianFeature+madFeature*3) & 
                       (nFeature_RNA>= lowerBound) & (percent.mt < 5))
  
  # PLOT: Violin Plot Before/After Filtering
  p1 = VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
               group.by = 'dataset', pt.size = 0.1)
  filename <- paste0("4a_filtering_ViolinPlot_before_",srat@project.name, ".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p1)
  
  p2 = VlnPlot(srat.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
               group.by = 'dataset', pt.size = 0.1)
  filename <- paste0("4b_filtering_ViolinPlot_after_",srat@project.name, ".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p2)
  
  # PLOT: Cells per gene
  p3=ggplot(srat@meta.data, aes(nFeature_RNA)) + geom_histogram(binwidth = 50) +
    geom_vline(xintercept = c(medianFeature), color='black') + coord_cartesian(xlim = c(0,4000)) +
    geom_vline(xintercept = c(lowerBound, upperBound), color='red')+
    ggtitle(paste0(srat@project.name,": median=",medianFeature,
                   ", lowerBound=",lowerBound,", upperBound=",upperBound,
                   ", Cells pre-filter=",nrow(srat@meta.data),
                   ", Cells post-filter=",nrow(srat.filter@meta.data)))
  filename <- paste0("4c_filtering_FeatureHistogram_",srat@project.name, ".pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p3)
  
  # Cluster
  srat.filter    <- RunPCA(srat.filter, npcs = 50, verbose = F)
  srat.filter    <- RunUMAP(srat.filter, dims = 1:50, verbose = F)
  srat.filter    <- RunTSNE(srat.filter, dims = 1:50, verbose = F)
  srat.filter    <- FindNeighbors(srat.filter, dims = 1:50, verbose = F)
  srat.filter    <- FindClusters(srat.filter, verbose = T)
  
  # PLOT: tSNE/UMAP
  p4 = DimPlot(srat.filter, reduction = "tsne", label = TRUE) +
    DimPlot(srat.filter, reduction = "umap", label = TRUE)
  filename <- paste0("4d_filtering_",srat.filter@project.name,"_tsne-umap.pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000, p4)
  
  # Save seurat file
  filename <- paste0("4e_filtering_",srat@project.name, "_clustered.rds")
  print("Saving RDS...")
  saveRDS(srat.filter, filename)
  
  rm(p1, p2, p3, p4, srat, srat.filter, lowerBound, upperBound, medianFeature, madFeature, filename)
  
}

rm(list = ls())   #remove all objects


######################################################################################
# 5_QUALITY CONTROL
######################################################################################

seurat.0819 <- readRDS("4e_filtering_antenna1_08192021_clustered.rds")
seurat.1006 <- readRDS("4e_filtering_antenna2_10062021_clustered.rds")

DefaultAssay(seurat.0819) <- "SCT"
DefaultAssay(seurat.1006) <- "SCT"

#....................................................................................
# PLOT: UMAP/TSNE
#....................................................................................

title1 <- "antenna1_0819"
title2 <- "antenna2_1006"

p01 = DimPlot(seurat.0819, reduction = "tsne", label = TRUE) + ggtitle(title1) + NoLegend()
p02 = DimPlot(seurat.0819, reduction = "umap", label = TRUE) + ggtitle(title1) + NoLegend()
p03 = DimPlot(seurat.1006, reduction = "tsne", label = TRUE) + ggtitle(title2) + NoLegend()
p04 = DimPlot(seurat.1006, reduction = "umap", label = TRUE) + ggtitle(title2) + NoLegend()

g = arrangeGrob(p01, p02, p03, p04, ncol = 2)
filename <- "5a_qualitycontrol_tsne-umap.pdf"
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000, g)

rm(p01, p02, p03, p04, g, filename)

#....................................................................................
# CALCULATE LOG(nCOUNT_RNA)
#....................................................................................

log_count = log10(seurat.0819@meta.data$nCount_RNA)
seurat.0819@meta.data$log_count = log_count
rm(log_count)

log_count = log10(seurat.1006@meta.data$nCount_RNA)
seurat.1006@meta.data$log_count = log_count
rm(log_count)

#....................................................................................
# FEATURE PLOT: FEATURE COUNT, TOTAL COUNT, PCT MT
#....................................................................................

p1 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T, features = 'nFeature_RNA')
p2 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T, features = 'log_count')
p3 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T, features = 'percent.mt')
p4 = FeaturePlot(seurat.0819, reduction = 'umap', label=T, features = 'nFeature_RNA')
p5 = FeaturePlot(seurat.0819, reduction = 'umap', label=T, features = 'log_count')
p6 = FeaturePlot(seurat.0819, reduction = 'umap', label=T, features = 'percent.mt')

g = arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3)
filename <- paste0("5b_qualitycontrol_feature-logcount-mt_",seurat.0819@project.name,".pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 4000,g)
rm(p1, p2, p3, p4, p5, p6, g, filename)

p1 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T, features = 'nFeature_RNA')
p2 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T, features = 'log_count')
p3 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T, features = 'percent.mt')
p4 = FeaturePlot(seurat.1006, reduction = 'umap', label=T, features = 'nFeature_RNA')
p5 = FeaturePlot(seurat.1006, reduction = 'umap', label=T, features = 'log_count')
p6 = FeaturePlot(seurat.1006, reduction = 'umap', label=T, features = 'percent.mt')

g = arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3)
filename <- paste0("5b_qualitycontrol_feature-logcount-mt_",seurat.1006@project.name,".pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 4000,g)
rm(p1, p2, p3, p4, p5, p6, g, filename)

#....................................................................................
# VIOLIN PLOT: FEATURE COUNT, TOTAL COUNT, PCT MT
#....................................................................................

p1 = VlnPlot(seurat.0819, pt.size=0, features = 'nFeature_RNA') + geom_boxplot()+
  ggtitle(paste0(title1,', nFeature_RNA')) + theme(legend.position = 'none')
p2 = VlnPlot(seurat.0819, pt.size=0, features = 'log_count') + geom_boxplot()+
  ggtitle(paste0(title1,', log_count')) + theme(legend.position = 'none')
p3 = VlnPlot(seurat.0819, pt.size=0, features = 'percent.mt') + geom_boxplot()+
  ggtitle(paste0(title1,', percent.mt')) + theme(legend.position = 'none')
p4 = VlnPlot(seurat.1006, pt.size=0, features = 'nFeature_RNA') + geom_boxplot()+
  ggtitle(paste0(title2,', nFeature_RNA')) + theme(legend.position = 'none')
p5 = VlnPlot(seurat.1006, pt.size=0, features = 'log_count') + geom_boxplot()+
  ggtitle(paste0(title2,', log_count')) + theme(legend.position = 'none')
p6 = VlnPlot(seurat.1006, pt.size=0, features = 'percent.mt') + geom_boxplot()+
  ggtitle(paste0(title2,', percent.mt')) + theme(legend.position = 'none')

g = arrangeGrob(p1,p2,p3,p4,p5,p6, ncol = 3)
filename <- paste0("5c_qualitycontrol_feature-logcount-mt_vlnplots.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 10000, height = 4000,g)
rm(p1, p2, p3, p4, p5, p6, g, filename)

#....................................................................................
# FEATURE PLOT: NEURAL MARKERS
#....................................................................................

p17 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T, features = c('LOC5564848'))+
  ggtitle(paste0(title2,', LOC5564848'))
p18 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T,  features = c('LOC5570381'))+
  ggtitle(paste0(title2,', LOC5570381'))
p19 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T,  features = c('LOC5565901'))+
  ggtitle(paste0(title2,', LOC5565901'))
p20 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T,  features = c('LOC5570204'))+
  ggtitle(paste0(title2,', LOC5570204'))

p21 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('LOC5564848'))+
  ggtitle(paste0(title2,', LOC5564848'))
p22 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('LOC5570381'))+
  ggtitle(paste0(title2,', LOC5570381'))
p23 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('LOC5565901'))+
  ggtitle(paste0(title2,', LOC5565901'))
p24 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('LOC5570204'))+
  ggtitle(paste0(title2,', LOC5570204'))

g = arrangeGrob(p17, p18, p19, p20, p21, p22, p23, p24, ncol = 4)
filename <- paste0("5d_qualitycontrol_neuralmarkers_featureplots.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 10000, height = 4000,g)
rm(p17, p18, p19, p20, p21, p22, p23, p24, g, filename)

#....................................................................................
# VIOLIN PLOT: NEURAL MARKERS
#....................................................................................

p1 = VlnPlot(seurat.0819, pt.size=0, features = 'LOC5564848') + 
  ggtitle(paste0(title1,', LOC5564848')) + theme(legend.position = 'none')
p2 = VlnPlot(seurat.0819, pt.size=0, features = 'LOC5570381') + 
  ggtitle(paste0(title1,', LOC5570381')) + theme(legend.position = 'none')
p3 = VlnPlot(seurat.0819, pt.size=0, features = 'LOC5565901') + 
  ggtitle(paste0(title1,', LOC5565901')) + theme(legend.position = 'none')
p4 = VlnPlot(seurat.0819, pt.size=0, features = 'LOC5570204') + 
  ggtitle(paste0(title1,', LOC5570204')) + theme(legend.position = 'none')

p5 = VlnPlot(seurat.1006, pt.size=0, features = 'LOC5564848') + 
  ggtitle(paste0(title2,', LOC5564848')) + theme(legend.position = 'none')
p6 = VlnPlot(seurat.1006, pt.size=0, features = 'LOC5570381') + 
  ggtitle(paste0(title2,', LOC5570381')) + theme(legend.position = 'none')
p7 = VlnPlot(seurat.1006, pt.size=0, features = 'LOC5565901') + 
  ggtitle(paste0(title2,', LOC5565901')) + theme(legend.position = 'none')
p8 = VlnPlot(seurat.1006, pt.size=0, features = 'LOC5570204') + 
  ggtitle(paste0(title2,', LOC5570204')) + theme(legend.position = 'none')

g = arrangeGrob(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 4)
filename <- paste0("5e_qualitycontrol_neuralmarkers_vlnplots.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 12000, height = 4000,g)
rm(p1,p2,p3,p4,p5,p6,p7,p8, g, filename)

#....................................................................................
# FEATURE PLOT: CORECEPTORS
#....................................................................................

p17 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T, features = c('Orco'))+
  ggtitle(paste0(title2,', Orco'))
p18 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T,  features = c('Ir25a'))+
  ggtitle(paste0(title2,', Ir25a'))
p19 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T,  features = c('Ir76b'))+
  ggtitle(paste0(title2,', Ir76b'))
p20 = FeaturePlot(seurat.0819, reduction = 'tsne', label=T,  features = c('Ir8a'))+
  ggtitle(paste0(title2,', Ir8a'))

p21 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('Orco'))+
  ggtitle(paste0(title2,', Orco'))
p22 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('Ir25a'))+
  ggtitle(paste0(title2,', Ir25a'))
p23 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('Ir76b'))+
  ggtitle(paste0(title2,', Ir76b'))
p24 = FeaturePlot(seurat.1006, reduction = 'tsne', label=T,  features = c('Ir8a'))+
  ggtitle(paste0(title2,', Ir8a'))

g = arrangeGrob(p17,p18,p19,p20,p21, p22, p23, p24, ncol = 4)
filename <- paste0("5f_qualitycontrol_coreceptors.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 10000, height = 4000,g)
rm(g, filename)

rm(list = ls())   #remove all objects


######################################################################################
# MERGE
######################################################################################

#....................................................................................
# LOAD OBJECTS
#....................................................................................

seurat.0819 <- readRDS("4e_filtering_antenna1_08192021_clustered.rds")
seurat.1006 <- readRDS("4e_filtering_antenna2_10062021_clustered.rds")

DefaultAssay(seurat.0819) <- "SCT"
DefaultAssay(seurat.1006) <- "SCT"

seurat.0819$dataset <- '0819'
seurat.1006$dataset <- '1006'


#....................................................................................
# MERGE
#....................................................................................

seurat.merge.filter =  merge(seurat.0819 , y =seurat.1006,
                             add.cell.ids = c("0819", "1006"), project = "antenna")

# Renormalize
seurat.merge.filter    <- SCTransform(seurat.merge.filter, verbose = T)

# Recluster
seurat.merge.filter    <- RunPCA(seurat.merge.filter, npcs = 50, verbose = F)
seurat.merge.filter    <- RunUMAP(seurat.merge.filter, dims = 1:50, verbose = F)
seurat.merge.filter    <- RunTSNE(seurat.merge.filter, dims = 1:50, verbose = F)
seurat.merge.filter    <- FindNeighbors(seurat.merge.filter, dims = 1:50, verbose = F)
seurat.merge.filter    <- FindClusters(seurat.merge.filter, verbose = T)

#....................................................................................
# PLOTS
#....................................................................................

# UMAP/tSNE
p1 <- DimPlot(seurat.merge.filter, reduction = "umap", label = TRUE) + 
  DimPlot(seurat.merge.filter, reduction = "tsne", label = TRUE)
filename <- paste0("6a_merge_tsne-umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)

# UMAP/tSNE - Colored by batch
p2 <- DimPlot(seurat.merge.filter, reduction = "umap", group.by = 'dataset') + 
  DimPlot(seurat.merge.filter, reduction = "tsne", group.by = 'dataset')
filename <- paste0("6b_merge_tsne-umap_coloredbybatch.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p2)

# UMAP/tSNE - batch v batch
p3 <- DimPlot(seurat.merge.filter, reduction = "umap", split.by = 'dataset') /
  DimPlot(seurat.merge.filter, reduction = "tsne", split.by = 'dataset')
filename <- paste0("6c_merge_tsne-umap_splitbybatch.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p3)

#....................................................................................
# SAVE SEURAT OBJECT
#....................................................................................

filename <- "6d_merge.rds"
print("Saving RDS...")
saveRDS(seurat.merge.filter, filename)

rm(list = ls())   #remove all objects

######################################################################################
# HARMONY
######################################################################################

seurat.merge.filter <- readRDS("6d_merge.rds")
seurat.merge.filter    <- RunPCA(seurat.merge.filter, npcs = 50, verbose = F)

#....................................................................................
# HARMONY
#....................................................................................

reduction='harmony'
pc1='harmony_1'
DefaultAssay(seurat.merge.filter) <- "SCT"

seurat.combined <- RunHarmony(seurat.merge.filter, "dataset", 
                              plot_convergence = TRUE, assay.use = "SCT")

harmony_embeddings <- Embeddings(seurat.combined, 'harmony')
# harmony_embeddings[1:5, 1:5]

seurat.combined    <- seurat.combined %>%
  RunUMAP(dims = 1:50, verbose = F, reduction = reduction) %>%
  RunTSNE(dims = 1:50, verbose = F, reduction = reduction) %>%
  FindNeighbors(dims = 1:50, verbose = F, reduction = reduction) %>%
  FindClusters(verbose = T)

# Check that batch correction worked
p1 <- DimPlot(seurat.merge.filter, reduction = "tsne", group.by = 'dataset')
p2 <- DimPlot(seurat.combined, reduction = "tsne", group.by = 'dataset')
p3 <- DimPlot(seurat.merge.filter, reduction = "umap", group.by = 'dataset')
p4 <- DimPlot(seurat.combined, reduction = "umap", group.by = 'dataset')

g = arrangeGrob(p1,p2,p3,p4, ncol = 2)
filename <- paste0("7a_harmony_batchcorrection.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
rm(p1, p2, p3, p4, g, filename)

#....................................................................................
# PCs WITH CORRECTION
#....................................................................................

p5 <- DimPlot(object = seurat.merge.filter, reduction = 'pca', pt.size = .1, group.by = "dataset")+
  VlnPlot(object = seurat.merge.filter, features = 'PC_1', group.by = "dataset", pt.size = .1)
filename <- paste0("7b_harmony_pcsbefore.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p5)

p6 <- DimPlot(object = seurat.combined, reduction = 'pca', pt.size = .1, group.by = "dataset")+
  VlnPlot(object = seurat.combined, features = 'PC_1', group.by = "dataset", pt.size = .1)
filename <- paste0("7c_harmony_pcsafter.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p6)

#....................................................................................
# PLOTS
#....................................................................................

# UMAP/tSNE
p7 <- DimPlot(seurat.combined, reduction = "umap", label = TRUE) + DimPlot(seurat.combined, reduction = "tsne", label = TRUE)
filename <- paste0("7d_harmony_tsne-umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p7)

# UMAP/tSNE - batch v batch
p8 <- DimPlot(seurat.merge.filter, reduction = "tsne", split.by = 'dataset') /
  DimPlot(seurat.merge.filter, reduction = "umap", split.by = 'dataset')
filename <- paste0("7e_harmony_before-splitbatch.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p8)

p9 <- DimPlot(seurat.combined, reduction = "tsne", split.by = 'dataset') /
  DimPlot(seurat.combined, reduction = "umap", split.by = 'dataset')
filename <- paste0("7f_harmony_after-splitbatch.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p9)

#QC Feature check (on uncorrected)
log_count = log10(seurat.merge.filter@meta.data$nCount_RNA)
seurat.merge.filter@meta.data$log_count = log_count
rm(log_count)

pointsize <- 1
p09 = FeaturePlot(seurat.merge.filter, reduction = 'tsne', label=T, features = 'nFeature_RNA', pt.size = pointsize)
p10 = FeaturePlot(seurat.merge.filter, reduction = 'tsne', label=T, features = 'log_count', pt.size = pointsize)
p11 = FeaturePlot(seurat.merge.filter, reduction = 'umap', label=T, features = 'nFeature_RNA', pt.size = pointsize)
p12 = FeaturePlot(seurat.merge.filter, reduction = 'umap', label=T, features = 'log_count', pt.size = pointsize)

g = arrangeGrob(p09, p10, p11, p12, ncol = 2)
filename <- paste0("7g_harmony_QCfeaturecheck-beforeBC.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 5000,g)
rm(p09, p10, p11, p12, filename, g)

#QC Feature check (on corrected)
log_count = log10(seurat.combined@meta.data$nCount_RNA)
seurat.combined@meta.data$log_count = log_count
rm(log_count)

p09 = FeaturePlot(seurat.combined, reduction = 'tsne', label=T, features = 'nFeature_RNA', pt.size = pointsize)
p10 = FeaturePlot(seurat.combined, reduction = 'tsne', label=T, features = 'log_count', pt.size = pointsize)
p11 = FeaturePlot(seurat.combined, reduction = 'umap', label=T, features = 'nFeature_RNA', pt.size = pointsize)
p12 = FeaturePlot(seurat.combined, reduction = 'umap', label=T, features = 'log_count', pt.size = pointsize)

g = arrangeGrob(p09, p10, p11, p12, ncol = 2)
filename <- paste0("7h_harmony_QCfeaturecheck-afterBC.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 5000,g)
rm(p09, p10, p11, p12, filename, g, pointsize)

#....................................................................................
# SAVE SEURAT OBJECT
#....................................................................................

filename <- "7i_harmony.rds"

print("Saving RDS...")
saveRDS(seurat.combined, filename)

rm(list = ls())   #remove all objects


######################################################################################
# IDENTIFY NEURON CLUSTERS / HARMONY
######################################################################################

seurat.combined <- readRDS("7i_harmony.rds")

reductionMethod='harmony'
DefaultAssay(seurat.combined) <- "SCT"
seurat.combined <- FindClusters(seurat.combined, resolution = 1, reduction = reductionMethod)

seurat.combined@meta.data$batch <- apply(seurat.combined@meta.data, 1,FUN = function(input) {
  dataDate = as.character(input['dataset'])
  # print(dataDate)
  if (dataDate == '0819') {return('batch1')}
  else {return('batch2')}
})

filename <- "7j_harmony_AllCellFinal.rds" #Uploaded on Zenodo as "SeuratObject1_Antenna_mergedBatches_AllCells.rds"
saveRDS(seurat.combined, filename)

reduction='harmony'
pc1='harmony_1'
DefaultAssay(seurat.combined) <- "SCT"

#....................................................................................
# PLOTS
#....................................................................................

DefaultAssay(seurat.combined) <- "SCT"

# UMAP/tSNE
p1 <- DimPlot(seurat.combined, reduction = "umap", label = TRUE) + 
  DimPlot(seurat.combined, reduction = "tsne", label = TRUE)
filename <- paste0("8a_IDneurons_tsne-umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)

# UMAP/tSNE - Colored by batch
p2 <- DimPlot(seurat.combined, reduction = "umap", group.by = 'dataset') + 
  DimPlot(seurat.combined, reduction = "tsne", group.by = 'dataset')
filename <- paste0("8b_IDneurons_coloredbatch.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p2)

# UMAP/tSNE - batch v batch
p3 <- DimPlot(seurat.combined, reduction = "umap", split.by = 'dataset') /
  DimPlot(seurat.combined, reduction = "tsne", split.by = 'dataset')
filename <- paste0("8c_IDneurons_splitbatch.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 8000,p3)

# UMAP/tSNE - coreceptors
p4 <- FeaturePlot(seurat.combined, reduction = 'tsne', 
                  features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Or82', 'LOC5575210'), ncol = 2)
filename <- paste0("8d_IDneurons_coreceptors.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 10000,p4)

# UMAP/tSNE - neuron clusters
p5 <- FeaturePlot(seurat.combined, reduction = 'tsne', 
                  features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'), ncol = 2)
filename <- paste0("8e_IDneurons_neuronmarkers.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 8000,p5)

# DotPlot neuron markers
p6 <- DotPlot(seurat.combined, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'))+ 
  scale_y_discrete(limits = rev)
filename <- paste0("8f_IDneurons_dotplot.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 4000,p6)

#....................................................................................
# NEURAL MARKERS
#....................................................................................

neuronM.plot <- DotPlot(seurat.combined, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204')) +
  scale_y_discrete(limits = rev)
neuronM.percentage  = neuronM.plot$data

neuronM.percentage.wider <- neuronM.percentage %>%
  pivot_wider(names_from = features.plot, values_from = c(avg.exp, pct.exp, avg.exp.scaled))

neuronM.percentage.wider <- neuronM.percentage %>%
  pivot_wider(names_from = features.plot, values_from = c(avg.exp, pct.exp, avg.exp.scaled))

neuronM.percentage.wider$ratio50N <- apply(neuronM.percentage.wider, 1, function(input){
  # print( as.numeric(input[6:9])  )
  ratio50counts=0
  for (ratio in as.numeric(input[6:9])) {
    if (ratio > 50) {ratio50counts= ratio50counts + 1}
  }
  # if (ratio50counts >=3) {print(ratio50counts)}
  return( ratio50counts )
})

neuronClusterIDs <- neuronM.percentage.wider %>% filter(ratio50N >= 3) %>% pull(id) %>% as.character()
antenna.neuron <- subset(seurat.combined, subset = seurat_clusters %in% neuronClusterIDs)

#....................................................................................
# PLOT
#....................................................................................

# UMAP/tSNE
p1 <- DimPlot(antenna.neuron, reduction = "umap", label = TRUE) + 
  DimPlot(antenna.neuron, reduction = "tsne", label = TRUE)
filename <- paste0("8g_IDneurons_neuralclusters_tsne-umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)

# DotPlot neuron markers
p6 <- DotPlot(antenna.neuron, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'))+ 
  scale_y_discrete(limits = rev)
filename <- paste0("8h_IDneurons_neuralclusters_dotplot.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 4000,p6)

#....................................................................................
# VIOLIN PLOT: FEATURE COUNT, TOTAL COUNT, PCT MT
#....................................................................................

p1 = VlnPlot(antenna.neuron, pt.size=0, features = 'nFeature_RNA') + 
  geom_boxplot()+ theme(legend.position = 'none')
p2 = VlnPlot(antenna.neuron, pt.size=0, features = 'log_count') + 
  geom_boxplot() + theme(legend.position = 'none')
p3 = VlnPlot(antenna.neuron, pt.size=0, features = 'percent.mt') + 
  geom_boxplot()+ theme(legend.position = 'none')

g = arrangeGrob(p1,p2,p3, ncol = 3)
filename <- paste0("8i_IDneurons_neuralclusters_QCviolin.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 10000, height = 4000,g)
rm(p1, p2, p3, p4, p5, p6, g, filename)

#....................................................................................
# FEATURE PLOT: FEATURE COUNT, TOTAL COUNT, PCT MT
#....................................................................................

p1 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'nFeature_RNA')
p2 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'log_count')
p3 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'percent.mt')
p4 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'nFeature_RNA')
p5 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'log_count')
p6 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'percent.mt')

g = arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3)
filename <- paste0("8j_IDneurons_neuralclusters_feature-logcount-mt.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 4000,g)
rm(p1, p2, p3, p4, p5, p6, g, filename)

#....................................................................................
# FEATURE PLOT: NEURAL MARKERS
#....................................................................................

p21 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5564848'))
p22 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5570381'))
p23 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5565901'))
p24 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5570204'))

g = arrangeGrob(p21, p22, p23, p24, ncol = 2)
filename <- paste0("8k_IDneurons_neuralclusters_neuralmarkers_featureplots.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
rm(p21, p22, p23, p24, g, filename)

saveRDS(antenna.neuron, '8l_IDneurons_neuralclusters.rds') 

#....................................................................................
# RECLUSTER & PLOT
#....................................................................................

# ID'd ambiguously neural cluster
cluster.list <- levels(Idents(antenna.neuron))
clusters_to_remove <- c(grep("\\b3\\b",cluster.list), grep("45",cluster.list))
cluster.list <- cluster.list[-c(clusters_to_remove)] #remove clusters 3 & 45
antenna.neuron <- subset(antenna.neuron,idents = c(cluster.list))
rm(cluster.list)

# UMAP/tSNE
p1 <- DimPlot(antenna.neuron, reduction = "umap", label = TRUE) + 
  DimPlot(antenna.neuron, reduction = "tsne", label = TRUE)
filename <- paste0("8m_IDneurons_rmclusts_tsne-umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)

# Recluster
reduction = 'harmony'
antenna.neuron    <- antenna.neuron %>%
  RunUMAP(dims = 1:50, verbose = T, reduction = reduction) %>%
  RunTSNE(dims = 1:50, verbose = T, reduction = reduction) %>%
  FindNeighbors(dims = 1:50, verbose = T, reduction = reduction) %>%
  FindClusters(verbose = T)

#....................................................................................
# PLOT
#....................................................................................

# UMAP/tSNE
p1 <- DimPlot(antenna.neuron, reduction = "umap", label = TRUE) + 
  DimPlot(antenna.neuron, reduction = "tsne", label = TRUE)
filename <- paste0("8n_IDneurons_rmclusts_tsne-umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)

# DotPlot neuron markers
p6 <- DotPlot(antenna.neuron, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'))+ 
  scale_y_discrete(limits = rev)
filename <- paste0("8o_IDneurons_rmclusts_neuralclusters_dotplot.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 4000,p6)

#....................................................................................
# VIOLIN PLOT: FEATURE COUNT, TOTAL COUNT, PCT MT
#....................................................................................

p1 = VlnPlot(antenna.neuron, pt.size=0, features = 'nFeature_RNA') + 
  geom_boxplot()+ theme(legend.position = 'none')
p2 = VlnPlot(antenna.neuron, pt.size=0, features = 'log_count') + 
  geom_boxplot() + theme(legend.position = 'none')
p3 = VlnPlot(antenna.neuron, pt.size=0, features = 'percent.mt') + 
  geom_boxplot()+ theme(legend.position = 'none')

g = arrangeGrob(p1,p2,p3, ncol = 3)
filename <- paste0("8p_IDneurons_rmclusts_QCviolin.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 10000, height = 4000,g)
rm(p1, p2, p3, p4, p5, p6, g, filename)

#....................................................................................
# FEATURE PLOT: FEATURE COUNT, TOTAL COUNT, PCT MT
#....................................................................................

p1 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'nFeature_RNA')
p2 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'log_count')
p3 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'percent.mt')
p4 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'nFeature_RNA')
p5 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'log_count')
p6 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'percent.mt')

g = arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3)
filename <- paste0("8q_IDneurons_rmclusts_feature-logcount-mt.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 4000,g)
rm(p1, p2, p3, p4, p5, p6, g, filename)

#....................................................................................
# FEATURE PLOT: NEURAL MARKERS
#....................................................................................

p21 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5564848'))
p22 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5570381'))
p23 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5565901'))
p24 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T,  features = c('LOC5570204'))

g = arrangeGrob(p21, p22, p23, p24, ncol = 2)
filename <- paste0("8r_IDneurons_rmclusts_neuralmarkers_featureplots.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
rm(p21, p22, p23, p24, g, filename)

saveRDS(antenna.neuron, '8s_IDneurons_rmclusts.rds')
saveRDS(neuronClusterIDs, '8t_neuronClusterIDs.rds')
rm(list = ls())   #remove all objects


######################################################################################
# TEST CLUSTERS / RESOLUTION TEST
######################################################################################

#....................................................................................
# RESOLUTION TEST
#....................................................................................

projectFolder= getwd()
getwd()
antenna.neuron <- readRDS('8s_IDneurons_rmclusts.rds')
dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- "SCT"

filename <- "resolution_test_harmony_rmclusts"
dir.create(filename)
setwd(paste(projectFolder,filename,sep='/'))
getwd()

#....................................................................................
# SET STATIC GENE LIST
#....................................................................................

dat <- antenna.neuron
dat <- FindClusters(antenna.neuron, verbose = T, resolution=1.0)

numberClusters <- length(levels(Idents(dat)))
clusterList <- levels(Idents(dat))
heatmap.genes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a")


suppressMessages({ 
for (i in 1:numberClusters) {
  u <- clusterList[i]
  
  cluster.subset <- subset(dat,idents = u)
  #print(paste0("cluster ",toString(u)," has ",dim(cluster.subset@meta.data)[1]," cells"))
  
  cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)                           #get avg expression of all genes in cluster      
  cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
  cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
  cl.ae <- as.data.frame(cl.ae)                             #convert to df      
  cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
  cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
  cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
  
  cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
  cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
  cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
  cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
  #cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
  #print(head(cl.ae$Gene, 5))
  
  heatmap.genes <- append(heatmap.genes, head(cl.ae$Gene, 5))
  
}
}) 

heatmap.genes <- unique(heatmap.genes)
res1_staticgenelist <- heatmap.genes

rm(u, cluster.subset, i, numberClusters, cl.ae, clusterList, heatmap.genes, dat)

#....................................................................................
theme_set(theme_grey()) #white backgrounds
test <- c('0.8', '1.0', '1.2', '1.5', '1.7', '2.0', '2.5', '3.0')

for (resolution in test) {
  #resolution <- '0.8'
  
  antenna.neuron <- FindClusters(antenna.neuron, verbose = T, resolution=as.double(resolution))
  DefaultAssay(antenna.neuron) <- "SCT"
  
  pipeline <- paste0("harmony_rmclusts_",resolution)
  
  DefaultAssay(antenna.neuron) <- "SCT"
  filetype=".png"
  
  pointsize = 1
  p1 <- DimPlot(antenna.neuron, reduction = "tsne", label = TRUE, pt.size = pointsize) +NoLegend()
  filename <- paste0("9a_heatmaps_",pipeline,"_tsne",filetype)
  ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p1)
  rm(p1, filename)
  
  heatmap.title1 = paste("All clusters, Method: ",pipeline,
                         "\n",toString(dim(antenna.neuron@meta.data)[1])," total neurons",sep="")
  
  #p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = antenna.neuron, 
  #                features = res1_staticgenelist, label=TRUE, group.bar.height =0.01, 
  #                draw.lines=FALSE, lines.width=4, size=4) + 
  #  theme(axis.text.y = element_text(size = 10)) +
  #  scale_fill_gradientn(colors = c("black", "yellow", "red"), a.value = "black")+ ggtitle(heatmap.title1)
  #filename <- paste0("11-2h_heatmaps_",pipeline,"_allclusters_staticgenelist",filetype)
  #ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 3000,p1)
  
  rm(dat, filename, heatmap.title1, p1)
  
  
  ######################################################################################
  # PRELIM VISUALIZATIONS
  ######################################################################################
  
  # UMAP/tSNE
  #p1 <- DimPlot(antenna.neuron, reduction = "umap", label = TRUE, pt.size = pointsize) + 
  #  DimPlot(antenna.neuron, reduction = "tsne", label = TRUE, pt.size = pointsize)
  #filename <- paste0("11-2a_heatmaps_",pipeline,"_tsne-umap",filetype)
  #ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)
  #
  #
  ## UMAP/tSNE - Colored by batch
  #p2 <- DimPlot(antenna.neuron, reduction = "umap", group.by = 'dataset', pt.size = pointsize) + 
  #  DimPlot(antenna.neuron, reduction = "tsne", group.by = 'dataset', pt.size = pointsize)
  #filename <- paste0("11-2b_heatmaps_",pipeline,"_coloredbatch",filetype)
  #ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 4000,p2)
  #
  ## UMAP/tSNE - batch v batch
  #p3 <- DimPlot(antenna.neuron, reduction = "umap", split.by = 'dataset', pt.size = pointsize) /
  #  DimPlot(antenna.neuron, reduction = "tsne", split.by = 'dataset', pt.size = pointsize)
  #filename <- paste0("11-2c_heatmaps_",pipeline,"_splitbatch",filetype)
  #ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 8000,p3)
  #
  ## UMAP/tSNE - coreceptors
  #p4 <- FeaturePlot(antenna.neuron, reduction = 'tsne', 
  #                  features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Or82', 'LOC5575210'), 
  #                  ncol = 3, pt.size = pointsize)
  #filename <- paste0("11-2d_heatmaps_",pipeline,"_coreceptors",filetype)
  #ggsave(filename, limitsize = FALSE, units = "px", width = 12000, height = 8000,p4)
  #
  ## UMAP/tSNE - neuron clusters
  #p5 <- FeaturePlot(antenna.neuron, reduction = 'tsne', 
  #                  features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'), 
  #                  ncol = 2, pt.size = pointsize)
  #filename <- paste0("11-2e_heatmaps_",pipeline,"_neuronmarkers",filetype)
  #ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 8000,p5)
  #
  #rm(p1,p2,p3,p4,p5)
  #
  ##....................................................................................
  ## PLOT: FEATURE COUNT VS TOTAL COUNT
  ##....................................................................................
  #
  #log_count = log10(antenna.neuron@meta.data$nCount_RNA)
  #antenna.neuron@meta.data$log_count = log_count
  #rm(log_count)
  #
  #p09 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'nFeature_RNA', pt.size = pointsize)
  #p10 = FeaturePlot(antenna.neuron, reduction = 'tsne', label=T, features = 'log_count', pt.size = pointsize)
  #p11 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'nFeature_RNA', pt.size = pointsize)
  #p12 = FeaturePlot(antenna.neuron, reduction = 'umap', label=T, features = 'log_count', pt.size = pointsize)
  #
  #g = arrangeGrob(p09, p10, p11, p12, ncol = 2)
  #filename <- paste0("11-2f_heatmaps_",pipeline,"_QC_feature-v-count",filetype)
  #ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 5000,g)
  #
  #rm(p1, p2, p3, p4, p5, p09, p10, p11, p12, filename, g)
  
  
  
  ######################################################################################
  # HEATMAP OF ALL CLUSTERS
  ######################################################################################
  
  
  dat <- antenna.neuron
  
  numberClusters <- length(levels(Idents(dat)))
  clusterList <- levels(Idents(dat))
  heatmap.genes <- c('LOC5575210',"Orco", "Ir25a", "Ir76b", "Ir8a")
  
  suppressMessages({
  for (i in 1:numberClusters) {
    u <- clusterList[i]
    
    cluster.subset <- subset(dat,idents = u)
    #print(paste0("cluster ",toString(u)," has ",dim(cluster.subset@meta.data)[1]," cells"))
    
    cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)                           #get avg expression of all genes in cluster      
    cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    #cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #print(head(cl.ae$Gene, 5))
    
    heatmap.genes <- append(heatmap.genes, head(cl.ae$Gene, 5))
    
  }
  })
  
  heatmap.genes <- unique(heatmap.genes)
  
  rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)
  
  heatmap.title1 = paste("All clusters, Method: ",pipeline,
                         "\n",toString(dim(antenna.neuron@meta.data)[1])," total neurons",sep="")
  
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = antenna.neuron, 
                  features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                  draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)
  #p1
  
  filename <- paste0("9h_heatmaps_",pipeline,"_allclusters",filetype)
  ggsave(filename, limitsize = FALSE, units = "px", width = 8000, height = 3000,p1)
  
  rm(dat, filename, heatmap.title1, p1)
  
  
}

#....................................................................................
# SAVE FINAL RDS
#....................................................................................

##
# res0.8 --> res1.0
# 0      --> 1, 26       valid
# 5      --> 6, 25       valid
# 2      --> 4, 19       valid


# res1.0 --> res1.2
# no change

# res1.2 --> res1.5
# 1      --> 1, 24       valid
# 2      --> 7, 14       not valid

# res1.5 --> res1.7
# 17     --> 29, 27      valid

# res1.7 --> res2.0
# 1, 23  --> 7, 6        valid
# 2      --> 4, 23       not valid

# res2.0 --> res2.5
# 0      --> 0, 14       valid
# 6      --> 15, 25      valid
# 11     --> 17, 32      unsure
# 22     --> 33, 34      valid

# res2.5 --> res3.0
# 9      --> 18, 29      valid
# 11     --> 14, 37      not valid

# res3.2 --> res3.5
# 0      --> 4, 23      not valid
# 2      --> 14, 19     unsure

# res3.2 --> res3.7
# 9      --> 19, 32     probably valid

# res3.7 --> res4.0
# 3      --> 6, 40      valid


antenna.neuron <- FindClusters(antenna.neuron, verbose = T, resolution=4.0)
DefaultAssay(antenna.neuron) <- "SCT"
saveRDS(antenna.neuron, "9_neurons_harmony_rmclusts_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
