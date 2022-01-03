suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
})

theme_set(theme_cowplot())
projectFolder='/Users/tzuchiao/ProjectAnalysis/Leslie_Vosshall'
setwd(projectFolder)



#####################################

### 6. Neuron subcluster

normalizationMethod='LogNormalize'
# normalizationMethod='SCTransform' # not good for visualization
print(c('normalizationMethod: ', normalizationMethod))

receptorFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/5.checkReceptors')
neuronFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/6.neuronClusters')
dir.create(neuronFolderPath, showWarnings = FALSE)
# picFolderPath=paste0(neuronFolderPath, '/pic')
# dir.create(picFolderPath, showWarnings = FALSE)


#####################################

load(file = file.path('Analysis/integratedData_decontX', 'antenna_5.antenna.neuron.RData') )
picFolderPath=paste0(neuronFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)

print(neuronFolderPath)
print(picFolderPath)


##############################################
### batch correction of the neuronal clusters
neuron.list <- SplitObject(antenna.neuron, split.by = "batch")
neuron.list

# normalize and identify variable features for each dataset independently
neuron.list <- lapply(X = neuron.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
# run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = neuron.list)
anchors <- FindIntegrationAnchors(object.list = neuron.list, anchor.features = features)
# this command creates an 'integrated' data assay
neuron.batchCorrect <- IntegrateData(anchorset = anchors)



### Perform an integrated analysis
DefaultAssay(neuron.batchCorrect) <- "integrated"

neuron.batchCorrect <- ScaleData(neuron.batchCorrect, vars.to.regress = c('nCount_RNA'))

# PCA, UMAP, tsne, clustering
neuron.batchCorrect    <- RunPCA(neuron.batchCorrect, npcs = 50, verbose = F)
neuron.batchCorrect    <- RunUMAP(neuron.batchCorrect, dims = 1:50, verbose = F)
neuron.batchCorrect    <- RunTSNE(neuron.batchCorrect, dims = 1:50, verbose = F)
neuron.batchCorrect    <- FindNeighbors(neuron.batchCorrect, dims = 1:50, verbose = F)
neuron.batchCorrect    <- FindClusters(neuron.batchCorrect, resolution = c(0.5, 0.8, 1, 1.5, 2, 3)[6])
# neuron.batchCorrect    <- FindClusters(neuron.batchCorrect, resolution = c(0.5, 0.8, 1, 1.2, 2, 3.3)[6])

pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_', 'batchCorrection_neuronOnly_batch.pdf')),
    width=13, height=6)
print(( DimPlot(neuron.batchCorrect, reduction = "umap", label = TRUE) + theme(legend.position = "none")) + 
          DimPlot(neuron.batchCorrect, reduction = "umap", group.by = 'batch'))
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_', 'batchCorrection_neuronOnly.pdf')),
    width=5, height=4.5)
print( DimPlot(neuron.batchCorrect, reduction = "umap", label = TRUE) + theme(legend.position = "none") )
dev.off()


pdf(file = file.path(picFolderPath, 
                     paste0('tSNE_', 'batchCorrection_neuronOnly_batch.pdf')),
    width=13, height=6)
print(( DimPlot(neuron.batchCorrect, reduction = "tsne", label = TRUE) + theme(legend.position = "none")) + 
       DimPlot(neuron.batchCorrect, reduction = "tsne", group.by = 'batch'))
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('tSNE_', 'batchCorrection_neuronOnly.pdf')),
    width=5, height=4.5)
print( DimPlot(neuron.batchCorrect, reduction = "tsne", label = TRUE) + theme(legend.position = "none") )
dev.off()


rm(anchors, antenna.neuron, neuron.list)

save.image(file = file.path('Analysis/integratedData_decontX', 'antenna_6.neuralClusters.RData') )






