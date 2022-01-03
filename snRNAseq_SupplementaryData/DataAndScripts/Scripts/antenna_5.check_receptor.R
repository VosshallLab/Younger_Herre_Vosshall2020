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

### 5. check clustering, , define neuron cells

# for (normalizationMethod in c('LogNormalize', 'SCTransform')) {

normalizationMethod='LogNormalize'
# normalizationMethod='SCTransform' # not good for visualization
print(c('normalizationMethod: ', normalizationMethod))

anchorFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/4.anchor')
receptorFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/5.checkReceptors')
dir.create(receptorFolderPath, showWarnings = FALSE)
picFolderPath=paste0(receptorFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)
picFolderPath=paste0(receptorFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)

print(anchorFolderPath)
print(picFolderPath)

seurat.combined <- readRDS(file = file.path('Analysis/integratedData_decontX', 'antenna_4.seurat.combined.rds') )

#
seurat.combined@meta.data$batch <- apply(seurat.combined@meta.data, 1, FUN = function(input) {
  dataDate = as.character(input['dataset'])
  print(dataDate)
  if (dataDate == '0819') {return('batch1')}
  else {return('batch2')}
})
# head(seurat.combined)
# tail(seurat.combined)



DefaultAssay(seurat.combined) <- "integrated"

### FindNeighbors & FindClusters
# seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:50)
seurat.combined <- FindClusters(seurat.combined, resolution = c(0.5, 0.8, 1, 2)[4])

### UMAP
seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:50, metric = c("cosine", 'correlation', 'euclidean')[1]  )
seurat.combined <- RunTSNE(seurat.combined, reduction = "pca", dims = 1:50)
# Visualization
# Idents(seurat.combined) <- seurat.combined@meta.data$seurat
# DimPlot(seurat.combined, reduction = "umap", label = TRUE)
# DimPlot(seurat.combined, reduction = "tsne", label = TRUE)

pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_tSNE_', 'batchCorrection.pdf')),
    width=15, height=6)
print(DimPlot(seurat.combined, reduction = "umap", label = TRUE) + DimPlot(seurat.combined, reduction = "tsne", label = TRUE))
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_tSNE_', 'batchCorrection_dataset.pdf')),
    width=13, height=6)
print(DimPlot(seurat.combined, reduction = "umap",  group.by = 'dataset') + 
        DimPlot(seurat.combined, reduction = "tsne", group.by = 'dataset'))
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_tSNE_', 'batchCorrection_dataset_split.pdf')),
    width=13, height=12)
print(DimPlot(seurat.combined, reduction = "umap", split.by = 'dataset', label = TRUE) / 
        DimPlot(seurat.combined, reduction = "tsne", split.by = 'dataset', label = TRUE))
dev.off()


# DimPlot(seurat.combined, reduction = "umap", label = TRUE)
# DimPlot(seurat.combined, reduction = "tsne", label = TRUE)


### Find markers
library(future)
plan("multiprocess", workers = 4)

DefaultAssay(seurat.combined) <- "RNA"
merge.markers <- FindAllMarkers(seurat.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(merge.markers, file = 'Analysis/integratedData_decontX/LogNormalize/5.checkReceptors/allMarkers.csv')

### coreceptor
pdf(file = file.path(picFolderPath, 
                     paste0('coreceptor_FeaturePlot_batchCorrect.pdf')),
    width=15, height=15)
print(FeaturePlot(seurat.combined, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b'), ncol = 2))
dev.off()

### Neuronal markers
pdf(file = file.path(picFolderPath, 
                     paste0('neuronMarkers_FeaturePlot_batchCorrect.pdf')),
    width=15, height=15)
print(FeaturePlot(seurat.combined, reduction = 'tsne', features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'), ncol = 2))
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('neuronMarkers_DotPlot_batchCorrect.pdf')),
    width=9, height=11)
print(DotPlot(seurat.combined, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'))+ scale_y_discrete(limits = rev))
dev.off()


### Plot major cell types
#############################################################################################################################################
#############################################################################################################################################
# LOC110678282	Repo	Glia
# LOC5565901	syt1	neuron
# LOC5564305	gah	Epithelia
# LOC5566990	Mhc1 Muscle
# LOC5575210  nompC 
FeaturePlot(seurat.combined, reduction = 'tsne',features = c('LOC5565901', 'LOC5564305', 'LOC110678282', 'LOC5566990'), ncol = 2)
#############################################################################################################################################
#############################################################################################################################################


####
### Subcluster neuron cells


### Plot neuronal markers

neuronM.plot <- DotPlot(seurat.combined, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204')) + scale_y_discrete(limits = rev)
# neuronM.plot
neuronM.percentage  = neuronM.plot$data
# neuronM.percentage$geneName <- row.names(neuronM.percentage)
# neuronM.percentage <- neuronM.percentage %>% filter(features.plot != 'LOC5578234')


# colnames(neuronM.percentage)

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
neuronClusterIDs
length(neuronClusterIDs)
antenna.neuron <- subset(seurat.combined, subset = seurat_clusters %in% neuronClusterIDs)
antenna.neuron
# antenna.nonNeuron <- subset(seurat.combined, subset = !(seurat_clusters %in% neuronClusterIDs))
# antenna.nonNeuron


# saveRDS(antenna.neuron, file = 
#           file.path('Analysis/integratedData_decontX', 'antenna_5.antenna.neuron.rds') )

save.image(file = file.path('Analysis/integratedData_decontX', 'antenna_5.antenna.neuron.RData') )






