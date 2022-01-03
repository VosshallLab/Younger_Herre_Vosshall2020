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

### 3. check clustering, , define neuron cells

# for (normalizationMethod in c('LogNormalize', 'SCTransform')) {

dataFolder='10282021_MaxPalp'
normalizationMethod='LogNormalize'
# normalizationMethod='SCTransform' # not good for visualization
print(c('normalizationMethod: ', normalizationMethod))

dfFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/2.doubletfinder')
receptorFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/3.checkReceptors')
dir.create(receptorFolderPath, showWarnings = FALSE)
picFolderPath=paste0(receptorFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)

print(receptorFolderPath)
print(picFolderPath)

sratDF <- readRDS( file.path('Analysis/10282021_MaxPalp', 'maxpalp_2.DoubletFinder.rds') )


### FindNeighbors & FindClusters

sratDF <- NormalizeData(sratDF)
sratDF <- FindVariableFeatures(sratDF, selection.method = "vst", nfeatures = 2000)
sratDF <- ScaleData(sratDF, vars.to.regress = c('nCount_RNA'))

sratDF <- RunPCA(sratDF, npcs = 50, verbose = F)
DimPlot(sratDF, reduction = "pca")
ElbowPlot(sratDF, ndims = 50)

# VariableFeatures
top10 <- head(VariableFeatures(sratDF), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sratDF)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


sratDF <- FindNeighbors(sratDF, reduction = "pca", dims = 1:50)
sratDF <- FindClusters(sratDF, resolution = c(0.5, 0.8, 1, 2.5)[4])

### UMAP
sratDF <- RunUMAP(sratDF, reduction = "pca", dims = 1:50, metric = c("cosine", 'correlation', 'euclidean')[1]  )
sratDF <- RunTSNE(sratDF, reduction = "pca", dims = 1:50)
# Visualization
# Idents(sratDF) <- sratDF@meta.data$seurat
# DimPlot(sratDF, reduction = "umap", label = TRUE)
# DimPlot(sratDF, reduction = "tsne", label = TRUE)

pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_tSNE', '.pdf')),
    width=10, height=5)
print( ( DimPlot(sratDF, reduction = "umap", label = TRUE) + theme(legend.position = "none")) + 
         ( DimPlot(sratDF, reduction = "tsne", label = TRUE) + theme(legend.position = "none")) )
dev.off()



### Find markers
library(future)
plan("multiprocess", workers = 4)

DefaultAssay(sratDF) <- "RNA"
merge.markers <- FindAllMarkers(sratDF, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(merge.markers, file = 'Analysis/10282021_MaxPalp/LogNormalize/3.checkReceptors/allMarkers.csv')

### coreceptor
pdf(file = file.path(picFolderPath, 
                     paste0('coreceptor_FeaturePlot.pdf')),
    width=15, height=10)
print(FeaturePlot(sratDF, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Gr3'), ncol = 3, order = TRUE, pt.size = 1))
dev.off()

### Neuronal markers
pdf(file = file.path(picFolderPath, 
                     paste0('neuronMarkers_FeaturePlot.pdf')),
    width=10, height=10)
print(FeaturePlot(sratDF, reduction = 'tsne', features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'), ncol = 2, order = TRUE, pt.size = 1))
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('neuronMarkers_DotPlot.pdf')),
    width=7, height=7)
print(DotPlot(sratDF, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'))+ scale_y_discrete(limits = rev))
dev.off()


### Plot major cell types
#############################################################################################################################################
#############################################################################################################################################
# Other markers ==> not finished !!!
# LOC110678282	Repo	Glia
# LOC5565901	syt1	neuron
# LOC5564305	gah	Epithelia
# LOC5566990	Mhc1 Muscle
# LOC5575210  nompC ?

pdf(file = file.path(picFolderPath, 
                     paste0('OtherMarkers_FeaturePlot.pdf')),
    width=10, height=10)
print(FeaturePlot(sratDF, reduction = 'tsne',features = c('LOC5565901', 'LOC5564305', 'LOC110678282', 'LOC5566990'), ncol = 2,  order = TRUE, pt.size = 1))
dev.off()

#############################################################################################################################################
#############################################################################################################################################


####
### Subcluster neuron cells


### Plot neuronal markers

neuronM.plot <- DotPlot(sratDF, features = c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204')) + scale_y_discrete(limits = rev)
# neuronM.plot
neuronM.percentage  = neuronM.plot$data
# neuronM.percentage$geneName <- row.names(neuronM.percentage)
# neuronM.percentage <- neuronM.percentage %>% filter(features.plot != 'LOC5578234')

pdf(file = file.path(picFolderPath, 
                     paste0('neuronMarkers_scatter_avg.exp_pct.exp.pdf')),
    width=5, height=5)
print(neuronM.percentage %>%
        ggplot(aes(avg.exp, pct.exp, color=id)) + geom_point())
dev.off()

# colnames(neuronM.percentage)

neuronM.percentage.wider <- neuronM.percentage %>%
  pivot_wider(names_from = features.plot, values_from = c(avg.exp, pct.exp, avg.exp.scaled))

# neuronM.percentage.wider$ratio45N <- apply(neuronM.percentage.wider, 1, function(input){
#   # print( as.numeric(input[6:9])  )
#   ratio45counts=0
#   for (ratio in as.numeric(input[6:9])) {
#     if (ratio > 45) {ratio45counts= ratio45counts + 1}
#   }
#   return( ratio45counts )
# })
neuronM.percentage.wider$ratio50N <- apply(neuronM.percentage.wider, 1, function(input){
  # print( as.numeric(input[6:9])  )
  ratio50counts=0
  for (ratio in as.numeric(input[6:9])) {
    if (ratio > 50) {ratio50counts= ratio50counts + 1}
  }
  return( ratio50counts )
})

#
neuronClusterIDs <- neuronM.percentage.wider %>% filter(ratio50N >= 3) %>% pull(id) %>% as.character()
neuronClusterIDs
length(neuronClusterIDs)
srat.neuron <- subset(sratDF, subset = seurat_clusters %in% neuronClusterIDs)
srat.neuron$originalClusters <- srat.neuron$seurat_clusters
srat.neuron
# antenna.nonNeuron <- subset(sratDF, subset = !(seurat_clusters %in% neuronClusterIDs))
# antenna.nonNeuron


# saveRDS(srat.neuron, file = 
#           file.path('Analysis/10282021_MaxPalp', 'maxpalp_3.neuron.rds') )

save.image(file = file.path('Analysis/10282021_MaxPalp', 'maxpalp_3.neuron.RData') )






