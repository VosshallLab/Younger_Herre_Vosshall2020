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


rm(list = ls())   #remove all objects
theme_set(theme_cowplot())

setwd("E:/GDrive Mirror/Lab/Lab Projects/2. scRNAseq/MOC/2022-03-14_pipelinerevisions/1_FINALPIPELINE")

antenna.allcells <- readRDS("7j_harmony_AllCellFinal.rds")  #Uploaded on Zenodo as "SeuratObject1_Antenna_mergedBatches_AllCells.rds"
neuronClusterIDs <- readRDS("8t_neuronClusterIDs.rds")

pipeline <- "harmony_rmclusts_res4"
filename <- '13_allcellsfigures'
dir.create(filename)
outputFolder = paste0(getwd(), '/', filename)
setwd(outputFolder)
getwd()

DefaultAssay(antenna.allcells) <- "SCT"

dim(antenna.allcells@meta.data)[1]


s######################################################################################
# tSNE
######################################################################################
pointsize = 1
# UMAP/tSNE
p1 <- DimPlot(antenna.allcells, reduction = "tsne", label = TRUE, pt.size = pointsize)
filename <- paste0("13a_tsne.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p1)

p1 <- DimPlot(antenna.allcells, reduction = "umap", label = TRUE, pt.size = pointsize)
filename <- paste0("13b_umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p1)

rm(p1, filename)

######################################################################################
# ALL CELLS HEATMAP
######################################################################################

###
### Heatmap for all cells
# Neuron: 'LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204'
# one more neuronal marker: para ortholog LOC5567355
# LOC110678282	Repo	Glia
# LOC5565901	syt1	neuron
# LOC5564305	grh	Epithelia

#FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5567355', label=T)

p01 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC110678282', label=T) + 
  ggtitle('repo (LOC110678282)')
p02 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5565901', label=T) + 
  ggtitle('syt1 (LOC5565901)')
p03 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5564305', label=T) + 
  ggtitle('grh (LOC5564305)')
p04 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5575210', label=T) + 
  ggtitle('nompC (LOC5575210)')

g = arrangeGrob(p01,p02,p03,p04, ncol = 2)
filename <- paste0("13c_clusterannotation_labels.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
rm(g)

p05 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC110678282', label=F) + 
  ggtitle('repo (LOC110678282)')
p06 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5565901', label=F) + 
  ggtitle('syt1 (LOC5565901)')
p07 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5564305', label=F) + 
  ggtitle('grh (LOC5564305)')
p08 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5575210', label=F) + 
  ggtitle('nompC (LOC5575210)')

g = arrangeGrob(p05,p06,p07,p08, ncol = 2)
filename <- paste0("13d_clusterannotation_nolabels.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
rm(g)


#.....................................................................................
# HEATMAP
# ....................................................................................

set.seed(96)

antenna.allcells$annotation <- apply(antenna.allcells@meta.data, 1, function(input){
  if (input['seurat_clusters'] %in% neuronClusterIDs) {return('Neuron')}
  else if (input['seurat_clusters'] %in% c('30', '41')) {return('Glia')}
  else if (input['seurat_clusters'] %in% c('4','19','29')) {return('Epithelia')}
  else {return('Others')}
})

antenna.allcells@meta.data$annotation <- factor(antenna.allcells@meta.data$annotation, 
                                               levels = c('Others', 'Epithelia', 'Glia', 'Neuron'))

mainTissueGeneL <- c(
  'LOC5564305',
  'LOC5570152',
  'LOC5576812',
  'LOC110678282',
  'LOC110676862',
  'LOC5566721',
  'LOC5564848',
  'LOC5570381',
  'LOC5565901',
  'LOC5570204',
  'LOC5567355')

#DefaultAssay(antenna.allcells) <- "RNA"
DefaultAssay(antenna.allcells) <- "SCT"


mainTissueMat <- antenna.allcells[["SCT"]]@data[mainTissueGeneL, ] %>% as.matrix()
# mainTissueMat<- t(scale(t(mainTissueMat)))
# head(mainTissueMat[1:5,1:5])
quantileL <- quantile(mainTissueMat, c(0.1, 0.93))
quantileL <- append(quantileL, mean(quantileL), after = 1)
col_fun = circlize::colorRamp2(quantileL, c('blue', '#CF92A2', 'yellow') ) 

cluster_anno<- antenna.allcells@meta.data$annotation
left_anno <- c(
  replicate(3, 'Epithelia Markers'),
  replicate(3, 'Glia Markers'),
  replicate(5, 'Neuron Markers')
)
# remove: replicate(3, 'Others'),

left_anno <- factor(left_anno,
                    levels = c('Epithelia Markers', 
                               'Glia Markers','Neuron Markers'))
# remove: 'Others', 


p09 <- Heatmap(mainTissueMat, name = "Expression",
               column_split = cluster_anno, # column_split = factor(cluster_anno, levels(antenna.allcells@meta.data$annotation2)),
               column_title_rot = 0,
               height = unit(0.75, "npc"),
               
               col = col_fun,
               # col=viridis(100),
               show_column_names = FALSE, use_raster = FALSE,
               
               cluster_columns = FALSE,
               show_column_dend = FALSE,
               cluster_column_slices = TRUE,
               column_title_gp = gpar(fontsize = 18),
               
               row_names_gp = gpar(fontsize = 12),
               column_gap = unit(0.5, "mm"),
               cluster_rows = FALSE,
               show_row_dend = FALSE,
               cluster_row_slices = TRUE,
               
               row_split = left_anno,
               row_title_gp = gpar(fontsize = 18),
               row_title_rot = 0,
               
               
               top_annotation =
                 HeatmapAnnotation(
                   foo = anno_block(gp = gpar(fill = scales::hue_pal()( length( levels(antenna.allcells@meta.data$annotation)) ) ))
                 ),
               
               # top_annotation = 
               #   HeatmapAnnotation(
               #     foo = anno_block(gp = gpar(fill = brewer.pal(length( levels(antenna.allcells@meta.data$annotation)), "Set3") ) )
               #   ),
               
               # left_annotation =
               #   rowAnnotation(
               #     foo = anno_block(gp = gpar(fill = scales::hue_pal()( length(levels(left_anno)) ) )),
               #     ),
               # 
)

filename <- paste0("13e_heatmap-cellmarkers_5x20_test.pdf")
pdf(file=filename, width=20, height=5)
print(p09)
dev.off()


######################################################################################
# FIGURE S5A
# 1_DECONTX: 
# Before & After Orco Feature Plots
######################################################################################
rm(list = ls())   #remove all objects
setwd()
projectFolder <- getwd()

for (dataset.name in c('antenna1_08192021', 'antenna2_10062021')) {
  #dataset.name <- 
  decontX_outFolder <- paste0(projectFolder,"/", dataset.name,"/decontX_out")
  Before_decontX_Folder <- paste0(projectFolder,"/", dataset.name,"/filtered_feature_bc_matrix")
  
  sratDecontx  <- 
    Read10X(decontX_outFolder) %>%
    CreateSeuratObject(project = dataset.name, min.cells = 12, min.features = 200)
  
  sratDecontx[["percent.mt"]] <- PercentageFeatureSet(sratDecontx, pattern = "^MT-")
  sratDecontx    <- SCTransform(sratDecontx, verbose = F)
  
  sratDecontx    <- RunPCA(sratDecontx, npcs = 50, verbose = F)
  sratDecontx    <- RunUMAP(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- RunTSNE(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindNeighbors(sratDecontx, dims = 1:50, verbose = F)
  sratDecontx    <- FindClusters(sratDecontx, verbose = T)
  
  
  beforeDecontx  <- 
    Read10X(Before_decontX_Folder) %>%
    CreateSeuratObject(project = dataset.name, min.cells = 0, min.features = 200)
  
  beforeDecontx[["percent.mt"]] <- PercentageFeatureSet(beforeDecontx, pattern = "^MT-")
  beforeDecontx    <- SCTransform(beforeDecontx, verbose = F)
  
  beforeDecontx    <- RunPCA(beforeDecontx, npcs = 50, verbose = F)
  beforeDecontx    <- RunUMAP(beforeDecontx, dims = 1:50, verbose = F)
  beforeDecontx    <- RunTSNE(beforeDecontx, dims = 1:50, verbose = F)
  beforeDecontx    <- FindNeighbors(beforeDecontx, dims = 1:50, verbose = F)
  beforeDecontx    <- FindClusters(beforeDecontx, verbose = T)
  
  #....................................................................................
  # PLOT: CORECEPTOR FEATURE PLOTS
  #....................................................................................
  
  p31 = FeaturePlot(beforeDecontx, reduction = 'tsne', features = c('Orco'))
  p32 = FeaturePlot(beforeDecontx, reduction = 'tsne', features = c('Ir25a'))
  p33 = FeaturePlot(beforeDecontx, reduction = 'tsne', features = c('Ir8a'))
  p34 = FeaturePlot(beforeDecontx, reduction = 'tsne', features = c('Ir76b'))
  p35 = FeaturePlot(beforeDecontx, reduction = 'tsne', features = c('Or82'))
  p36 = FeaturePlot(beforeDecontx, reduction = 'tsne', features = c('LOC5575210'))
  
  g = arrangeGrob(p31, p32, p33, p34, p35, p36, ncol = 3)
  filename <- paste0("1c_decontX_coreceptors_",dataset.name,"_beforedecontX.pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 4000,g)
  rm(g, filename)
  
  
  p25 = FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Orco'))
  p26 = FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Ir25a'))
  p27 = FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Ir8a'))
  p28 = FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Ir76b'))
  p29 = FeaturePlot(sratDecontx, reduction = 'tsne', features = c('Or82'))
  p30 = FeaturePlot(sratDecontx, reduction = 'tsne', features = c('LOC5575210'))
  
  g = arrangeGrob(p25, p26, p27, p28, p29, p30, ncol = 3)
  filename <- paste0("1d_decontX_coreceptors_",dataset.name,"_afterdecontX.pdf")
  ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height = 4000,g)
  rm(g, filename)

  
  # Save workspace
  filename <- paste0("1e_decontX_comparison-tsne_",dataset.name, ".rdata")
  print("Saving workspace...")
  save.image(filename)
  
  rm(decontX_outFolder, Before_decontX_Folder, sratDecontx, beforeDecontx, 
     p31, p32, p34, p35, p36, p25, p26, p27, p29, p30)
  
}

#FIG S9-A
load("1e_decontX_comparison-tsne_antenna1_08192021.rdata")
b1_p31 <- p31
b1_p25 <- p25

load("1e_decontX_comparison-tsne_antenna2_10062021.rdata")
b2_p31 <- p31
b2_p25 <- p25

g = arrangeGrob(b1_p31, b2_p31, b1_p25, b2_p25, ncol = 2)
filename <- paste0("14_decontX-before-after.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
rm(g, filename)


rm(list = ls())   #remove all objects




######################################################################################
# FIGURE S5, C-E
######################################################################################
rm(list = ls())   #remove all objects

setwd("E:/GDrive Mirror/Lab/Lab Projects/2. scRNAseq/MOC/2022-03-14_pipelinerevisions/1_FINALPIPELINE")
antenna.allcells <- readRDS("7j_harmony_AllCellFinal.rds")

Idents(object = antenna.allcells) <- "all_cells"

p1 = VlnPlot(antenna.allcells, pt.size=.1, features = 'nFeature_RNA',  split.by = 'batch')+
  theme(legend.position = 'none')
p2 = VlnPlot(antenna.allcells, pt.size=.1, features = 'log_count',  split.by = 'batch')+
  theme(legend.position = 'none')
p3 = VlnPlot(antenna.allcells, pt.size=.1, features = 'percent.mt',  split.by = 'batch')

g = arrangeGrob(p1,p2,p3, ncol = 3)
filename <- paste0("15_violinplots_filtering.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 3000, height = 2000,g)


######################################################################################
# S9-E
######################################################################################


p05 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC110678282', label=F) + 
  ggtitle('repo')
p06 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5565901', label=F) + 
  ggtitle('syt1')
p07 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5564305', label=F) + 
  ggtitle('grh')
p08 <- FeaturePlot(antenna.allcells, reduction='tsne', features='LOC5575210', label=F) + 
  ggtitle('nompC')

g = arrangeGrob(p05,p06,p07,p08, ncol = 4)
filename <- paste0("16_cellmarkers.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 1000,g)
rm(g)


