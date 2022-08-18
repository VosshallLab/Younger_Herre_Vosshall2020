suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(future)
  library(scales)
  library(data.table)
  library(ggplot2)
  library(stats)
  library(gridExtra)
  library(scater)
  library(ggrepel)
  library(ggplot2)
  library(reshape2)
})



######################################################################################
# CELLS BY CHEMORECEPTOR
######################################################################################

rm(list = ls())   #remove all objects
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
pipeline <- "harmony_rmclusts_res4"
dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

#.....................................................................................
# SET UP
# ....................................................................................

outputFeaturePlot <- 'Antenna_SuppData2_FeaturePlots'

dir.create(outputFeaturePlot)

#.....................................................................................
# GENERATE GENE LIST
# ....................................................................................

dat <- antenna.neuron
Idents(object = dat) <- "all_cells"
DefaultAssay(dat) <- dat.assay

ae <- AverageExpression(dat, assays=dat.assay)
ae <- ae[[dat.assay]]
ae <- ae[order(ae[,1],decreasing=TRUE),]      
ae <- as.data.frame(ae)                       
ae <- tibble::rownames_to_column(ae, "Gene")  
ae <- ae[!grepl("LOC*", ae$Gene),]            
ae <- ae[!grepl("MT*", ae$Gene),]
ae <- ae[!grepl("Orco", ae$Gene),]
ae <- ae[!grepl("Ir25a", ae$Gene),]
ae <- ae[!grepl("Ir76b", ae$Gene),]
ae <- ae[!grepl("Ir8a", ae$Gene),]

GeneList <- ae
GeneList$ae <- seq(1:nrow(GeneList))

names(GeneList)[names(GeneList) == 'ae'] <- "Number"
GeneList$Number[1:100] <- sprintf('%0.3d', 1:100)

coreceptors <- c('LOC5575210','Orco','Ir25a','Ir76b','Ir8a')
receptorList <- append(coreceptors, GeneList$Gene) 

rm(ae, dat, GeneList)

#.....................................................................................
# GENE OF INTEREST
# ....................................................................................

for (geneName in receptorList) {

  #geneName <- "Or47"
  geneNumber <- sprintf('%0.3d', grep(paste0(geneName,"\\b"), receptorList))
  expression.cutoff <- 1
  dat1 <- antenna.neuron
  totalcells1 <- dim(antenna.neuron@meta.data)[1]
  
  ################################################################## Heatmap 1, no recluster
  
  p1 <- FeaturePlot(dat1,  reduction = 'tsne', features = c(geneName) , pt.size = .1, slot = "data")
  p2 <- FeaturePlot(dat1,  reduction = 'umap', features = c(geneName) , pt.size = .1, slot = "data")
  
  g = arrangeGrob(p1,p2, ncol = 2)
  filename <- paste0("19a_",geneNumber, "_", geneName,"_",pipeline,"_featurePlot.pdf")
  ggsave(filename, path=outputFeaturePlot, limitsize = FALSE, units = "px", width = 4000, height = 2000, g)
  
  #DimPlot(dat1.cl, reduction = "umap", group.by = 'dataset') + 
  #  DimPlot(dat1.cl, reduction = "tsne", group.by = 'dataset')
  
}
