
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



set.seed(96)
rm(list = ls())   #remove all objects

setwd()
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
pipeline <- "harmony_rmclusts_res4"

dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolder <- "Antenna_Figures3_NeuronExampleHeatmaps_Output"
dir.create(outputFolder)

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
#saveRDS(receptorList, "receptorList.rds")



######################################################################################
# Or82
######################################################################################

#.....................................................................................
# SET UP
# ....................................................................................

DefaultAssay(antenna.neuron) <- dat.assay

geneName <- "Or82"
geneNumber <- grep(paste0(geneName,"\\b"), receptorList)
print(geneName)

dat1 <- antenna.neuron
expression.cutoff <- 1
totalcells1 <- dim(antenna.neuron@meta.data)[1]

#.....................................................................................
# GENE OF INTEREST
# ....................................................................................

for (geneName in "Or82") {
  geneNumber <- sprintf('%0.3d', grep(paste0(geneName,"\\b"), receptorList))
  
  expression.cutoff <- 1
  
  dat1 <- antenna.neuron
  totalcells1 <- dim(antenna.neuron@meta.data)[1]
  
  ############################################################################################################### Heatmap 1, no recluster
  
  expr <- FetchData(object = dat1[[dat.assay]], vars = geneName)
  dat1 <- dat1[, which(x = expr > expression.cutoff)]
  subsetcellnumber <- dim(dat1@meta.data)[1]
  print(paste0(toString(subsetcellnumber), " cells expressing ",geneName))
  
  heatmap.title1 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Original Clusters",
                          "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  numberClusters <- length(levels(Idents(dat1)))
  clusterList <- levels(Idents(dat1))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  allgenes <- unique(allgenes)
  allgenes.og <- allgenes
  
  rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)
  
  
  ############################################################################################################### Heatmap 2, reclustered
  
  dat1.cl <- dat1
  
  if (grepl("harm*", pipeline)) {
    #  dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #  print("Variable features imported from 'integrated' assay.")
    reduction <- "harmony"
    dat.assay<-"SCT"
    
  }
  
  if (grepl("anch*", pipeline)) {
    #dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #print("Variable features imported from 'integrated' assay.")
    
    dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    reduction <- "pca"
    
    DefaultAssay(dat1.cl) <- "integrated"
    
    dat1.cl    <- RunPCA(dat1.cl, npcs = 50, verbose = F)
  }
  
  recluster.code <- "("
  
  try({
    dat1.cl <- RunUMAP(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] UMAP successful"))
    recluster.code <- append(recluster.code, "U")
  })
  
  try({
    dat1.cl <- RunTSNE(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] tSNE successful"))
    recluster.code <- append(recluster.code, "S")
  })
  
  
  try({
    dat1.cl <- FindNeighbors(dat1.cl, assay=dat.assay, 
                             dims = 1:50, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindNeighbors successful"))
    recluster.code <- append(recluster.code, "F")
  })
  
  try({
    dat1.cl <- FindClusters(dat1.cl, resolution = 1, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindClusters successful"))
    recluster.code <- append(recluster.code, "C")
  })
  
  
  heatmap.title2 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Reclustered ",paste(unlist(recluster.code),collapse=""),")",
                          "\n",toString(dim(dat1.cl@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  
  numberClusters <- length(levels(Idents(dat1.cl)))
  clusterList <- levels(Idents(dat1.cl))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1.cl,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  #allgenes <- append(allgenes, allgenes.og)
  allgenes <- unique(allgenes)
  
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, features = allgenes.og, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)
  
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1.cl, features = allgenes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("14a_",geneNumber, "_", geneName,"_",pipeline,"_rawHeatmaps.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
  
  rm(g)
  
  rm(u, cluster.subset, i, numberClusters, filename, cl.ae, allgenes, clusterList)
  
  
  #..............................................................................................UMI violin
  
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  p4 <- VlnPlot(dat1.cl, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  
  g = arrangeGrob(p3,p4, ncol = 2)
  filename <- paste0("14b_",geneNumber, "_", geneName,"_",pipeline,"_UMIviolin.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000, g)
  
  p5 <- DimPlot(dat1, reduction = "tsne", group.by = 'dataset')
  p6 <- DimPlot(dat1, reduction = "umap", group.by = 'dataset')
  p7 <- DimPlot(dat1.cl, reduction = "tsne", group.by = 'dataset')
  p8 <- DimPlot(dat1.cl, reduction = "umap", group.by = 'dataset')
  g = arrangeGrob(p5,p6,p7,p8, ncol = 2)
  filename <- paste0("14c_",geneNumber, "_", geneName,"_",pipeline,"_batchcheck.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000, g)
  
  
}


#.....................................................................................
# FIGURE PRESENTATION
# ....................................................................................

dat2 <- dat1

levels(Idents(dat2))

my_levels <- c("37","19","32","38","1",  "3" , "5",  "13" ,"14" ,
               "17" ,"21" ,"23", "25" ,"30", 
               "34","40")

dat2@active.ident <- factor(x = dat2@active.ident, levels = my_levels)

levels(Idents(dat2))

heatmap.genes <- c("LOC5575210","Orco","Ir25a","Ir76b","Ir8a","Or82",
                   "Or6","Or71",
                   "Or4","Or47","Or3", "Or50","Gr15",
                   "Ir41l","Ir41m", 
                   "Or81"     , "Or80")


p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat2, features = heatmap.genes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)

p1
filename <- paste0("14d_",geneNumber, "_", geneName,"_",pipeline,"_figurePresentation.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000,p1)

#.....................................................................................
# Checking cluster 38, that Or82+/Ir41l+ appear in both batches
# ....................................................................................


p4 <- DimPlot(antenna.neuron, reduction = "tsne", group.by = 'dataset')+
  DimPlot(antenna.neuron, reduction = "tsne", label=T)
p4
filename <- paste0("18_antennaneurons_batchcheck.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 2000,p4)


dat_Vln <- subset(dat2,idents = c(38))
p3 <- VlnPlot(dat_Vln, assay=dat.assay, slot="data", features=c("Orco","Ir25a", "Or82","Ir41l"), pt.size = 1, 
              split.by = "dataset", ncol=2) + 
  scale_y_continuous(breaks=seq(0,50,1))
filename <- paste0("18a_cl38_batchcheck.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000,p3)

p3 <- VlnPlot(antenna.neuron, assay=dat.assay, slot="data", features=c("Orco","Ir25a","Ir76b","Ir8a"), pt.size = 1, 
              split.by = "dataset", ncol=1) + 
  scale_y_continuous(breaks=seq(0,50,1))
p3
filename <- paste0("18a_allclusters_batchcheck.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 8000,p3)




###########################################################################################
# CLUSTER 0
###########################################################################################

setwd()
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
pipeline <- "harmony_rmclusts_res4"

dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolder <- "Antenna_Figures3_NeuronExampleHeatmaps_Output"


dat <- antenna.neuron
numberClusters <- length(levels(Idents(dat)))
clusterList <- levels(Idents(dat))

for (u in c("0")) {
  
  dat1 <- subset(dat,idents = u)
  print(paste0(dim(dat1@meta.data)[1]," cells in cluster ",toString(u)))
  
  if (grepl("harm*", pipeline)) {
    #  dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #  print("Variable features imported from 'integrated' assay.")
    #print("harmony")
    reduction <- "harmony"
    dat.assay<-"SCT"
    
  }
  
  if (grepl("anch*", pipeline)) {
    #dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #print("Variable features imported from 'integrated' assay.")
    print("anchor")
    dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    reduction <- "pca"
    
    DefaultAssay(dat1.cl) <- "integrated"
    
    dat1    <- RunPCA(dat1.cl, npcs = 50, verbose = F)
  }
  
  recluster.code <- "("
  
  try({
    dat1 <- RunUMAP(dat1, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0("Cluster ",u," UMAP successful"))
    recluster.code <- append(recluster.code, "U")
  })
  
  try({
    dat1 <- RunTSNE(dat1, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0("Cluster ",u," tSNE successful"))
    recluster.code <- append(recluster.code, "S")
  })
  
  try({
    dat1 <- FindNeighbors(dat1, assay=dat.assay, 
                          dims = 1:50, verbose = F)
    print(paste0("Cluster ",u," FindNeighbors successful"))
    recluster.code <- append(recluster.code, "F")
  })
  
  try({
    dat1 <- FindClusters(dat1, resolution = 4, verbose = F)
    print(paste0("Cluster ",u," FindClusters successful"))
    recluster.code <- append(recluster.code, "C")
  })
  
  
  heatmap.genes <- c("LOC5575210","Orco", "Ir25a", "Ir76b", "Ir8a")
  
  dat.nocl <- dat1
  Idents(object = dat.nocl) <- "all_cells"
  DefaultAssay(dat.nocl) <- dat.assay
  
  
  if (dim(dat.nocl@meta.data)[1] >= 1) {
    cl.ae <- as.matrix(dat1@assays[[dat.assay]]@data)
  }
  
  if (dim(dat.nocl@meta.data)[1] > 1) {
    suppressMessages({
      cl.ae <- AverageExpression(dat.nocl, assays=dat.assay)})                          #get avg expression of all genes in cluster      
    cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
  }
  
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
  
  average.genes <- append(heatmap.genes, head(cl.ae$Gene, 20))
  
  rm(cl.ae, dat.nocl)  
  
  
  
  numberClusters1 <- length(levels(Idents(dat1)))
  clusterList1 <- levels(Idents(dat1))
  heatmap.genes <- c("LOC5575210","Orco", "Ir25a", "Ir76b", "Ir8a")
  
  for (b in 1:numberClusters1) {
    a <- clusterList1[b]
    
    cluster.subset <- subset(dat1,idents = a)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    
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
    #head(cl.ae$Gene, 20)
    
    heatmap.genes <- append(heatmap.genes, head(cl.ae$Gene, 5))
    
  }
  heatmap.genes <- unique(heatmap.genes)
  
  heatmap.title1 = paste("Cluster ", toString(u),", Batch Correction: ",pipeline,
                         "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(dim(antenna.neuron@meta.data)[1])," total neurons, reclustered ",paste(unlist(recluster.code),collapse=""),")",
                         "\n","Gene List: Total Average Cluster Expression, Top 20",sep="")
  
  heatmap.title2 = paste("Cluster ", toString(u),", Batch Correction: ",pipeline,
                         "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(dim(antenna.neuron@meta.data)[1])," total neurons, reclustered ",paste(unlist(recluster.code),collapse=""),")",
                         "\n","Gene List: Sub Cluster Average Expression, Top 5",sep="")
  
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                  features = average.genes, label=TRUE, group.bar.height =0.01, 
                  draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
    ggtitle(heatmap.title1)
  
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                  features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                  draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
    ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("16a_cluster",toString(u),"_rawHeatmaps.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
  
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = average.genes[1:12], pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  #p3
  filename <- paste0("16c_cluster",toString(u),"_violinSCTcounts.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000, p3)
  
  p3 <- VlnPlot(dat1, assay="RNA", slot="counts", features = average.genes[1:12], pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  #p3
  filename <- paste0("16d_cluster",toString(u),"_violinRNAcounts.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000, p3)
  
  
  p5 <- DimPlot(dat1, reduction = "tsne", group.by = 'dataset', label=T)+ggtitle(paste0("cluster ",toString(u)," by batch"))
  p6 <- DimPlot(dat1, reduction = "umap", group.by = 'dataset', label=T)+ggtitle(paste0("cluster ",toString(u)," by batch"))
  #p5+p6
  g = arrangeGrob(p5,p6, ncol = 2)
  filename <- paste0("16e_cluster",toString(u),"_batchcheck.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 2000, g)
  
}
#............................................................................................................................
# FIGURE PRESENTATION

#mini.heatmap.genes <- average.genes[1:10]

mini.heatmap.genes <-  c( "LOC5575210", "Orco"   ,    "Ir25a"   ,   "Ir76b",     
                          "Ir8a"  ,  "Ir41o" ,"Ir41p" ,  "Ir41k","Ir41j",   "Ir41a")

my_levels <- c( "1",  "0"  ,"4"  ,"2",  "3" , "5",  "6"  ,"7",  "8"  ,"9",  "10" ,"11",
                "12" ,"13", "14", "15" ,"16", "17" ,"18", "19" ,"20", "21", "22", "23",
                "24")

dat1@active.ident <- factor(x = dat1@active.ident, levels = my_levels)

levels(Idents(dat1))


heatmap.title3 = paste("Cluster ", toString(u),", Batch Correction: ",pipeline,
                       "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(dim(antenna.neuron@meta.data)[1])," total neurons, reclustered ",paste(unlist(recluster.code),collapse=""),")",
                       "\n","Gene List: Sub Cluster Average Expression, Top 6",sep="")



p8 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                features = mini.heatmap.genes, label=T, group.bar.height = 0.01, 
                draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
  ggtitle(heatmap.title3)
p8
filename <- paste0("16g_cluster",toString(u),"_figurePresentation.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 6000, height = 4000,p8)





###########################################################################################
# CLUSTER 3
###########################################################################################

setwd()
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
pipeline <- "harmony_rmclusts_res4"

dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolder <- "Antenna_Figures3_NeuronExampleHeatmaps_Output"

dat <- antenna.neuron
numberClusters <- length(levels(Idents(dat)))
clusterList <- levels(Idents(dat))

for (u in c("3")) {
  
  dat1 <- subset(dat,idents = u)
  print(paste0(dim(dat1@meta.data)[1]," cells in cluster ",toString(u)))
  
  if (grepl("harm*", pipeline)) {
    #  dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #  print("Variable features imported from 'integrated' assay.")
    #print("harmony")
    reduction <- "harmony"
    dat.assay<-"SCT"
    
  }
  
  if (grepl("anch*", pipeline)) {
    #dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #print("Variable features imported from 'integrated' assay.")
    print("anchor")
    dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    reduction <- "pca"
    
    DefaultAssay(dat1.cl) <- "integrated"
    
    dat1    <- RunPCA(dat1.cl, npcs = 50, verbose = F)
  }
  
  recluster.code <- "("
  
  try({
    dat1 <- RunUMAP(dat1, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0("Cluster ",u," UMAP successful"))
    recluster.code <- append(recluster.code, "U")
  })
  
  try({
    dat1 <- RunTSNE(dat1, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0("Cluster ",u," tSNE successful"))
    recluster.code <- append(recluster.code, "S")
  })
  
  try({
    dat1 <- FindNeighbors(dat1, assay=dat.assay, 
                          dims = 1:50, verbose = F)
    print(paste0("Cluster ",u," FindNeighbors successful"))
    recluster.code <- append(recluster.code, "F")
  })
  
  try({
    dat1 <- FindClusters(dat1, resolution = 4, verbose = F)
    print(paste0("Cluster ",u," FindClusters successful"))
    recluster.code <- append(recluster.code, "C")
  })
  
  
  heatmap.genes <- c("LOC5575210","Orco", "Ir25a", "Ir76b", "Ir8a")
  
  dat.nocl <- dat1
  Idents(object = dat.nocl) <- "all_cells"
  DefaultAssay(dat.nocl) <- dat.assay
  
  
  if (dim(dat.nocl@meta.data)[1] >= 1) {
    cl.ae <- as.matrix(dat1@assays[[dat.assay]]@data)
  }
  
  if (dim(dat.nocl@meta.data)[1] > 1) {
    suppressMessages({
      cl.ae <- AverageExpression(dat.nocl, assays=dat.assay)})                          #get avg expression of all genes in cluster      
    cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
  }
  
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
  
  average.genes <- append(heatmap.genes, head(cl.ae$Gene, 20))
  
  rm(cl.ae, dat.nocl)  
  
  
  
  numberClusters1 <- length(levels(Idents(dat1)))
  clusterList1 <- levels(Idents(dat1))
  heatmap.genes <- c("LOC5575210","Orco", "Ir25a", "Ir76b", "Ir8a")
  
  for (b in 1:numberClusters1) {
    a <- clusterList1[b]
    
    cluster.subset <- subset(dat1,idents = a)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    
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
    #head(cl.ae$Gene, 20)
    
    heatmap.genes <- append(heatmap.genes, head(cl.ae$Gene, 5))
    
  }
  heatmap.genes <- unique(heatmap.genes)
  
  heatmap.title1 = paste("Cluster ", toString(u),", Batch Correction: ",pipeline,
                         "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(dim(antenna.neuron@meta.data)[1])," total neurons, reclustered ",paste(unlist(recluster.code),collapse=""),")",
                         "\n","Gene List: Total Average Cluster Expression, Top 20",sep="")
  
  heatmap.title2 = paste("Cluster ", toString(u),", Batch Correction: ",pipeline,
                         "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(dim(antenna.neuron@meta.data)[1])," total neurons, reclustered ",paste(unlist(recluster.code),collapse=""),")",
                         "\n","Gene List: Sub Cluster Average Expression, Top 5",sep="")
  
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                  features = average.genes, label=TRUE, group.bar.height =0.01, 
                  draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
    ggtitle(heatmap.title1)
  
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                  features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                  draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
    ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("16a_cluster",toString(u),"_rawHeatmaps.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
  
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = average.genes[1:12], pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  #p3
  filename <- paste0("16c_cluster",toString(u),"_violinSCTcounts.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000, p3)
  
  p3 <- VlnPlot(dat1, assay="RNA", slot="counts", features = average.genes[1:12], pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  #p3
  filename <- paste0("16d_cluster",toString(u),"_violinRNAcounts.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000, p3)
  
  
  p5 <- DimPlot(dat1, reduction = "tsne", group.by = 'dataset', label=T)+ggtitle(paste0("cluster ",toString(u)," by batch"))
  p6 <- DimPlot(dat1, reduction = "umap", group.by = 'dataset', label=T)+ggtitle(paste0("cluster ",toString(u)," by batch"))
  #p5+p6
  g = arrangeGrob(p5,p6, ncol = 2)
  filename <- paste0("16e_cluster",toString(u),"_batchcheck.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 2000, g)
  
}

#............................................................................................................................
# FIGURE PRESENTATION

#mini.heatmap.genes <- heatmap.genes

mini.heatmap.genes <-  c("LOC5575210", "Orco"  ,"Ir25a" ,"Ir76b",
                         "Ir8a"      , "Or54"  ,"Or30" ,
                         "Or57"      , "Or51"  ,  "Or47" , "Or55","Or27"  , "Or32" ,
                         "Or84" ,   "Or79" , "Or26", 
                         "Or78",
                         "Or38"  ,
                         "Or29",
                         "Or14"     ,  "Or15" , "Ir75g", "Or33", "Ir7j" ,
                         "Or25"     )

my_levels <- c( "0" , "5"  ,"6", "1" , "2" , "3"  ,"4",   "7"  ,"8",  "9" , "10" ,"11",
                "12", "13", "14" ,"15", "16" ,"17", "18" ,"19", "20" ,"21", "22", "23")

dat1@active.ident <- factor(x = dat1@active.ident, levels = my_levels)

levels(Idents(dat1))


heatmap.title3 = paste("Cluster ", toString(u),", Batch Correction: ",pipeline,
                       "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(dim(antenna.neuron@meta.data)[1])," total neurons, reclustered ",paste(unlist(recluster.code),collapse=""),")",
                       sep="")



p8 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                features = mini.heatmap.genes, label=T, group.bar.height = 0.01, 
                draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
  ggtitle(heatmap.title3)
p8
filename <- paste0("16g_cluster",toString(u),"_figurePresentation.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 6000, height = 4000,p8)






######################################################################################
# Or36
######################################################################################

#.....................................................................................
# SET UP
# ....................................................................................
setwd()
rm(list = ls())
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
receptorList <- readRDS("Antenna_receptorList.rds")
pipeline <- "harmony_rmclusts_res4"

dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolder <- "Antenna_Figures3_NeuronExampleHeatmaps_Output"


geneName <- "Or36"
geneNumber <- grep(paste0(geneName,"\\b"), receptorList)
print(geneName)

dat1 <- antenna.neuron
expression.cutoff <- 1
totalcells1 <- dim(antenna.neuron@meta.data)[1]

#.....................................................................................
# GENE OF INTEREST
# ....................................................................................

for (geneName in geneName) {
  #geneName <- "Or4"
  geneNumber <- sprintf('%0.3d', grep(paste0(geneName,"\\b"), receptorList))
  
  expression.cutoff <- 1
  
  dat1 <- antenna.neuron
  totalcells1 <- dim(antenna.neuron@meta.data)[1]
  
  ############################################################################################################### Heatmap 1, no recluster
  
  expr <- FetchData(object = dat1[[dat.assay]], vars = geneName)
  dat1 <- dat1[, which(x = expr > expression.cutoff)]
  subsetcellnumber <- dim(dat1@meta.data)[1]
  print(paste0(toString(subsetcellnumber), " cells expressing ",geneName))
  
  heatmap.title1 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Original Clusters",
                          "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  numberClusters <- length(levels(Idents(dat1)))
  clusterList <- levels(Idents(dat1))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  allgenes <- unique(allgenes)
  allgenes.og <- allgenes
  
  rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)
  
  
  ############################################################################################################### Heatmap 2, reclustered
  
  dat1.cl <- dat1
  
  if (grepl("harm*", pipeline)) {
    #  dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #  print("Variable features imported from 'integrated' assay.")
    reduction <- "harmony"
    dat.assay<-"SCT"
    
  }
  
  if (grepl("anch*", pipeline)) {
    #dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #print("Variable features imported from 'integrated' assay.")
    
    dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    reduction <- "pca"
    
    DefaultAssay(dat1.cl) <- "integrated"
    
    dat1.cl    <- RunPCA(dat1.cl, npcs = 50, verbose = F)
  }
  
  recluster.code <- "("
  
  try({
    dat1.cl <- RunUMAP(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] UMAP successful"))
    recluster.code <- append(recluster.code, "U")
  })
  
  try({
    dat1.cl <- RunTSNE(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] tSNE successful"))
    recluster.code <- append(recluster.code, "S")
  })
  
  
  try({
    dat1.cl <- FindNeighbors(dat1.cl, assay=dat.assay, 
                             dims = 1:50, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindNeighbors successful"))
    recluster.code <- append(recluster.code, "F")
  })
  
  try({
    dat1.cl <- FindClusters(dat1.cl, resolution = 1, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindClusters successful"))
    recluster.code <- append(recluster.code, "C")
  })
  
  
  heatmap.title2 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Reclustered ",paste(unlist(recluster.code),collapse=""),")",
                          "\n",toString(dim(dat1.cl@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  
  numberClusters <- length(levels(Idents(dat1.cl)))
  clusterList <- levels(Idents(dat1.cl))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1.cl,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  #allgenes <- append(allgenes, allgenes.og)
  allgenes <- unique(allgenes)
  
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, features = allgenes.og, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)
  
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1.cl, features = allgenes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("14a_",geneNumber, "_", geneName,"_",pipeline,"_rawHeatmaps.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
  
  rm(g)
  
  rm(u, cluster.subset, i, numberClusters, filename, cl.ae, clusterList)
  
  
  #..............................................................................................UMI violin
  
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  p4 <- VlnPlot(dat1.cl, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  
  g = arrangeGrob(p3,p4, ncol = 2)
  filename <- paste0("14b_",geneNumber, "_", geneName,"_",pipeline,"_UMIviolin.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000, g)
  
  #DimPlot(dat1.cl, reduction = "umap", group.by = 'dataset') + 
  #  DimPlot(dat1.cl, reduction = "tsne", group.by = 'dataset')
  
}


#.....................................................................................
# FIGURE PRESENTATION
# ....................................................................................

dat2 <- dat1
my_levels <- c("4" , "22", "10" ,"0",  "1" , "3" ,"5",  "6",  "8", "11",
               "12", "13", "14", "15" ,"16" ,"17" ,"21","25", "26" ,"27" ,"28",
               "29", "30", "31", "34", "35", "38")

dat2@active.ident <- factor(x = dat2@active.ident, levels = my_levels)
levels(Idents(dat2))

heatmap.genes <- c("LOC5575210", "Orco" ,"Ir25a","Ir76b" ,
                   "Ir8a"      , "Or36" ,"Or87","Or88","Or52","Or72", 
                   "Or84")


p5 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat2, features = heatmap.genes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)

p5
filename <- paste0("14d_",geneNumber, "_", geneName,"_",pipeline,"_figurePresentation.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000,p5)





######################################################################################
# Or47
######################################################################################

#.....................................................................................
# SET UP
# ....................................................................................
setwd()
rm(list = ls())
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
receptorList <- readRDS("Antenna_receptorList.rds")
pipeline <- "harmony_rmclusts_res4"

dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolder <- "Antenna_Figures3_NeuronExampleHeatmaps_Output"


geneName <- "Or47"
geneNumber <- grep(paste0(geneName,"\\b"), receptorList)
print(geneName)

dat1 <- antenna.neuron
expression.cutoff <- 1
totalcells1 <- dim(antenna.neuron@meta.data)[1]

#.....................................................................................
# GENE OF INTEREST
# ....................................................................................

for (geneName in geneName) {
  #geneName <- "Or4"
  geneNumber <- sprintf('%0.3d', grep(paste0(geneName,"\\b"), receptorList))
  
  expression.cutoff <- 1
  
  dat1 <- antenna.neuron
  totalcells1 <- dim(antenna.neuron@meta.data)[1]
  
  ############################################################################################################### Heatmap 1, no recluster
  
  expr <- FetchData(object = dat1[[dat.assay]], vars = geneName)
  dat1 <- dat1[, which(x = expr > expression.cutoff)]
  subsetcellnumber <- dim(dat1@meta.data)[1]
  print(paste0(toString(subsetcellnumber), " cells expressing ",geneName))
  
  heatmap.title1 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Original Clusters",
                          "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  numberClusters <- length(levels(Idents(dat1)))
  clusterList <- levels(Idents(dat1))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  allgenes <- unique(allgenes)
  allgenes.og <- allgenes
  
  rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)
  
  
  ############################################################################################################### Heatmap 2, reclustered
  
  dat1.cl <- dat1
  
  if (grepl("harm*", pipeline)) {
    #  dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #  print("Variable features imported from 'integrated' assay.")
    reduction <- "harmony"
    dat.assay<-"SCT"
    
  }
  
  if (grepl("anch*", pipeline)) {
    #dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #print("Variable features imported from 'integrated' assay.")
    
    dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    reduction <- "pca"
    
    DefaultAssay(dat1.cl) <- "integrated"
    
    dat1.cl    <- RunPCA(dat1.cl, npcs = 50, verbose = F)
  }
  
  recluster.code <- "("
  
  try({
    dat1.cl <- RunUMAP(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] UMAP successful"))
    recluster.code <- append(recluster.code, "U")
  })
  
  try({
    dat1.cl <- RunTSNE(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] tSNE successful"))
    recluster.code <- append(recluster.code, "S")
  })
  
  
  try({
    dat1.cl <- FindNeighbors(dat1.cl, assay=dat.assay, 
                             dims = 1:50, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindNeighbors successful"))
    recluster.code <- append(recluster.code, "F")
  })
  
  try({
    dat1.cl <- FindClusters(dat1.cl, resolution = 1, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindClusters successful"))
    recluster.code <- append(recluster.code, "C")
  })
  
  
  heatmap.title2 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Reclustered ",paste(unlist(recluster.code),collapse=""),")",
                          "\n",toString(dim(dat1.cl@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  
  numberClusters <- length(levels(Idents(dat1.cl)))
  clusterList <- levels(Idents(dat1.cl))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1.cl,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  #allgenes <- append(allgenes, allgenes.og)
  allgenes <- unique(allgenes)
  
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, features = allgenes.og, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)
  
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1.cl, features = allgenes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("14a_",geneNumber, "_", geneName,"_",pipeline,"_rawHeatmaps.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
  
  rm(g)
  
  rm(u, cluster.subset, i, numberClusters, filename, cl.ae, clusterList)
  
  
  #..............................................................................................UMI violin
  
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  p4 <- VlnPlot(dat1.cl, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  
  g = arrangeGrob(p3,p4, ncol = 2)
  filename <- paste0("14b_",geneNumber, "_", geneName,"_",pipeline,"_UMIviolin.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000, g)
  
  #DimPlot(dat1.cl, reduction = "umap", group.by = 'dataset') + 
  #  DimPlot(dat1.cl, reduction = "tsne", group.by = 'dataset')
  
}


#.....................................................................................
# FIGURE PRESENTATION
# ....................................................................................

dat2 <- dat1
my_levels <- c("19"  ,"32" ,"0"  ,"1",  "2"  ,"3" , "4" , "5"  ,"6" , "7",  "8",  "9",  "10" ,"11",
               "13" ,"14" ,"16", "17", "18",  "23" ,"25" ,"29", "33", "35",
               "37" ,"40")

dat2@active.ident <- factor(x = dat2@active.ident, levels = my_levels)
levels(Idents(dat2))

heatmap.genes <-  c("LOC5575210", "Orco"  ,"Ir25a" ,"Ir76b",
                    "Ir8a"  ,"Or47" ,"Or4"    ,
                    "Or71"  ,"Or82"  ,"Or3"    ,
                    "Or50","Gr15")

heatmap.title3 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                        "\n","Method: ",toString(pipeline), ", Original Clusters",
                        "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(totalcells1)," total neurons")

p5 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat2, features = heatmap.genes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title3)

p5
filename <- paste0("14d_",geneNumber, "_", geneName,"_",pipeline,"_figurePresentation.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000,p5)







######################################################################################
# Or84
######################################################################################

#.....................................................................................
# SET UP
# ....................................................................................
setwd()
rm(list = ls())
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
receptorList <- readRDS("Antenna_receptorList.rds")
pipeline <- "harmony_rmclusts_res4"

dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolder <- "Antenna_Figures3_NeuronExampleHeatmaps_Output"


geneName <- "Or84"
geneNumber <- grep(paste0(geneName,"\\b"), receptorList)
print(geneName)

dat1 <- antenna.neuron
expression.cutoff <- 1
totalcells1 <- dim(antenna.neuron@meta.data)[1]

#.....................................................................................
# GENE OF INTEREST
# ....................................................................................

for (geneName in geneName) {
  #geneName <- "Or4"
  geneNumber <- sprintf('%0.3d', grep(paste0(geneName,"\\b"), receptorList))
  
  expression.cutoff <- 1
  
  dat1 <- antenna.neuron
  totalcells1 <- dim(antenna.neuron@meta.data)[1]
  
  ############################################################################################################### Heatmap 1, no recluster
  
  expr <- FetchData(object = dat1[[dat.assay]], vars = geneName)
  dat1 <- dat1[, which(x = expr > expression.cutoff)]
  subsetcellnumber <- dim(dat1@meta.data)[1]
  print(paste0(toString(subsetcellnumber), " cells expressing ",geneName))
  
  heatmap.title1 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Original Clusters",
                          "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  numberClusters <- length(levels(Idents(dat1)))
  clusterList <- levels(Idents(dat1))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  allgenes <- unique(allgenes)
  allgenes.og <- allgenes
  
  rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)
  
  
  ############################################################################################################### Heatmap 2, reclustered
  
  dat1.cl <- dat1
  
  if (grepl("harm*", pipeline)) {
    #  dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #  print("Variable features imported from 'integrated' assay.")
    reduction <- "harmony"
    dat.assay<-"SCT"
    
  }
  
  if (grepl("anch*", pipeline)) {
    #dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    #print("Variable features imported from 'integrated' assay.")
    
    dat1.cl@assays[["SCT"]]@var.features <- dat1.cl@assays[["integrated"]]@var.features
    reduction <- "pca"
    
    DefaultAssay(dat1.cl) <- "integrated"
    
    dat1.cl    <- RunPCA(dat1.cl, npcs = 50, verbose = F)
  }
  
  recluster.code <- "("
  
  try({
    dat1.cl <- RunUMAP(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] UMAP successful"))
    recluster.code <- append(recluster.code, "U")
  })
  
  try({
    dat1.cl <- RunTSNE(dat1.cl, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0(geneName,"[",geneNumber,"] tSNE successful"))
    recluster.code <- append(recluster.code, "S")
  })
  
  
  try({
    dat1.cl <- FindNeighbors(dat1.cl, assay=dat.assay, 
                             dims = 1:50, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindNeighbors successful"))
    recluster.code <- append(recluster.code, "F")
  })
  
  try({
    dat1.cl <- FindClusters(dat1.cl, resolution = 1, verbose = F)
    print(paste0(geneName,"[",geneNumber,"] FindClusters successful"))
    recluster.code <- append(recluster.code, "C")
  })
  
  
  heatmap.title2 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                          "\n","Method: ",toString(pipeline), ", Reclustered ",paste(unlist(recluster.code),collapse=""),")",
                          "\n",toString(dim(dat1.cl@meta.data)[1])," out of ",toString(totalcells1)," total neurons")
  
  
  numberClusters <- length(levels(Idents(dat1.cl)))
  clusterList <- levels(Idents(dat1.cl))
  allgenes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a", geneName)
  
  for (i in 1:numberClusters) {
    #i <- 1
    u <- clusterList[i]
    
    cluster.subset <- subset(dat1.cl,idents = u)
    #print(dim(cluster.subset@meta.data)[1])
    
    if (dim(cluster.subset@meta.data)[1] >= 1) {
      cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
    }
    
    if (dim(cluster.subset@meta.data)[1] > 1) {
      suppressMessages({
        cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
      cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
    }
    
    cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
    cl.ae <- as.data.frame(cl.ae)                             #convert to df      
    cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
    cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MC genes                     
    cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
    
    cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
    cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
    cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  #allgenes <- append(allgenes, allgenes.og)
  allgenes <- unique(allgenes)
  
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, features = allgenes.og, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)
  
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1.cl, features = allgenes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("14a_",geneNumber, "_", geneName,"_",pipeline,"_rawHeatmaps.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
  
  rm(g)
  
  rm(u, cluster.subset, i, numberClusters, filename, cl.ae, clusterList)
  
  
  #..............................................................................................UMI violin
  
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  p4 <- VlnPlot(dat1.cl, assay=dat.assay, slot="counts", features = geneName, pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  
  g = arrangeGrob(p3,p4, ncol = 2)
  filename <- paste0("14b_",geneNumber, "_", geneName,"_",pipeline,"_UMIviolin.pdf")
  ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000, g)
  
  #DimPlot(dat1.cl, reduction = "umap", group.by = 'dataset') + 
  #  DimPlot(dat1.cl, reduction = "tsne", group.by = 'dataset')
  
}


#.....................................................................................
# FIGURE PRESENTATION
# ....................................................................................

dat2 <- dat1

heatmap.genes <-  c("LOC5575210", "Orco"  ,"Ir25a" ,"Ir76b",
                    "Ir8a"  , "Or84" ,      "Or87" ,      "Ir41b"  ,   
                    "Or85")

my_levels <- c(   "1", "21",
                  "0" ,  "2" , "3",  "4" , "5",  "6" , "7",  "8",  "10", "11", "12",
                  "13", "14" ,"15", "16" ,"17", "18" ,"19", "20",  "22", "24", "25",
                  "26", "28" ,"29", "30" ,"31", "33" ,"34", "35", "38" ,"40" ,"41")

dat2@active.ident <- factor(x = dat2@active.ident, levels = my_levels)
levels(Idents(dat2))


heatmap.title3 = paste0("Cells with ", toString(geneName)," above normalized expression > ", toString(expression.cutoff),
                        "\n","Method: ",toString(pipeline), ", Original Clusters",
                        "\n",toString(dim(dat1@meta.data)[1])," out of ",toString(totalcells1)," total neurons")

p5 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat2, features = heatmap.genes, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title3)

p5
filename <- paste0("14d_",geneNumber, "_", geneName,"_",pipeline,"_figurePresentation.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 4000,p5)







######################################################################################
# ALL CLUSTERS
######################################################################################
#.....................................................................................
# SET UP
# ....................................................................................
setwd()
rm(list = ls())
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"
receptorList <- readRDS("Antenna_receptorList.rds")
pipeline <- "harmony_rmclusts_res4"

dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolder <- "Antenna_Figures3_NeuronExampleHeatmaps_Output"

#.....................................................................................
# TOP 5 GENES PER CLUSTER
# ....................................................................................

dat <- antenna.neuron
numberClusters <- length(levels(Idents(dat)))
clusterList <- levels(Idents(dat))
heatmap.genes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a")

for (i in 1:numberClusters) {
  #i <- 1
  u <- clusterList[i]
  
  cluster.subset <- subset(dat,idents = u)
  print(paste0("cluster ",toString(u)," has ",dim(cluster.subset@meta.data)[1]," cells"))
  
  if (dim(cluster.subset@meta.data)[1] >= 1) {
    cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
  }
  
  if (dim(cluster.subset@meta.data)[1] > 1) {
    
    suppressMessages({
      cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
    cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
  }
  
  cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
  cl.ae <- as.data.frame(cl.ae)                             #convert to df      
  cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
  cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MT genes                     
  cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
  
  cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
  cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
  cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
  cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
  #cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
  #print(head(cl.ae$Gene, 5))
  
  heatmap.genes <- append(heatmap.genes, head(cl.ae$Gene, 5))
  
}
heatmap.genes <- unique(heatmap.genes)
rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)

heatmap.title1 = paste("All clusters, Method: ",pipeline,
                       "\n",toString(dim(antenna.neuron@meta.data)[1])," total neurons",sep="")

p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = antenna.neuron, 
                features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)
p1

filename <- paste0("10a_heatmaps_",pipeline,"_allclusters_top5genespercluster.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 3000,p1)

rm(dat, filename, heatmap.title1, p1)



#.....................................................................................
# TOP 7 GENES PER CLUSTER
# ....................................................................................

dat <- antenna.neuron
numberClusters <- length(levels(Idents(dat)))
clusterList <- levels(Idents(dat))
heatmap.genes <- c("LOC5575210", "Orco", "Ir25a", "Ir76b", "Ir8a")

for (i in 1:numberClusters) {
  #i <- 1
  u <- clusterList[i]
  
  cluster.subset <- subset(dat,idents = u)
  print(paste0("cluster ",toString(u)," has ",dim(cluster.subset@meta.data)[1]," cells"))
  
  if (dim(cluster.subset@meta.data)[1] >= 1) {
    cl.ae <- as.matrix(cluster.subset@assays[[dat.assay]]@data)
  }
  
  if (dim(cluster.subset@meta.data)[1] > 1) {
    
    suppressMessages({
      cl.ae <- AverageExpression(cluster.subset, assays=dat.assay)})                          #get avg expression of all genes in cluster      
    cl.ae <- cl.ae[[dat.assay]]                                        #pull out only relevant data                               
  }
  
  cl.ae <- cl.ae[order(cl.ae[,1],decreasing=TRUE),]        #sort from largest to smallest                          
  cl.ae <- as.data.frame(cl.ae)                             #convert to df      
  cl.ae <- tibble::rownames_to_column(cl.ae, "Gene")        #make gene names a column                           
  cl.ae <- cl.ae[!grepl("LOC*", cl.ae$Gene),]              #get rid of LOC & MT genes                     
  cl.ae <- cl.ae[!grepl("MT*", cl.ae$Gene),]
  
  cl.ae <- cl.ae[!grepl("Orco", cl.ae$Gene),]              #Delete co-receptors from main variable                     
  cl.ae <- cl.ae[!grepl("Ir25a", cl.ae$Gene),]
  cl.ae <- cl.ae[!grepl("Ir76b", cl.ae$Gene),]
  cl.ae <- cl.ae[!grepl("Ir8a", cl.ae$Gene),]
  #cl.ae <- cl.ae[!grepl(geneName, cl.ae$Gene),]
  #print(head(cl.ae$Gene, 5))
  
  heatmap.genes <- append(heatmap.genes, head(cl.ae$Gene, 7))
  
}
heatmap.genes <- unique(heatmap.genes)
rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)

heatmap.title1 = paste("All clusters, Method: ",pipeline,
                       "\n",toString(dim(antenna.neuron@meta.data)[1])," total neurons",sep="")

p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = antenna.neuron, 
                features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                draw.lines=FALSE, lines.width=4, size=4) +  
  theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)
p1

filename <- paste0("10b_heatmaps_",pipeline,"_allclusters_top7genespercluster.pdf")
ggsave(filename,path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 3000,p1)

rm(dat, filename, heatmap.title1, p1)


#.....................................................................................
# FIGURE PRESENTATION
# ....................................................................................
dat2 <- antenna.neuron

my_levels <- c( "1" ,"3"  ,"4" ,"22","6" ,"7" ,"8"  ,"10","11",
                "13","18","34", "14" ,"19","20","21",
                "24" ,"37" ,"39" ,"40",
                "27" ,"32", "12","41","9" ,"5","38",
                "0" ,"17","30","31","33","35", "28",
                
                "15" ,"16","23","29","36",
                "2","25", "26")

dat2@active.ident <- factor(x = dat2@active.ident, levels = my_levels)

heatmap.title1 = paste("All clusters, Method: ",pipeline,
                       "\n",toString(dim(antenna.neuron@meta.data)[1])," total neurons",sep="")


heatmap.genes <- c("LOC5575210","Orco"   ,"Ir25a"  ,"Ir76b"  ,
                   "Ir8a"      ,"Or84"   ,"Or87"   ,
                   
                   "Or38","Or25","Or33", "Or29",  "Or36", 
                   
                   "Or88"     ,   "Or16"  , 
                   "Or41"     , "Or42"  ,  "Or63"  , 
                   "Or64"     , "Or45"  ,  "Or122" , 
                   "Or121"    , "Or91"  , "Or130" , "Or116" , 
                   "Or52"     , "Or72"  , "Or132" , "Or10"  , 
                   "Or94"  , 
                   "Or97"  ,"Or100",  "Gr77", "Or110" , "Or103" , "Or115" , "Or114" , 
                   "Or125"    ,
                   "Or47"  ,"Or4"   , "Or71"  , "Or82"  ,
                   
                   "Or69"   , "Or85"  , "Or23"  , "Or6", 
                   
                   "Or79"     ,"Or27"  , "Or32"  ,"Or78"  ,  "Or26"  ,"Or77", "Or76",
                   
                   "Or44"     , "Or133" , "Or67"  , "Or43"  , 
                   "Or59"  , "Or31"  ,
                   "Or3"      ,"Gr15"  , "Or50"  ,   "Or11"  ,
                   "Or113" , "Or112","Or111"    , "Or104" , "Or105" ,
                   "Ir41l"    , "Or80"  , "Or81"  ,"Ir41m" ,
                   "Ir41k","Ir41p","Ir41a", "Ir41o" ,
                   "Ir41j"    ,"Ir41c" , "Ir41f" ,"Ir41b" ,  "Ir41e" , "Ir87a1"   , "Ir87a2",
                   "Ir93a"    , "Ir21a" , "Ir31a1", "Ir31a2", 
                   "Ir75k"    ,"Ir75l", "Ir75b" , "Ir75c" , "Ir75a" , 
                   "Ir75h"    , "Ir75g" , "Ir75f" , 
                   "Ir75e"    , "Ir75i" ,"Ir64a", "Ir75d"
)



p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat2, 
                features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                draw.lines=FALSE, lines.width=4, size=2, angle=0) +  
  theme(axis.text.y = element_text(size = 4)) + 
  scale_fill_gradientn(colors = c("black", "yellow", "red"))+ ggtitle(heatmap.title1)
p1

filename <- paste0("10c_allclustersheatmap_withlegend.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 8000, height = 5000,p1)

p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat2, 
                features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                draw.lines=FALSE, lines.width=4, size=1) +  
  theme(axis.text.y = element_text(size = 4)) + NoLegend()+
  scale_fill_gradientn(colors = c("black", "yellow", "red"))

filename <- paste0("10d_allclustersheatmap_nolegend.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "in", width = 8.5, height = 8,p1)

p3 <- VlnPlot(antenna.neuron, assay="SCT", slot="data", ncol=1,
              features = c("Orco","Ir25a","Ir76b","Ir8a"), pt.size = 0) + 
  scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()+
  geom_violin(position = position_dodge(0.8), color = NA)
p3
filename <- paste0("10e_allclustersviolin.pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "in", width = 8.5, height = 6,p3)





