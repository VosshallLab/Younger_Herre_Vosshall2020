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
# CELLS BY CLUSTER: TOP CHEMORECEPTORS IN EACH CLUSTER
######################################################################################

rm(list = ls())   #remove all objects
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")
pipeline <- "harmony_rmclusts_res4"
dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

outputFolderHeatmap <- "10b_byCluster_heatmaps_avggenelist_reclustered_res4"
dir.create(outputFolderHeatmap)
outputFolderVLN <- "10c_byCluster_violinplot_topgenes_reclustered_res4"
dir.create(outputFolderVLN)
outputFolderBatchCheck <- "10d_byCluster_batchcheck_dimensionreduction"
dir.create(outputFolderBatchCheck)


dat <- antenna.neuron
numberClusters <- length(levels(Idents(dat)))
clusterList <- levels(Idents(dat))

for (i in 1:numberClusters) {
  #i <- 42
  u <- clusterList[i]
  
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
    print(paste0("Cluster ",clusterList[i]," UMAP successful"))
    recluster.code <- append(recluster.code, "U")
  })
  
  try({
    dat1 <- RunTSNE(dat1, dims = 1:50, verbose = F, reduction = reduction)
    print(paste0("Cluster ",clusterList[i]," tSNE successful"))
    recluster.code <- append(recluster.code, "S")
  })
  
  try({
    dat1 <- FindNeighbors(dat1, assay=dat.assay, 
                             dims = 1:50, verbose = F)
    print(paste0("Cluster ",clusterList[i]," FindNeighbors successful"))
    recluster.code <- append(recluster.code, "F")
  })
  
  try({
    dat1 <- FindClusters(dat1, resolution = 4, verbose = F)
    print(paste0("Cluster ",clusterList[i]," FindClusters successful"))
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
  
  #Top genes, total average expression in entire cluster
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                  features = average.genes, label=TRUE, group.bar.height =0.01, 
                  draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
    ggtitle(heatmap.title1)

  #Top genes, average expression of subcluster
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, 
                  features = heatmap.genes, label=TRUE, group.bar.height =0.01, 
                  draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ 
    ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("10b_heatmaps_cluster",toString(u),".pdf")
  ggsave(filename, path=outputFolderHeatmap, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)

  #Violin plot of top genes (entire cluster) within subclusters
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = average.genes[1:12], pt.size = 1) + scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  filename <- paste0("10b_heatmaps_cluster",toString(u),"_violinUMI_",dat.assay,".pdf")
  ggsave(filename, path=outputFolderVLN, limitsize = FALSE, units = "px", width = 4000, height = 4000, p3)
  
  #Seeing distribution of cells within cluster by batch
  p5 <- DimPlot(dat1, reduction = "tsne", group.by = 'dataset')+ggtitle(paste0("cluster ",toString(u)," by batch"))
  p6 <- DimPlot(dat1, reduction = "umap", group.by = 'dataset')+ggtitle(paste0("cluster ",toString(u)," by batch"))
  g = arrangeGrob(p5,p6, ncol = 2)
  filename <- paste0("10d_heatmaps_cluster",toString(u),"_batchcheck.pdf")
  ggsave(filename, path=outputFolderBatchCheck, limitsize = FALSE, units = "px", width = 4000, height = 2000, g)
  
}

rm(list=setdiff(ls(), "antenna.neuron"))


######################################################################################
# CELLS BY CHEMORECEPTOR
######################################################################################

rm(list = ls())   #remove all objects
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")
pipeline <- "harmony_rmclusts_res4"
dat.assay <- "SCT"
DefaultAssay(antenna.neuron) <- dat.assay

#.....................................................................................
# SET UP
# ....................................................................................


outputGeneHeatmaps <- '11a_byGene_heatmaps'
outputGeneVLN <- '11b_byGene_violinplot'

dir.create(outputGeneHeatmaps)
dir.create(outputGeneVLN)

dat.assay <- 'SCT'
DefaultAssay(antenna.neuron) <- dat.assay

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

for (geneName in receptorList[18:191]) {

  #geneName <- "Or47"
  geneNumber <- sprintf('%0.3d', grep(paste0(geneName,"\\b"), receptorList))
  expression.cutoff <- 1
  dat1 <- antenna.neuron
  totalcells1 <- dim(antenna.neuron@meta.data)[1]
  
  ################################################################## Heatmap 1, no recluster
  
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
    cl.ae <- cl.ae[!grepl((paste0(geneName,"\\b")), cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  allgenes <- unique(allgenes)
  allgenes.originalclusters <- allgenes
  
  rm(u, cluster.subset, i, numberClusters, cl.ae, allgenes, expr, clusterList)
  
  
  ################################################################## Heatmap 2, reclustered
  
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
    cl.ae <- cl.ae[!grepl((paste0(geneName,"\\b")), cl.ae$Gene),]
    #head(cl.ae$Gene, 20)
    
    allgenes <- append(allgenes, head(cl.ae$Gene, 5))
    #print(head(cl.ae))
    
  }
  allgenes.reclustered <- unique(allgenes)
  
  #All cells expressing chemoreceptor above threshold, original clusters
  p1 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1, features = allgenes.originalclusters, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title1)

  #All cells expressing chemoreceptor above threshold, reclustered  
  p2 <- DoHeatmap(slot = "data", assay = dat.assay, object = dat1.cl, features = allgenes.reclustered, label=TRUE, group.bar.height =0.01, draw.lines=FALSE, lines.width=4, size=4) +  
    theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = c("black", "yellow", "red"), na.value = "black")+ ggtitle(heatmap.title2)
  
  #p1/p2
  g = arrangeGrob(p1,p2, ncol = 1)
  filename <- paste0("11a_",geneNumber, "_", geneName,"_",pipeline,".pdf")
  ggsave(filename, path=outputGeneHeatmaps, limitsize = FALSE, units = "px", width = 4000, height = 4000,g)
  
  rm(g)
  
  rm(u, cluster.subset, i, numberClusters, filename, cl.ae, allgenes, clusterList)
  
  
  #..............................................................................................UMI violin
  
  #violin plot of sctransform-adjusted UMIs on for original clusters
  p3 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = geneName, 
                pt.size = 1) + ggtitle(paste0(geneName,", all neurons \n (sctransform-adjusted UMIs)")) + NoLegend() #+ scale_y_continuous(breaks=seq(0,50,1))
  p4 <- VlnPlot(dat1, assay=dat.assay, slot="counts", features = geneName, 
                pt.size = 1) + ggtitle(paste0(geneName," cells above SCT value ",expression.cutoff,", original clusters \n (sctransform-adjusted UMIs)")) + NoLegend() #+ scale_y_continuous(breaks=seq(0,50,1))
  p5 <- VlnPlot(dat1.cl, assay=dat.assay, slot="counts", features = geneName, 
                pt.size = 1)+ ggtitle(paste0(geneName," cells above SCT value ",expression.cutoff,", reclustered \n (sctransform-adjusted UMIs)"))+ NoLegend() #+ scale_y_continuous(breaks=seq(0,50,1)) + NoLegend()
  
  #p3+p4+p5
  g = arrangeGrob(p3,p4,p5, ncol = 1)
  filename <- paste0("11c_",geneNumber, "_", geneName,"_",pipeline,"_UMIviolin.pdf")
  ggsave(filename, path=outputGeneVLN, limitsize = FALSE, units = "px", width = 4000, height = 6000, g)
  
  #DimPlot(dat1.cl, reduction = "umap", group.by = 'dataset') + 
  #  DimPlot(dat1.cl, reduction = "tsne", group.by = 'dataset')
  
}
