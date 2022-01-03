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

### 4. Neuron subcluster

normalizationMethod='LogNormalize'
# normalizationMethod='SCTransform' # not good for visualization
print(c('normalizationMethod: ', normalizationMethod))

receptorFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/3.checkReceptors')
neuronFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/4.neuronClusters')
dir.create(neuronFolderPath, showWarnings = FALSE)
# picFolderPath=paste0(neuronFolderPath, '/pic')
# dir.create(picFolderPath, showWarnings = FALSE)


#####################################

load( file = file.path('Analysis/10282021_MaxPalp', 'maxpalp_3.neuron.RData') )
picFolderPath=paste0(neuronFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)

print(neuronFolderPath)
print(picFolderPath)


### Plot UMAP/tSNE

### (1)
sratDF <- NormalizeData(sratDF)
sratDF <- FindVariableFeatures(sratDF, selection.method = "vst", nfeatures = 2000)
sratDF <- ScaleData(sratDF, vars.to.regress = c('nCount_RNA'))
sratDF <- RunPCA(sratDF, npcs = 50, verbose = F)

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



srat.neuron <- NormalizeData(srat.neuron)
srat.neuron <- FindVariableFeatures(srat.neuron, selection.method = "vst", nfeatures = 2000)
srat.neuron <- ScaleData(srat.neuron, vars.to.regress = c('nCount_RNA'))

srat.neuron    <- RunPCA(srat.neuron, npcs = 50, verbose = F)
DimPlot(srat.neuron, reduction = "pca")
ElbowPlot(srat.neuron, ndims = 50)


# VariableFeatures
top10 <- head(VariableFeatures(srat.neuron), 10)
# VariableFeatures(srat.neuron)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(srat.neuron)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


srat.neuron    <- RunUMAP(srat.neuron, dims = 1:10, verbose = F)
srat.neuron    <- RunTSNE(srat.neuron, dims = 1:10, verbose = F)
srat.neuron    <- FindNeighbors(srat.neuron, dims = 1:10, verbose = F)
srat.neuron    <- FindClusters(srat.neuron, verbose = T, resolution = c(0.5, 0.8, 2, 4)[4])

head(srat.neuron@meta.data)


pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_tSNE_', 'neuronOnly.pdf')),
    width=10, height=5)
print( ( DimPlot(srat.neuron, reduction = "umap", label = TRUE) + theme(legend.position = "none")) + 
        ( DimPlot(srat.neuron, reduction = "tsne", label = TRUE) + theme(legend.position = "none")) )
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('nFeature_nCount_ptMt_ploidy.pdf')),
    width=12, height=3.5)
print(VlnPlot(sratDF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)) # ploidy difference in cluster #2?
dev.off()

print(VlnPlot(srat.neuron, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))


save.image(file = file.path(neuronFolderPath, 'neuralClusters.RData') )








################################################################################################################################
load(file = file.path(neuronFolderPath, 'neuralClusters.RData'))

DefaultAssay(srat.neuron) <- "RNA"

### coreceptor
pdf(file = file.path(picFolderPath, 
                     paste0('coreceptor_FeaturePlot_neuronOnly_batchCorrect.pdf')),
    width=15, height=15)
print(FeaturePlot(srat.neuron, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Gr3'), ncol = 3))
dev.off()


print(FeaturePlot(sratDF, reduction = 'tsne', features = c('Or8', 'Or49', 'Ir100a', 'Ir93a'), ncol = 2))
print(FeaturePlot(srat.neuron, reduction = 'tsne', features = c('Or8', 'Or49', 'Ir100a', 'Ir93a'), ncol = 2))



  # Receptor list: for plotting all receptors
allGeneID.df<- as.character(read.table('SeuratFile/features.tsv', header = FALSE)$V1)
head(allGeneID.df)

receptorL=c()
lapply(allGeneID.df, function(x) {
  # print(x[[1]])
  if ( isFALSE(startsWith(x, 'LOC'))  ) {
    if ( isFALSE(startsWith(x, 'MT'))  ) {
      receptorL <<- append(receptorL, x)
    }
  }
})

# Plot receptor 
receptor.dot <- DotPlot(srat.neuron, features = receptorL) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
# receptor.dot
receptor.dot.value <- receptor.dot$data %>% as.data.frame()


receptor.dot.value.noCoR %>%
  ggplot(aes(avg.exp, pct.exp)) + 
  geom_point()+scale_x_log10()+
  geom_hline(yintercept=c(25,50), color=scatterColors[3])+
  geom_vline(xintercept=1, color=scatterColors[3])

receptor.dot.value %>%
  ggplot(aes(pct.exp)) + 
  geom_density()+scale_x_log10()
receptor.dot.value %>%
  ggplot(aes(avg.exp)) + 
  geom_density()+scale_x_log10()


### avg.exp vs. pct.exp: scatter plot 
library(scales)
hex_codes1 <- hue_pal()(10)  
scatterColors <- c('#66c2a5', '#fc8d62', '#8da0cb')

pmain <- receptor.dot.value.noCoR %>%
  ggplot(aes(avg.exp, pct.exp)) + 
  geom_point()+scale_x_log10()+
  geom_hline(yintercept=c(25,50), color=scatterColors[3])+
  geom_vline(xintercept=c(1,5), color=scatterColors[3])
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = receptor.dot.value.noCoR, aes(x = avg.exp), fill='grey')+
  geom_vline(xintercept = c(1,5), color=scatterColors[3])+scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = receptor.dot.value.noCoR, aes(x = pct.exp), fill='grey') + 
  geom_vline(xintercept = c(25,50), color=scatterColors[3]) +
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Orco_Ir25a_normalizedExp.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()





srat.neuron@meta.data %>%
  ggplot(aes(Orco_Exp)) + 
  geom_density()


receptor.dot.value.max25 %>%
  filter(features.plot == 'Orco' ) %>%
  ggplot(aes(avg.exp, color=id)) + 
  geom_density()

class(receptor.dot.value.max25$id)








receptor.dot.value.noCoR <- receptor.dot.value %>% filter(! features.plot %in% c('Orco', 'Ir25a', 'Ir76b', 'Ir8a'))
coReceptorL <- c('Orco', 'Ir25a', 'Ir76b', 'Ir8a')

receptor.dot.value.noCoR <- receptor.dot.value.noCoR %>% arrange(desc(avg.exp))
receptorExpSortL <- unique(receptor.dot.value.noCoR$features.plot)
receptorExpSortL <- c(as.character(receptorExpSortL), coReceptorL)
head(receptorExpSortL)
tail(receptorExpSortL)

clusterExpSortL <- unique(receptor.dot.value.noCoR$id) %>% as.character()
# clusterExpSortL


# tree
# BiocManager::install("treeio")
# devtools::install_github("YuLab-SMU/ggtree")
mat <- receptor.dot.value %>%
  select(-avg.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = pct.exp) %>% as.data.frame()
row.names(mat) <- mat$id
row.names(mat)
mat <- mat[,-1] #drop gene column as now in rows

library(data.table)
mat.T <- transpose(mat)
# t(mat)
colnames(mat.T) <- rownames(mat)
rownames(mat.T) <- colnames(mat)





clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot
clust$labels
v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
v_clust
ddgram_col <- as.dendrogram(v_clust)
ggtree_plot_col <- ggtree::ggtree(ddgram_col) + ggtree::layout_dendrogram()
ggtree_plot_col



# id = factor(id, levels = clust$labels[clust$order])
# factor(features.plot, levels = v_clust$labels[v_clust$order] )

dotplot <- receptor.dot.value %>%
  mutate(`% Expressing` = pct.exp,
         features.plot = factor(features.plot, levels =  receptorExpSortL),
         id = factor(id, levels = clusterExpSortL)
  ) %>% 
  ggplot(aes(x=features.plot, y = id, color = avg.exp, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,20), oob = scales::squish, name = 'avg.exp')+
  scale_y_discrete(position = "right", limits = rev)
dotplot

# png(filename = 'Analysis/Olivia_10062021/plot/dotplot_allReceptors.png', width = 7000, height = 1800, res = 150)
# dotplot
# dev.off()


### 25% of the cell

mat.T$max.pct <- apply(mat.T, 1, function(input){
  # print(max(as.numeric(input)))
  return(max(as.numeric(input)))
})
mat.T$geneName <- row.names(mat.T)
receptroMax25L <- mat.T %>% filter(max.pct >= 25) %>% pull(geneName)


receptorExpSort_max25L <- receptorExpSortL[receptorExpSortL %in% receptroMax25L]
receptor.dot.value.max25 <- receptor.dot.value %>% filter(features.plot %in% receptorExpSort_max25L)

# no: limits = c(0,4)
# oob = scales::squish,
dotplot <- receptor.dot.value.max25 %>%
  mutate(`% Expressing` = pct.exp,
         features.plot = factor(features.plot, levels =  receptorExpSortL),
         id = factor(id, levels = clusterExpSortL)
  ) %>% 
  ggplot(aes(x=features.plot, y = id, color = avg.exp, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('') +ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'avg.exp', oob = scales::squish, limits = c(0,20))+
  scale_y_discrete(position = "right", limits = rev)
dotplot

dotplot <- receptor.dot.value.max25 %>%
  mutate(`% Expressing` = pct.exp,
         features.plot = factor(features.plot, levels =  receptorExpSortL),
         id = factor(id, levels = clusterExpSortL)
  ) %>% 
  ggplot(aes(x=features.plot, y = id, color = avg.exp, size = `% Expressing`)) + 
  geom_point() + 
  theme_half_open() +
  background_grid() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('') +ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(5,20))+
  scale_y_discrete(position = "right", limits = rev) +
  scale_size(limits = c(50, 100)) 

pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_25receptor_batchCorrect.pdf')),
    width=20, height=11)
print(dotplot)
aqwdev.off()


# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[1:6], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[7:12], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[13:18], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[19:24], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[25:30], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[31:36], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[37:42], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[43:48], ncol = 3)
# FeaturePlot(srat.neuron, reduction = 'tsne', features = receptorExpSort_max25L[49:54], ncol = 3)
# 
# FeaturePlot(orcoLow.seurat, reduction = 'tsne',features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
# FeaturePlot(orcoLow.seurat, reduction = 'tsne',features = c('nFeature_RNA', 'Orco'), ncol = 2)
# FeaturePlot(orcoLow.seurat, reduction = 'tsne',features = c('nCount_RNA', 'Orco'), ncol = 2)
# FeaturePlot(orcoLow.seurat, reduction = 'tsne',features = c('Orco', 'Ir25a'), ncol = 2)
# VlnPlot(orcoLow.seurat, features = c('Orco', 'Ir25a', 'Ir76b'), ncol = 3)



##########################################################
### Plotting scatter plot

library(scales)
hex_codes1 <- hue_pal()(10)  
scatterColors <- c('#66c2a5', '#fc8d62', '#8da0cb')

### Orco vs. Ir25a
head(srat.neuron@assays$RNA@counts[c('Orco','Ir25a'),])

orcoIl25a.mtx <- srat.neuron@assays$RNA@counts[c('Orco','Ir25a'),]

srat.neuron$Orco_UMIs <- srat.neuron@assays$RNA@counts['Orco',]
srat.neuron$Ir25a_UMIs <- srat.neuron@assays$RNA@counts['Ir25a',]

srat.neuron$Orco_Exp <- srat.neuron@assays$RNA@data['Orco',]
srat.neuron$Ir25a_Exp <- srat.neuron@assays$RNA@data['Ir25a',]

pmain <- srat.neuron@meta.data %>%
  ggplot( aes(Orco_Exp, Ir25a_Exp) ) + 
  geom_point(size=0.3) +
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 2, color=scatterColors[3])
  # coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = srat.neuron@meta.data, aes(x = Orco_Exp), fill=scatterColors[1])+
  geom_vline(xintercept = 2, color=scatterColors[3])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = srat.neuron@meta.data, aes(x = Ir25a_Exp), fill=scatterColors[2]) + 
  geom_vline(xintercept = 1, color=scatterColors[3]) +
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Orco_Ir25a_normalizedExp.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()

### Orco vs. Ir25a in cluster level

receptor.dot.value.max25 %>%
  filter(features.plot %in% c('Orco', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() %>%
  ggplot(aes(Orco, Ir25a, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0,) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()

library(ggrepel)
receptor.dot.value.max25 %>%
  filter(features.plot %in% c('Orco', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() %>%
  filter(Orco > 50 & Ir25a > 10)


# + scale_x_log10() + scale_y_log10()


### Plotting all coreceptors

p0=FeatureScatter(object = srat.neuron, feature1 = 'Orco', feature2 = 'Ir25a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p1=FeatureScatter(object = srat.neuron, feature1 = 'Orco', feature2 = 'Ir76b') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p2=FeatureScatter(object = srat.neuron, feature1 = 'Orco', feature2 = 'Ir8a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p3=FeatureScatter(object = srat.neuron, feature1 = 'Ir25a', feature2 = 'Ir76b') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p4=FeatureScatter(object = srat.neuron, feature1 = 'Ir25a', feature2 = 'Ir8a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p5=FeatureScatter(object = srat.neuron, feature1 = 'Ir76b', feature2 = 'Ir8a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_allCoreceptors_normalizedExp.pdf')),
    width=15, height=10)
print((p0 + p1 + p2) /
  (p3 + p4 + p5))
dev.off()

# plotting coexpression of receptors
p0=FeatureScatter(object = srat.neuron, feature1 = 'Ir41f', feature2 = 'Or103') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p1= FeatureScatter(object = srat.neuron, feature1 = 'Ir41f', feature2 = 'Ir75h') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p2=FeatureScatter(object = srat.neuron, feature1 = 'Or4', feature2 = 'Ir31a1') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p3=FeatureScatter(object = srat.neuron, feature1 = 'Or4', feature2 = 'Or47') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p4=FeatureScatter(object = srat.neuron, feature1 = 'Or84', feature2 = 'Or85') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p5=FeatureScatter(object = srat.neuron, feature1 = 'Ir41k', feature2 = 'Ir41j') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_receptors_normalizedExp.pdf')),
    width=15, height=10)
print((p0 + p1 + p2) /
        (p3 + p4 + p5))
dev.off()



##########################################################
### chord plot
library(chorddiag)

highReceptorRNA = srat.neuron@assays$RNA@counts[receptorExpSortL[1:20],]
# str(highReceptorRNA)
# diff(highReceptorRNA@p)
# highReceptorRNA != 0
# summary(highReceptorRNA)$j


highReceptorRNA.expCol=highReceptorRNA[,colSums(highReceptorRNA) != 0]
# for (recptor in rownames(highReceptorRNA.expCol)) {
# }


# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = 20, nrow = 20))
colnames(coexp.df) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df) <- rownames(highReceptorRNA.expCol)

coexp.df2 <- data.frame(matrix(ncol = 20, nrow = 20))
colnames(coexp.df2) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df2) <- rownames(highReceptorRNA.expCol)


for (i in 1:20) {
  for (j in 1:20) {
    # print(c('i = ', i))
    # print(c('j = ', j))
    if (i==j) { coexp.df[i,j] <- 0 }
    else {
      lapply(1:length(ncol(highReceptorRNA.expCol)), function(x){
        print(c('i = ', i))
        print(c('j = ', j))
        iExpressedCells.mx = highReceptorRNA.expCol[,highReceptorRNA.expCol[i,] >5]
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] > 5])
        # print(jExpressedN)
        coexp.df[i,j]<<-jExpressedN
        
        iExpressedCells.mx2 = highReceptorRNA.expCol[,highReceptorRNA.expCol[i,] >10]
        jExpressedN2=NCOL(iExpressedCells.mx2[, iExpressedCells.mx2[j,] > 10])
        # print(jExpressedN)
        coexp.df2[i,j]<<-jExpressedN2
      })
    }
  }
}

library(scales)
hex_codes1 <- hue_pal()(20)  
# names(hex_codes1) <- colnames(coexp.df)



print(chorddiag(as.matrix(coexp.df), groupColors = hex_codes1, groupnamePadding = 20))
chorddiag(as.matrix(coexp.df2), groupColors = hex_codes1, groupnamePadding = 20)


##########################################################
### Grouping co-receptor

library(ComplexHeatmap)
receptor.dot.value.max25 = receptor.dot.value.max25[receptorExpSortL,]


mat.avg.exp <- receptor.dot.value %>%
  select(-avg.exp.scaled, -pct.exp) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame()
row.names(mat.avg.exp) <- mat.avg.exp$id
row.names(mat.avg.exp)
mat.avg.exp <- mat.avg.exp[,-1] #drop gene column as now in rows

mat.avg.exp.scaled <- receptor.dot.value %>%
  select(-avg.exp, -pct.exp) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp.scaled) %>% as.data.frame()
row.names(mat.avg.exp.scaled) <- mat.avg.exp.scaled$id
row.names(mat.avg.exp.scaled)
mat.avg.exp.scaled <- mat.avg.exp.scaled[,-1] #drop gene column as now in rows

mat.pct.exp <- receptor.dot.value %>%
  select(-avg.exp.scaled, -avg.exp) %>%
  pivot_wider(names_from = features.plot, values_from = pct.exp) %>% as.data.frame()
row.names(mat.pct.exp) <- mat.pct.exp$id
row.names(mat.pct.exp)
mat.pct.exp <- mat.pct.exp[,-1] #drop gene column as now in rows


Heatmap(as.matrix(mat.pct.exp[,coReceptorL]))
Heatmap(as.matrix(mat.avg.exp[,coReceptorL]))
Heatmap(as.matrix(mat.avg.exp.scaled[,coReceptorL]))


Heatmap(as.matrix(mat.pct.exp[,coReceptorL[1:3]]), cluster_columns = FALSE)

Heatmap(as.matrix(transpose(mat.avg.exp[,coReceptorL])))


# Try value range
Heatmap(as.matrix(mat.avg.exp[,coReceptorL[1:3]]),
        name = "Average\nexpression",
        col = circlize::colorRamp2(c(0,25, 50), c("blue", "white", "red")))

Heatmap(scale( as.matrix(mat.avg.exp[,coReceptorL[1:3]]) ),
        name = "Average\nexpression",
        )
Heatmap(scale( as.matrix(mat.avg.exp[,coReceptorL[1:3]]) ),
        name = "Average\nexpression",
        col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
Heatmap(as.matrix(mat.avg.exp.scaled[,coReceptorL[1:3]]),
        name = "Scaled\naverage\nexpression",
)
Heatmap(as.matrix(mat.avg.exp.scaled[,coReceptorL[1:3]]),
        name = "Scaled\naverage\nexpression",
        col = circlize::colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
)


save.image(file = file.path('Analysis/10282021_MaxPalp', 'maxpalp_4.neuralClusters.RData') )



