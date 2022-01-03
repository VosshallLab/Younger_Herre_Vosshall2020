suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(ComplexHeatmap)
  library(viridis)
  library(dichromat)
  library(RColorBrewer)
})

theme_set(theme_cowplot())
projectFolder='/Users/tzuchiao/ProjectAnalysis/Leslie_Vosshall'
setwd(projectFolder)


#####################################

### 7. expressed receptors

normalizationMethod='LogNormalize'
# normalizationMethod='SCTransform' # not good for visualization
print(c('normalizationMethod: ', normalizationMethod))

neuronFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/6.neuronClusters')

################################################################################################################################
load(file = file.path('Analysis/integratedData_decontX', 'antenna_6.neuralClusters.RData') )

receptorFolderPath=paste0('Analysis/integratedData_decontX/', normalizationMethod, '/7.expressedReceptors')
dir.create(receptorFolderPath, showWarnings = FALSE)
picFolderPath=paste0(receptorFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)

DefaultAssay(neuron.batchCorrect) <- "integrated"


### Check nFeature_RNA mtTx in the cluster level
# Idents(seurat.combined) <- seurat.combined$seurat_clusters
p2 = VlnPlot(neuron.batchCorrect, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = 'seurat_clusters', pt.size = 0)

pdf(file = file.path(picFolderPath, 
                     paste0('integrate_nFeature_nCount_ptMt_', normalizationMethod, '_afterFiler.pdf')),
    width=10, height=5)
print(p2)
dev.off()



### Main tissue types:
set.seed(96)

seurat.combined <- ScaleData(seurat.combined, vars.to.regress = c('nCount_RNA'))

seurat.combined$annotation <- apply(seurat.combined@meta.data, 1, function(input){
  if (input['seurat_clusters'] %in% neuronClusterIDs) {return('Neuron')}
  else if (input['seurat_clusters'] %in% c('35', '55')) {return('Glia')}
  else if (input['seurat_clusters'] %in% c('6', '14', '15', '23')) {return('Epithelia')}
  else {return('Others')}
})
seurat.combined@meta.data$annotation <- factor(seurat.combined@meta.data$annotation, 
                                               levels = c('Others', 'Epithelia', 'Glia', 'Neuron'))

mainTissueGeneL <- c('LOC5564305','LOC5570152','LOC5576812','LOC110678282','LOC110676862','LOC5566721','LOC5564848','LOC5570381','LOC5565901','LOC5570204','LOC5567355')


mainTissueMat <- seurat.combined[["RNA"]]@data[mainTissueGeneL, ] %>% as.matrix()
quantileL <- quantile(mainTissueMat, c(0.1, 0.93))
# quantileL <- append(quantileL, mean(quantileL), after = 1)
# col_fun = circlize::colorRamp2(quantileL, colorRampPalette(brewer.pal(9, "YlGnBu"))(3) ) 
col_fun = circlize::colorRamp2(quantileL, c('blue', 'yellow') ) 


cluster_anno<- seurat.combined@meta.data$annotation
left_anno <- c(
               replicate(3, 'Epithelia Markers'),
               replicate(3, 'Glia Markers'),
               replicate(5, 'Neuron Markers')
)

left_anno <- factor(left_anno,
                    levels = c('Epithelia Markers', 
                               'Glia Markers','Neuron Markers'))


pdf(file = file.path(picFolderPath, 
                     paste0('heatmap_allCells_mainMarkers.pdf')),
    width=35, height=8)
print(Heatmap(mainTissueMat, name = "Expression",
              column_split = cluster_anno, # column_split = factor(cluster_anno, levels(seurat.combined@meta.data$annotation2)),
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
                  foo = anno_block(gp = gpar(fill = scales::hue_pal()( length( levels(seurat.combined@meta.data$annotation)) ) ))
                ),
              
              # left_annotation =
              #   rowAnnotation(
              #     foo = anno_block(gp = gpar(fill = scales::hue_pal()( length(levels(left_anno)) ) )),
              #     ),
              # 
))
dev.off()


### Exclude others

seurat.combined.noOthers <- subset(seurat.combined, subset = annotation != 'Others')

mainTissueMat <- seurat.combined.noOthers[["RNA"]]@data[mainTissueGeneL, ] %>% as.matrix()


quantileL <- quantile(mainTissueMat, c(0.1, 0.93))
quantileL <- append(quantileL, mean(quantileL), after = 1)
col_fun = circlize::colorRamp2(quantileL, colorRampPalette(brewer.pal(9, "YlGnBu"))(3) ) 
col_fun = circlize::colorRamp2(quantileL, c('blue', 'yellow') ) 


cluster_anno<- seurat.combined.noOthers@meta.data$annotation
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


pdf(file = file.path(picFolderPath, 
                     paste0('heatmap_allCells_mainMarkers_noOthers_blue-yellow.pdf')),
    width=35, height=8)
print(Heatmap(mainTissueMat, name = "Expression",
              column_split = cluster_anno, # column_split = factor(cluster_anno, levels(seurat.combined@meta.data$annotation2)),
              column_title_rot = 0,
              height = unit(0.75, "npc"),
              
              col = col_fun, show_column_names = FALSE, use_raster = FALSE,
              
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 18),
              
              row_names_gp = gpar(fontsize = 10),
              column_gap = unit(0.5, "mm"),
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              cluster_row_slices = TRUE,
              
              
              row_split = left_anno,
              row_title_gp = gpar(fontsize = 18),
              row_title_rot = 0,
              
              
              top_annotation = 
                HeatmapAnnotation(
                  foo = anno_block(gp = gpar(fill = scales::hue_pal()( length( levels(seurat.combined@meta.data$annotation)) ) ))
                ),
              # left_annotation =
              #   rowAnnotation(
              #     foo = anno_block(gp = gpar(fill = scales::hue_pal()( length(levels(left_anno)) ) )),
              #     ),
              # 
))
dev.off()





##############################
### all main markers

# mainMarkerGenes <- unique(c('LOC5576180','LOC5567584','LOC5564497','LOC5572971','LOC5570013','LOC5570015','LOC5570873','LOC5575988','LOC5563869','LOC5566583','LOC110678401','LOC5569395','LOC5577119','LOC5569912','LOC5563727','LOC5571968','LOC5577691','LOC110676950','LOC5580054','LOC5564903','LOC5579912','LOC5573368','LOC5572260','LOC5568897','LOC5576779','LOC5566412','LOC5564863','LOC110674032','LOC5573137','LOC5568897','LOC110679166','LOC5567149','LOC110678184','LOC5576990','LOC5576977','LOC5575889','LOC5566583','LOC110678401','LOC110675128','MT-AaegL5-0439288','LOC5565553','LOC5563786','LOC5563869','LOC110678401','LOC5566583','LOC5570013','LOC5578493','LOC5568491','MT-AaegL5-0439282','LOC5569243','LOC5569853','LOC5569665','LOC5570081','LOC5576180','LOC5574769','LOC110678523','LOC5564265','LOC5569083','LOC5564263','LOC5569822','LOC5578664','LOC5565410','LOC5579336','LOC5575140','LOC110676586','LOC5574907','LOC5564305','LOC5570152','LOC5576812','LOC5570721','LOC5564305','LOC5574686','LOC5578194','LOC5565277','LOC5568005','LOC5564305','LOC5568514','LOC5574271','LOC5565260','LOC5568988','LOC5564305','LOC5573135','LOC5564386','LOC5567363','LOC5570162','LOC110678282','LOC110676862','LOC5566721','LOC5566438','LOC5571030','LOC5579799','LOC5571097','LOC5572489','LOC110678282','LOC5564848','LOC5570381','LOC5565901','LOC5570204','LOC5567355')
#                           )# no receptor gene
seurat.combined$annotation2 <- apply(seurat.combined@meta.data, 1, function(input){
  if (input['seurat_clusters'] %in% neuronClusterIDs) {return('Neuron')}
  else if (input['seurat_clusters'] %in% c('35', '55')) {return('Glia')}
  else if (input['seurat_clusters'] %in% c('6', '14', '15', '23')) {return('Epithelia')}
  else {return(input['seurat_clusters'])}
})
# factor(seurat.combined@meta.data$annotation2 )
# levels(factor(seurat.combined@meta.data$annotation2 ))

seurat.combined@meta.data$annotation2 <- 
  factor(cluster_anno<- seurat.combined@meta.data$annotation2,
         levels = c( 
           as.character(as.numeric(levels(factor(seurat.combined@meta.data$annotation2))[1:17])), 
           c("Epithelia", "Glia", "Neuron" )))
# levels(factor(seurat.combined@meta.data$annotation2 ))
# length(levels(factor(seurat.combined@meta.data$annotation2 )))



# mainMarkerMat <- seurat.combined[["RNA"]]@data[mainMarkerGenes, ] %>% as.matrix()
# mainMarkerMat<- t(scale(t(mainMarkerMat)))
mainMarkerMat <- seurat.combined[["RNA"]]@data[mainTissueGeneL, ] %>% as.matrix()
mainMarkerMat<- t(scale(t(mainMarkerMat)))

quantileL <- quantile(mainMarkerMat, c(0.10, 0.95))
quantileL
# write.table(mainMarkerMat,file="mainMarkerGenes.mtx")


col_fun = circlize::colorRamp2(quantileL, c("black", "yellow")) # for scaled values
colours=colorRampPalette( c('black', 'yellow'))(24)



cluster_anno<- seurat.combined@meta.data$annotation2
left_anno <- c(
  replicate(6, 'Cluster Markers'),
  replicate(3, 'Epithelia Markers'),
  replicate(3, 'Glia Markers'),
  replicate(5, 'Neuron Markers')
)
left_anno <- factor(left_anno, 
                    levels = c('Cluster Markers', 'Epithelia Markers', 'Glia Markers', 'Neuron Markers'))

pdf(file = file.path(picFolderPath, 
                     paste0('heatmap_allCells_mainMarkers_withClusterMarkers_reorder.pdf')),
    width=35, height=10)
print(Heatmap(mainMarkerMat, name = "Expression",
        column_split = cluster_anno, # column_split = factor(cluster_anno, levels(seurat.combined@meta.data$annotation2)),
        column_title_rot = 90,
        height = unit(0.75, "npc"),
        
        # row_order = mainMarkerGenes,
        row_order = mainTissueGeneL, 
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        
        row_split = left_anno,
        row_title_gp = gpar(fontsize = 14),
        row_title_rot = 0,
        
        top_annotation = 
          HeatmapAnnotation(
            foo = anno_block(gp = gpar(fill = scales::hue_pal()( length( levels(seurat.combined@meta.data$annotation2)) ) )),
            annotation_name_rot=45
            ),
        use_raster = FALSE,
        # col = c(col_fun, colours)[[2]],
        col = col_fun,
        # col = colours,
        show_column_names = FALSE))
dev.off()



### markers in each cluster
### Find markers
library(future)
plan("multiprocess", workers = 4)

neuron.markers <- FindAllMarkers(neuron.batchCorrect, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(merge.markers, file = file.path(receptorFolderPath, 'neuronMarkers.csv'))


##############################
### coreceptor featureplot
pdf(file = file.path(picFolderPath, 
                     paste0('coreceptor_FeaturePlot_neuronOnly_batchCorrect.pdf')),
    width=10, height=10)
print(FeaturePlot(neuron.batchCorrect, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b'), ncol = 2, pt.size = 0.5))
dev.off()


##############################
### scatter plot & DotPlot

DefaultAssay(seurat.combined) <- "RNA"
DefaultAssay(neuron.batchCorrect) <- "RNA"
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
receptor.dot <- DotPlot(neuron.batchCorrect, features = receptorL) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
# receptor.dot
receptor.dot.value <- receptor.dot$data %>% as.data.frame()

coReceptorL <- c('Orco', 'Ir25a', 'Ir76b', 'Ir8a')

receptor.dot.value.noCoR <- receptor.dot.value %>% filter(! features.plot %in% coReceptorL)
receptor.dot.value.noCoR <- receptor.dot.value.noCoR %>% arrange(desc(avg.exp))

receptorExpSortL <- unique(receptor.dot.value.noCoR$features.plot)
receptorExpSortL <- c(as.character(receptorExpSortL), coReceptorL)
head(receptorExpSortL)
tail(receptorExpSortL)

clusterExpSortL <- unique(receptor.dot.value.noCoR$id) %>% as.character()
# clusterExpSortL

### matrix values
mat <- receptor.dot.value %>%
  select(-avg.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = pct.exp) %>% as.data.frame()
row.names(mat) <- mat$id
row.names(mat)
mat <- mat[,-1] #drop gene column as now in rows

mat_avg.exp <- receptor.dot.value %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame()
row.names(mat_avg.exp) <- mat_avg.exp$id
row.names(mat_avg.exp)
mat_avg.exp <- mat_avg.exp[,-1] #drop gene column as now in rows


mat.T <- transpose(mat)
# t(mat)
colnames(mat.T) <- rownames(mat)
rownames(mat.T) <- colnames(mat)

mat.T$max.pct <- apply(mat.T, 1, function(input){
  # print(max(as.numeric(input)))
  return(max(as.numeric(input)))
})
mat.T$geneName <- row.names(mat.T)
head(mat.T)

#######################################
### avg.exp vs. pct.exp: scatter plot 
hex_codes1 <- hue_pal()(20)
scatterColors <- c('#66c2a5', '#fc8d62', '#8da0cb')

# relaxed criteria
pmain <- receptor.dot.value.noCoR %>%
  ggplot(aes(avg.exp, pct.exp)) + 
  geom_point(size=0.5)+scale_x_log10()+
  geom_hline(yintercept=c(25,50)[1], color=scatterColors[3])+
  geom_vline(xintercept=c(1,5)[1], color=scatterColors[3]) +
  xlab('Average expression') + ylab('Percentage of expressing cell')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = receptor.dot.value.noCoR, aes(x = avg.exp), fill='grey')+
  geom_vline(xintercept = c(1,5)[1], color=scatterColors[3])+scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = receptor.dot.value.noCoR, aes(x = pct.exp), fill='grey') + 
  geom_vline(xintercept = c(25,50)[1], color=scatterColors[3]) +
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_avg.exp_pct.exp_allchemoreceptors_relaxed.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()

# stringent criteria
pmain <- receptor.dot.value.noCoR %>%
  ggplot(aes(avg.exp, pct.exp)) + 
  geom_point(size=0.5)+scale_x_log10()+
  geom_hline(yintercept=c(25, 35)[2], color=scatterColors[3])+
  geom_vline(xintercept=c(1,5)[1], color=scatterColors[3]) +
  xlab('Average expression') + ylab('Percentage of expressing cell')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = receptor.dot.value.noCoR, aes(x = avg.exp), fill='grey')+
  geom_vline(xintercept = c(1,5)[1], color=scatterColors[3])+scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = receptor.dot.value.noCoR, aes(x = pct.exp), fill='grey') + 
  geom_vline(xintercept = c(25,35)[2], color=scatterColors[3]) +
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_avg.exp_pct.exp_allchemoreceptors_stringent.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()


#######################################
### dotplot for chemreceptors

### No filter

dotplot <- receptor.dot.value %>%
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
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_y_discrete(position = "right", limits = rev) 

pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_receptor_batchCorrect_noFilter.pdf')),
    width=20, height=11)
print(dotplot)
dev.off()



## 25% & 35% of the cell
receptroMax25L <- mat.T %>% filter(max.pct >= 25) %>% pull(geneName)
receptorExpSort_max25L <- receptorExpSortL[receptorExpSortL %in% receptroMax25L]
receptor.dot.value.max25 <- receptor.dot.value %>% filter(features.plot %in% receptorExpSort_max25L)

receptroMax35L <- mat.T %>% filter(max.pct >= 35) %>% pull(geneName)
receptorExpSort_max35L <- receptorExpSortL[receptorExpSortL %in% receptroMax35L]
receptor.dot.value.max35 <- receptor.dot.value %>% filter(features.plot %in% receptorExpSort_max35L)


# relaxed criteria: 25% ==> show all signals
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
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_y_discrete(position = "right", limits = rev) 

pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_receptor_batchCorrect_relaxedCriteria_allSignals.pdf')),
    width=22, height=13)
print(dotplot)
dev.off()


# relaxed criteria: 25% ==> remove background
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
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_y_discrete(position = "right", limits = rev) +
  scale_size(limits = c(25, 100)) 

pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_receptor_batchCorrect_relaxedCriteria.pdf')),
    width=20, height=11)
print(dotplot)
dev.off()


# stringent criteria: 35%

receptor.dot.value.max35 <- receptor.dot.value.max35 %>%
  mutate(`% Expressing` = pct.exp,
         features.plot = factor(features.plot, levels =  receptorExpSort_max35L),
         id = factor(id, levels = clusterExpSortL))


clusterExpSortL_rmNoChemoL <- vector()
for (clusterID in levels(receptor.dot.value.max35$id) ) {
  clusterChemL <- receptor.dot.value.max35 %>%
    filter(`% Expressing` >=35 & avg.exp >= 1) %>%
    filter(id == clusterID & (! features.plot %in% coReceptorL)) %>%
    pull(features.plot) %>% as.character()
  # print(clusterChemL)
  
  if (length(clusterChemL) > 0) {
    clusterExpSortL_rmNoChemoL <<- c(clusterExpSortL_rmNoChemoL, clusterID)
  }
}


dotplot <- receptor.dot.value.max35 %>% 
  filter(`% Expressing` >=35 & avg.exp >= 1 & id %in% clusterExpSortL_rmNoChemoL) %>%
  ggplot(aes(x=features.plot, y = id, color = avg.exp, size = `% Expressing`)) + 
  geom_point() + 
  theme_half_open() +
  background_grid() +
  theme(axis.line  = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank()) +
  xlab('') +ylab('') +
  # scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_color_gradientn(colours = c('#008080', '#FF00FF'),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_y_discrete(position = "right", limits = rev) +
  scale_size(limits = c(35, 100)) 

pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_receptor_batchCorrect_stringent_rmNoChemoL_teal-magenta .pdf')),
    width=25, height=11)
print(dotplot)
dev.off()


### use tree order from Olivia's paper, y axis use clustering results

receptorTreeOrederL <- c(
  c('Or6','Or4','Or28','Or8','Or33','Or132N','Or11','Or2','Or10','Or9','Or38','Or31','Or30','Or47','Or32','Or51P','Or50','Or26','Or27','Or55','Or54','Or57P','Or48','Or3','Or15','Or14','Or71','Or49F','Or29','Or13','Or25','Or22','Or67','Or42','Or41','Or23','Or24','Or21','Or20','Or16','Or19','Or18P','Or45','Or64','Or63','Or133N','Or43','Or44','Or66','Or70','Or68P','Or69','Or40','Or34','Or37','Or59','Or128','Or62','Or73','Or74','Or58','Or75','Or76','Or77','Or79','Or78','Or52','Or72','Or36','Or61','Or60','Or85','Or84','Or80','Or81','Or82','Or94','Or95','Or97','Or96','Or88','Or87','Or86','Or90','Or89','Or90','Or89','Or91','Or93','Or92P','Or99','Or123','Or125','Or111F','Or113','Or112','Or105','Or104','Or100','Or115','Or114','Or103','Or110','Or101','Or102','Or106','Or108','Or107','Or120P','Or122','Or121','Or118','Or119','Or117','Or116'),
  c('Ir93a','Ir75k','Ir75lN','Ir64a','Ir75d','Ir75jNP','Ir75i','Ir75h','Ir75h','Ir75f','Ir75e','Ir31a2N','Ir31a1','Ir75a','Ir75c','Ir75b','Ir40a','Ir7a','Ir7d','Ir7c','Ir7b','Ir7i','Ir7k','Ir7j','Ir7f','Ir7e','Ir7h','Ir7g','Ir7l','Ir7rN','Ir7qN','Ir7pN','Ir7o','Ir7n','Ir7m','Ir41p','Ir41o','Ir41n','Ir41mNC','Ir41l','Ir41k','Ir41j','Ir41iC','Ir41hC','Ir41g','Ir41bC','Ir41a','Ir41f','Ir41e','Ir41dP','Ir41c','Ir68a','Ir21a','Ir60a','Ir87a1','Ir87a2','Ir104N','Ir103N','Ir101','Ir102N','Ir100d','Ir100c','Ir100bC','Ir100a','Ir106','Ir105NF','Ir107N','Ir108N','Ir110NP','Ir109N','Ir111N','Ir112NP','Ir117NF','Ir116N','Ir113NF','Ir114N','Ir115NP','Ir119','Ir125','Ir124','Ir123','Ir122N','Ir121','Ir120','Ir118N','Ir127N','Ir126N','Ir128C','Ir130N','Ir129','Ir132','Ir133N','Ir131N','Ir135NP','Ir134NP','Ir137','Ir136','Ir138','Ir139','Ir141','Ir140','Ir142N','Ir143N','Ir144','Ir146NP','Ir154N','Ir147N','Ir148C','Ir150','Ir149','Ir151','Ir154','Ir153','Ir152','Ir156N','Ir155N','Ir157N','Ir159N','Ir158N','Ir161NP','Ir160N','Ir162NC','Ir164NC','Ir163N','Ir166NP','Ir165N','Ir168N','Ir167NP','Ir172N','Ir170N','Ir169N'),
  c('Gr77')
)

length(unique(receptorTreeOrederL))

receptorTreeOrederL <- receptorTreeOrederL[receptorTreeOrederL %in% receptorExpSort_max35L] %>%
  unique()
receptorTreeOreder_wCoreceptorL <- c(receptorTreeOrederL, coReceptorL)

### use cell % to cluster the y axis, coreceptors are not included
mat35 = mat[rownames(mat) %in% clusterExpSortL_rmNoChemoL,
            (colnames(mat) %in% receptorExpSort_max35L & 
               (! colnames(mat) %in% coReceptorL))]
mat35 <- scale(mat35)
clust <- hclust(dist(mat35 %>% as.matrix())) # hclust with distance matrix
# colnames(mat)
ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot
# order of the tree: clust$labels[clust$order]

### use avg.exp to cluster the y axis, coreceptors are not included

mat35_avg.exp = mat_avg.exp[rownames(mat_avg.exp) %in% clusterExpSortL_rmNoChemoL,
            (colnames(mat_avg.exp) %in% receptorExpSort_max35L & 
               (! colnames(mat_avg.exp) %in% coReceptorL))]
mat35_avg.exp <- scale(mat35_avg.exp)



length(receptorTreeOreder_wCoreceptorL)
length(levels(receptor.dot.value.max35$features.plot))


receptor.dot.value.max35.tree <- receptor.dot.value.max35 %>%
  filter(id %in% clusterExpSortL_rmNoChemoL & features.plot %in% receptorTreeOreder_wCoreceptorL) %>%
  mutate(`% Expressing` = pct.exp,
         features.plot = factor(features.plot, levels =  receptorTreeOreder_wCoreceptorL),
         id = factor(id, levels = clust$labels[clust$order])) 

dotplot <- receptor.dot.value.max35.tree %>% 
  filter(id %in% clusterExpSortL_rmNoChemoL) %>%
  ggplot(aes(x=features.plot, y = id, color = avg.exp, size = `% Expressing`)) + 
  geom_point() + 
  theme_half_open() +
  background_grid() +
  theme(axis.line  = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank()) +
  xlab('') +ylab('') +
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_y_discrete(position = "right", limits = rev) +
  scale_size(limits = c(35, 100)) 
pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_receptor_batchCorrect_stringent_tree.pdf')),
    width=25, height=11)
print(dotplot)
dev.off()




### use cell % to cluster the y axis, coreceptors are not included
clust <- hclust(dist(mat35 %>% as.matrix())) # hclust with distance matrix
# colnames(mat)
ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot
# order of the tree: clust$labels[clust$order]
v_clust <- hclust(dist(mat35 %>% as.matrix() %>% t()))
v_clust
ddgram_col <- as.dendrogram(v_clust)
ggtree_plot_col <- ggtree::ggtree(ddgram_col) + ggtree::layout_dendrogram()
ggtree_plot_col

receptor.dot.value.max35.unsuptree <- receptor.dot.value.max35 %>%
  filter(id %in% clusterExpSortL_rmNoChemoL & features.plot %in% receptorExpSort_max35L) %>%
  mutate(`% Expressing` = pct.exp,
         features.plot = factor(features.plot, levels =  c(v_clust$labels[v_clust$order], coReceptorL)),
         id = factor(id, levels = clust$labels[clust$order])) 

length(receptorExpSort_max35L)
length(levels(receptor.dot.value.max35$features.plot))

dotplot <- receptor.dot.value.max35.unsuptree %>% 
  filter(id %in% clusterExpSortL_rmNoChemoL) %>%
  ggplot(aes(x=features.plot, y = id, color = avg.exp, size = `% Expressing`)) + 
  geom_point() + 
  theme_half_open() +
  background_grid() +
  theme(axis.line  = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank()) +
  xlab('') +ylab('') +
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_y_discrete(position = "right", limits = rev) +
  scale_size(limits = c(35, 100)) 
pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_receptor_batchCorrect_stringent_unsupervisedCluster.pdf')),
    width=25, height=11)
print(dotplot)
dev.off()





##########################################################
### Plotting scatter plot

### Orco vs. Ir25a

# head(neuron.batchCorrect@assays$RNA@counts[c('Orco','Ir25a'),])
# orcoIl25a.mtx <- neuron.batchCorrect@assays$RNA@counts[c('Orco','Ir25a'),]

neuron.batchCorrect$Orco_UMIs <- neuron.batchCorrect@assays$RNA@counts['Orco',]
neuron.batchCorrect$Ir25a_UMIs <- neuron.batchCorrect@assays$RNA@counts['Ir25a',]

neuron.batchCorrect$Orco_Exp <- neuron.batchCorrect@assays$RNA@data['Orco',]
neuron.batchCorrect$Ir25a_Exp <- neuron.batchCorrect@assays$RNA@data['Ir25a',]



### UMI: Orco vs. Ir25a
pmain <- neuron.batchCorrect@meta.data %>%
  ggplot( aes(Orco_UMIs, Ir25a_UMIs) ) + 
  geom_point(size=0.3) + 
  scale_x_log10() + scale_y_log10()
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = neuron.batchCorrect@meta.data, aes(x = Orco_UMIs), fill=scatterColors[1]) +scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = neuron.batchCorrect@meta.data, aes(x = Ir25a_UMIs), fill=scatterColors[2]) + scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

print(ggdraw(p2))

head(neuron.batchCorrect@meta.data)

### normalized exp: Orco vs. Ir25a (Orco cutoff 1)
pmain <- neuron.batchCorrect@meta.data %>%
  ggplot( aes(Orco_Exp, Ir25a_Exp) ) + 
  geom_point(size=0.3) +
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = neuron.batchCorrect@meta.data, aes(x = Orco_Exp), fill=scatterColors[1]) +
  geom_vline(xintercept = 1, color=scatterColors[3])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = neuron.batchCorrect@meta.data, aes(x = Ir25a_Exp), fill=scatterColors[2]) + 
  geom_vline(xintercept = 1, color=scatterColors[3])+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Orco_Ir25a_normalizedExp_Orco1.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()


### normalized exp: Orco vs. Ir25a (Orco cutoff 2)
pmain <- neuron.batchCorrect@meta.data %>%
  ggplot( aes(Orco_Exp, Ir25a_Exp) ) + 
  geom_point(size=0.3) +
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 2, color=scatterColors[3])
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = neuron.batchCorrect@meta.data, aes(x = Orco_Exp), fill=scatterColors[1]) +
  geom_vline(xintercept = 2, color=scatterColors[3])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = neuron.batchCorrect@meta.data, aes(x = Ir25a_Exp), fill=scatterColors[2]) + 
  geom_vline(xintercept = 1, color=scatterColors[3])+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Orco_Ir25a_normalizedExp_Orco2.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()





### Orco vs. Ir25a in cluster level
orcoIr25a_avgExp.clusters.df <- 
  receptor.dot.value %>%
  filter(features.plot %in% c('Orco', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  orcoIr25a_avgExp.clusters.df %>%
  ggplot(aes(Orco, Ir25a, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Orco') + ylab('Average expression of Ir25a')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = orcoIr25a_avgExp.clusters.df, aes(x = Orco), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = orcoIr25a_avgExp.clusters.df, aes(x = Ir25a), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Orco_Ir25a_normalizedExp_clusterLevel.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()


#################################
### Scatter plot: different combination


### Plotting all coreceptors

p0=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Orco', feature2 = 'Ir25a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p1=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Orco', feature2 = 'Ir76b') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p2=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Orco', feature2 = 'Ir8a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p3=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir25a', feature2 = 'Ir76b') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p4=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir25a', feature2 = 'Ir8a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p5=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir76b', feature2 = 'Ir8a') + 
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
p0=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir41f', feature2 = 'Or103') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p1= FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir41f', feature2 = 'Ir75h') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p2=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Or4', feature2 = 'Ir31a1') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p3=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Or4', feature2 = 'Or47') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p4=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Or84', feature2 = 'Or85') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p5=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir41k', feature2 = 'Ir41j') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_receptors_normalizedExp.pdf')),
    width=15, height=10)
print((p0 + p1 + p2) /
        (p3 + p4 + p5))
dev.off()


### 10262021 genes
p0=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Or4', feature2 = 'Or47') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p1= FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir41f', feature2 = 'Or103') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p2=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir41l', feature2 = 'Or82') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p3=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir41f', feature2 = 'Or84') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p4=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir31a1', feature2 = 'Or4') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_receptors1026_normalizedExp.pdf')),
    width=15, height=10)
print((p0 + p1 + p2) /
        (p3 + p4 + plot_spacer()))
dev.off()

neuron.batchCorrect@meta.data$chemoreceptorAnno
Idents(neuron.batchCorrect)
### 11012021 genes
p0=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir41l', feature2 = 'Or82',  clos = 'chemoreceptorAnno') + 
  geom_point() + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p1= FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir31a1', feature2 = 'Or82', clos = 'chemoreceptorAnno') + 
  geom_point() + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p2=FeatureScatter(object = neuron.batchCorrect, feature1 = 'Ir75g', feature2 = 'Or84',  clos = 'chemoreceptorAnno') + 
  geom_point() + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])


pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_receptors1101_normalizedExp.pdf')),
    width=15, height=10)
print((p0 + p1 + p2))
dev.off()




##########################################################
### chord plot
library(chorddiag)

highReceptorRNA = neuron.batchCorrect@assays$RNA@counts[receptorExpSortL[1:20],]
highReceptorRNA.data = neuron.batchCorrect@assays$RNA@data[receptorExpSortL[1:20],]


highReceptorRNA.expCol=highReceptorRNA[,colSums(highReceptorRNA) != 0]
highReceptorRNA.data.expCol=highReceptorRNA.data[,colSums(highReceptorRNA.data) != 0]
# for (recptor in rownames(highReceptorRNA.expCol)) {
# }
# highReceptorRNA.data.expCol

# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = 20, nrow = 20))
colnames(coexp.df) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df) <- rownames(highReceptorRNA.expCol)

coexp.df2 <- data.frame(matrix(ncol = 20, nrow = 20))
colnames(coexp.df2) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df2) <- rownames(highReceptorRNA.expCol)

coexp.df3 <- data.frame(matrix(ncol = 20, nrow = 20))
colnames(coexp.df3) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df3) <- rownames(highReceptorRNA.expCol)


# normalized expression
for (i in 1:20) {
  for (j in 1:20) {
    # print(c('i = ', i))
    # print(c('j = ', j))
    if (i==j) { coexp.df[i,j] <- 0 }
    else {
      lapply(1:length(ncol(highReceptorRNA.data.expCol)), function(x){
        print(c('i = ', i))
        print(c('j = ', j))
        iExpressedCells.mx = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1]
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        # print(jExpressedN)
        coexp.df[i,j]<<-jExpressedN
        
        iExpressedCells.mx2 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 2]
        jExpressedN2=NCOL(iExpressedCells.mx2[, iExpressedCells.mx2[j,] >= 2])
        # print(jExpressedN)
        coexp.df2[i,j]<<-jExpressedN2
        
        iExpressedCells.mx3 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1.5]
        jExpressedN3=NCOL(iExpressedCells.mx3[, iExpressedCells.mx3[j,] >= 1.5])
        # print(jExpressedN)
        coexp.df3[i,j]<<-jExpressedN3
      })
    }
  }
}



# names(hex_codes1) <- colnames(coexp.df)

print(chorddiag(as.matrix(coexp.df), groupColors = hex_codes1, groupnamePadding = 10, showTicks = FALSE))
chorddiag(as.matrix(coexp.df2), groupColors = hex_codes1, groupnamePadding = 20)
chorddiag(as.matrix(coexp.df3), groupColors = hex_codes1, groupnamePadding = 20)



### Use different cell# cutoff for plotting chord plot
### Remove coexpression <= 10 or 20


# groupnamePadding: shift of group name
chorddiag(as.matrix(coexp.df), groupColors = hex_codes1, 
          groupnamePadding = c( 10 ), 
          showTicks = TRUE, groupPadding=1, margin=120)

coexp.df.filter10<- coexp.df
coexp.df.filter10[coexp.df <= 10] <- 0
chorddiag(as.matrix(coexp.df.filter10), groupColors = hex_codes1, 
          groupnamePadding = c( replicate(18, 10), 50, 10 ), 
          showTicks = TRUE, groupPadding=1, margin=120, showZeroTooltips=FALSE, fadeLevel=1)

coexp.df.filter20 <- coexp.df
coexp.df.filter20[coexp.df <= 20] <- 0
chorddiag(as.matrix(coexp.df.filter20), groupColors = hex_codes1, 
          groupnamePadding = c( replicate(16, 10), 50, 10, 50, 10 ), 
          showTicks = TRUE, groupPadding=1, margin=120, showZeroTooltips=FALSE)
chorddiag(as.matrix(coexp.df.filter20), groupColors = hex_codes1, margin=120, showZeroTooltips=FALSE)


### Remove non-coexpressed genes
coexp.df.filter10 <- coexp.df
coexp.df.filter10[coexp.df <= 10] <- 0
# check genes with 0 coexpression
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter10)) {
  if (sum(coexp.df.filter10[,i]) == 0) {
    print(colnames(coexp.df.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter10 <- coexp.df.filter10[-zeroGenes_V, -zeroGenes_V]
chorddiag(as.matrix(coexp.df.filter10), groupColors = hue_pal()(ncol(coexp.df.filter10)), 
          groupnamePadding = 30, 
          showTicks = TRUE, groupPadding=1, margin=120)

coexp.df.filter20 <- coexp.df
coexp.df.filter20[coexp.df <= 20] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter20)) {
  if (sum(coexp.df.filter20[,i]) == 0) {
    print(colnames(coexp.df.filter20)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter20 <- coexp.df.filter20[-zeroGenes_V, -zeroGenes_V]
chorddiag(as.matrix(coexp.df.filter20), groupColors = hue_pal()(ncol(coexp.df.filter20)), 
          groupnamePadding = 30, 
          showTicks = TRUE, groupPadding=1, margin=120)





##########################################################
###  Histogram summarizing expressing receptors in each cell (? or cluster)

receptorRNA = neuron.batchCorrect@assays$RNA@counts[receptorExpSortL,]
# receptorRNA
receptorRNA.data = neuron.batchCorrect@assays$RNA@data[receptorExpSortL,]

# receptorRNA.t.df <- as.data.frame(t(as.matrix(receptorRNA)))



## randomly coexpressed genes

random_L <- c('Ir41f', 'Or103', 'Ir41f', 'Ir75h', 'Or4', 'Ir31a1')
hisGenes.df <- as.data.frame(
  matrix(ncol = 4, nrow = 3))
rownames(hisGenes.df) <- c('Gene A only', 'Gene B only', 'Gene A and B') #remove : 'Neither gene A nor B'
colnames(hisGenes.df)[4] <- 'Group'
hisGenes.df$Group <- rownames(hisGenes.df)
colN=0
# i=1
for (i in seq(1,6,2)) {
  colN = colN +1
  colnames(hisGenes.df)[colN] <- paste(random_L[i], random_L[i+1],sep = '\n')

  hisGenes.df[,colN] <- c(
    length(
      receptorRNA.data[
        random_L[i], receptorRNA.data[random_L[i],] >= 1 & receptorRNA.data[random_L[i+1],] < 1
      ]),
    length(
      receptorRNA.data[
        random_L[i], receptorRNA.data[random_L[i],] < 1 & receptorRNA.data[random_L[i+1],] >= 1
      ]),
    length(
      receptorRNA.data[
        random_L[i], receptorRNA.data[random_L[i],] >= 1 & receptorRNA.data[random_L[i+1],] >= 1
      ])
  )
}



p1 <- hisGenes.df %>%
  pivot_longer(cols=!Group, names_to = 'GroupName', values_to = 'Count') %>%
  group_by(GroupName) %>%
  mutate(GroupName= factor(GroupName, levels = colnames(hisGenes.df)[1:3])) %>%
  ggplot(aes(GroupName, Count, fill=Group, labels)) + geom_bar(stat="identity", position = "fill")+
  ylab('Percentage of neurons') + 
  scale_fill_brewer(palette = "Pastel1") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal_hgrid()+
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        plot.margin = unit(c(1,1,2,5), "lines")) +
  annotate("text", x = c(0.5, 1:3), y = -0.05, 
           label = c('A\nB', colnames(hisGenes.df)[1:3]))
p1

ggplot(aes(criterion, Ratio, fill=ExpGroup)) + geom_bar(stat="identity", position = "fill")+
  scale_y_continuous(labels = scales::percent_format())


## coexpressed genes
coExpGenes.df <- as.data.frame(
  matrix(ncol = 6, nrow = 3)
)

coExpGenesL <- c('Or4', 'Or47', 'Or84', 'Or85', 'Ir41k', 'Ir41j')
coExpGenes.df <- as.data.frame(
  matrix(ncol = 4, nrow = 3))
rownames(coExpGenes.df) <- c('Gene A only', 'Gene B only', 'Gene A and B') #remove : 'Neither gene A nor B'
colnames(coExpGenes.df)[4] <- 'Group'
coExpGenes.df$Group <- rownames(coExpGenes.df)
colN=0
# i=1
for (i in seq(1,6,2)) {
  colN = colN +1
  colnames(coExpGenes.df)[colN] <- paste(coExpGenesL[i], coExpGenesL[i+1],sep = '\n')
  
  coExpGenes.df[,colN] <- c(
    length(
      receptorRNA.data[ 
        coExpGenesL[i], receptorRNA.data[coExpGenesL[i],] >= 1 & receptorRNA.data[coExpGenesL[i+1],] < 1
      ]),
    length(
      receptorRNA.data[ 
        coExpGenesL[i], receptorRNA.data[coExpGenesL[i],] < 1 & receptorRNA.data[coExpGenesL[i+1],] >= 1
      ]),
    length(
      receptorRNA.data[ 
        coExpGenesL[i], receptorRNA.data[coExpGenesL[i],] >= 1 & receptorRNA.data[coExpGenesL[i+1],] >= 1
      ])
  )
}


p2 <- coExpGenes.df %>%
  pivot_longer(cols=!Group, names_to = 'GroupName', values_to = 'Count') %>%
  group_by(GroupName) %>%
  mutate(GroupName= factor(GroupName, levels = colnames(coExpGenes.df)[1:3])) %>%
  ggplot(aes(GroupName, Count, fill=Group, labels)) + geom_bar(stat="identity", position = "fill")+
  ylab('Percentage of neurons') + 
  scale_fill_brewer(palette = "Pastel1") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal_hgrid()+
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        plot.margin = unit(c(1,1,2,5), "lines")) +
  annotate("text", x = c(0.5, 1:3), y = -0.05, 
           label = c('A\nB', colnames(coExpGenes.df)[1:3]))

p2



## coexpressed coreceptors genes
corecptorExpGenes.df <- as.data.frame(
  matrix(ncol = 7, nrow = 4))
rownames(corecptorExpGenes.df) <- c('Gene A only', 'Gene B only', 'Gene A and B', 'Neither gene A nor B')
colnames(corecptorExpGenes.df)[7] <- 'Group'
corecptorExpGenes.df$Group <- rownames(corecptorExpGenes.df)
colN=0
for (i in 1: (length(coReceptorL)-1)){
  for (j in 1:(length(coReceptorL)-i)){
    if (j <= length(coReceptorL)){
      colN =colN + 1 
      # print(c(i,j, colN))
      geneA=coReceptorL[i]
      geneB=coReceptorL[i+j]
      # print(c(geneA, geneB))
      abColumnName <- paste(geneA, geneB, sep = '\n')
      colnames(corecptorExpGenes.df)[colN] <- abColumnName
      
      corecptorExpGenes.df[,colN] <-
        c(
          length(receptorRNA.data[geneA, receptorRNA.data[geneA,] >= 1 & receptorRNA.data[geneB,] < 1]),
          length(receptorRNA.data[geneA, receptorRNA.data[geneA,] < 1 & receptorRNA.data[geneB,] >= 1]),
          length(receptorRNA.data[geneA, receptorRNA.data[geneA,] >= 1 & receptorRNA.data[geneB,] >= 1]),
          length(receptorRNA.data[geneA, receptorRNA.data[geneA,] < 1 & receptorRNA.data[geneB,] < 1])
        )
      
    }
  }
}


p3 <- corecptorExpGenes.df %>%
  pivot_longer(cols=!Group, names_to = 'GroupName', values_to = 'Count') %>%
  group_by(GroupName) %>%
  mutate(GroupName= factor(GroupName, levels = colnames(corecptorExpGenes.df)[1:6])) %>%
  ggplot(aes(GroupName, Count, fill=Group, labels)) + geom_bar(stat="identity", position = "fill")+
  ylab('Percentage of neurons') + 
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal_hgrid()+
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(1,1,3,5), "lines")) +
  annotate("text", x = c(0.5, 1:6), y = -0.05, 
           label = c('A\nB', colnames(corecptorExpGenes.df)[1:6]))



pdf(file = file.path(picFolderPath, 
                     paste0('barGraph_percentageOfCoexpression.pdf')),
    width=20, height=6)
print(p1 + p2 + p3 + plot_layout(widths = c(1, 1, 2)))
dev.off()








corecptorExpGenes.df[1,1:2] <-  c('Orco', 'Ir25a')
corecptorExpGenes.df[2,1:2] <-  c('Orco', 'Ir76b')
corecptorExpGenes.df[3,1:2] <-  c('Orco', 'Ir8a')
corecptorExpGenes.df[4,1:2] <-  c('Ir25a', 'Ir76b')
corecptorExpGenes.df[5,1:2] <-  c('Ir25a', 'Ir8a')
corecptorExpGenes.df[6,1:2] <-  c('Ir76b', 'Ir8a')
for (i in 1:nrow(corecptorExpGenes.df)){
  rownames(corecptorExpGenes.df)[i] <- paste(corecptorExpGenes.df[i,1], corecptorExpGenes.df[i,2], sep= '_')
}
for (i in 1:nrow(corecptorExpGenes.df)) {
  corecptorExpGenes.df[i,3:6] <- c(
    length(
      receptorRNA.data[ 
        corecptorExpGenes.df[i,1], receptorRNA.data[corecptorExpGenes.df[i,1],] >= 1 & receptorRNA.data[corecptorExpGenes.df[i,2],] < 1
      ]),
    length(
      receptorRNA.data[ 
        corecptorExpGenes.df[i,1], receptorRNA.data[corecptorExpGenes.df[i,1],] < 1 & receptorRNA.data[corecptorExpGenes.df[i,2],] >= 1
      ]),
    length(
      receptorRNA.data[ 
        corecptorExpGenes.df[i,1], receptorRNA.data[corecptorExpGenes.df[i,1],] >= 1 & receptorRNA.data[corecptorExpGenes.df[i,2],] >= 1
      ]),
    length(
      receptorRNA.data[ 
        corecptorExpGenes.df[i,1], receptorRNA.data[corecptorExpGenes.df[i,1],] < 1 & receptorRNA.data[corecptorExpGenes.df[i,2],] < 1
      ])
  )
}
# corecptorExpGenes.df$sumCells <- apply(corecptorExpGenes.df, 1, function(input){
#   return(sum(as.numeric(input[3:6])))
# })

corecptorExpGenes.df <- as.data.frame(t(as.matrix(corecptorExpGenes.df)))

corecptorExpGenes.df %>% 
  pivot_longer(!criterion, names_to = 'ExpGroup', values_to = 'Ratio') %>%
  mutate(
    ExpGroup = factor(ExpGroup, levels =  c('Orco+/Ir25a+', 'Orco+', 'Ir25a+', 'Orco-/Ir25-')),
  ) %>%
  group_by(criterion) %>%
  ggplot(aes(criterion, Ratio, fill=ExpGroup)) + geom_bar(stat="identity", position = "fill")+
  scale_y_continuous(labels = scales::percent_format())






ratio.df <- data.frame(matrix(ncol = 4, nrow = 2))
colnames(ratio.df) <- c('Orco+', 'Ir25a+', 'Orco+/Ir25a+', 'Orco-/Ir25-')

ratio.df[1,] <- c(3152,528,2593,372)
ratio.df[2,] <- c(4627,571,588,859)
rownames(ratio.df) <- c('nor.exp_1', 'nor.exp_2')
ratio.df$criterion <- c('nor.exp_1', 'nor.exp_2')




nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp >=1 & neuron.batchCorrect@meta.data$Ir25a_Exp < 1, ])
nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp < 1 & neuron.batchCorrect@meta.data$Ir25a_Exp >= 1, ])
nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp >=1 & neuron.batchCorrect@meta.data$Ir25a_Exp >= 1, ])
nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp <1 & neuron.batchCorrect@meta.data$Ir25a_Exp < 1, ])

nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp >=2 & neuron.batchCorrect@meta.data$Ir25a_Exp < 2, ])
nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp < 2 & neuron.batchCorrect@meta.data$Ir25a_Exp >= 2, ])
nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp >=2 & neuron.batchCorrect@meta.data$Ir25a_Exp >= 2, ])
nrow(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_Exp <2 & neuron.batchCorrect@meta.data$Ir25a_Exp < 2, ])


ratio.df %>% 
  pivot_longer(!criterion, names_to = 'ExpGroup', values_to = 'Ratio') %>%
  mutate(
    ExpGroup = factor(ExpGroup, levels =  c('Orco+/Ir25a+', 'Orco+', 'Ir25a+', 'Orco-/Ir25-')),
  ) %>%
  group_by(criterion) %>%
  ggplot(aes(criterion, Ratio, fill=ExpGroup)) + geom_bar(stat="identity", position = "fill")+
  scale_y_continuous(labels = scales::percent_format())



cluster.seurat <- subset(neuron.batchCorrect, subset = seurat_clusters == 44)
cluster.seurat2 <- subset(neuron.batchCorrect, subset = seurat_clusters == 36)

# plotting coexpression of receptors
p0=FeatureScatter(object = cluster.seurat, feature1 = 'Orco', feature2 = 'Or6') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p1= FeatureScatter(object = cluster.seurat, feature1 = 'Orco', feature2 = 'Ir25a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p2=FeatureScatter(object = cluster.seurat, feature1 = 'Orco', feature2 = 'Ir76b') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p3=FeatureScatter(object = cluster.seurat2, feature1 = 'Orco', feature2 = 'Or6',  slot = 'data') +  
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p4=FeatureScatter(object = cluster.seurat2, feature1 = 'Orco', feature2 = 'Ir25a') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
p5=FeatureScatter(object = cluster.seurat2, feature1 = 'Orco', feature2 = 'Ir76b') + 
  geom_point(color='black') + theme(plot.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_receptors_normalizedExp.pdf')),
    width=15, height=10)
print((p0 + p1 + p2) /
        (p3 + p4 + p5))


##############################################################################

head(neuron.batchCorrect@meta.data, n=20)
tail(neuron.batchCorrect@meta.data, n=20)

# neuron.batchCorrect@meta.data[,18:21] <- NULL

dim(neuron.batchCorrect@meta.data)

length(receptorRNA.data['Orco',][receptorRNA.data["Orco",] >= 2] )
length(receptorRNA.data['Ir25a',][receptorRNA.data["Ir25a",] >= 1] )


coReceptorL

for (coreceptor in coReceptorL) {
  print(coreceptor)
  positiveL <- names( receptorRNA.data[coreceptor,][receptorRNA.data[coreceptor,] >= 1] )
  positiveL2 <- names( receptorRNA.data[coreceptor,][receptorRNA.data[coreceptor,] >= 2] )
  # print(positiveL)
  
  neuron.batchCorrect@meta.data[paste0(coreceptor, '_norExp1')] <- 
    as.numeric(
      lapply(rownames(neuron.batchCorrect@meta.data), function(x){
        # print(x)
        if (as.character(x) %in% positiveL) {
          return(1)
        } else {return(0)}
      })
    )
  
  if (coreceptor == 'Orco') {
    neuron.batchCorrect@meta.data[paste0(coreceptor, '_norExp2')] <- 
      as.numeric(
        lapply(rownames(neuron.batchCorrect@meta.data), function(x){
          # print(x)
          if (as.character(x) %in% positiveL2) {
            return(1)
          } else {return(0)}
        })
      )
  }
}


write_lines( row.names(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_norExp1 == 1,]), 
            file = paste0(receptorFolderPath, '/Orco_norExp1.list'))
write_lines( row.names(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Orco_norExp2 == 1,]), 
             file = paste0(receptorFolderPath, '/Orco_norExp2.list'))
write_lines( row.names(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Ir25a_norExp1 == 1,]), 
             file = paste0(receptorFolderPath, '/Ir25a_norExp1.list'))
write_lines( row.names(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Ir76b_norExp1 == 1,]), 
             file = paste0(receptorFolderPath, '/Ir76b_norExp1.list'))
write_lines( row.names(neuron.batchCorrect@meta.data[neuron.batchCorrect@meta.data$Ir8a_norExp1 == 1,]), 
             file = paste0(receptorFolderPath, '/Ir8a_norExp1.list'))


neuron.batchCorrect@meta.data$chemoreceptorAnno
neuronClusterIDs












##############################################################################
### Figure D: expressed # of Or/Ir/Gr at cluster level
### Histogram

# AverageExpression(neuron.batchCorrect, features = "Ir41k", slot = "data")
# AverageExpression(neuron.batchCorrect, features = "Ir41k", slot = "data")$RNA[53]
# log1p(AverageExpression(neuron.batchCorrect, features = "Ir41k", slot = "data")$RNA[53])


clusterReceptor.df <-data.frame(matrix(nrow=55, ncol = 4))
colnames(clusterReceptor.df) <- c('Cluster', 'Receptor', 'ReceptorType', 'ReceptorTypeN')
receptorType_L = vector()

# cluster level
for (clusterN in levels(receptor.dot.value$id)) {
  clusterRecptors <- receptor.dot.value %>% 
    filter(id == clusterN & pct.exp >=35 & avg.exp >= 1 & (! features.plot %in% coReceptorL)) %>%
    pull(features.plot) %>% as.character()
  # print(c('Cluster: ', clusterN))
  # print(clusterRecptors)
  # print(length(clusterRecptors))
  
  if (length(clusterRecptors) >= 1) {
    
    prefix_L = vector()
    for (receptor in clusterRecptors) {
      prefix=substring(receptor, 1, 1)
      # print(c(receptor, prefix))
      if (!prefix %in% prefix_L) {prefix_L=c(prefix_L, prefix)}
    }
    print(clusterRecptors)
    print(prefix_L)
    
    if (length(prefix_L) > 1) {
      receptorType='Mixed'
      receptorTypeN = 'Mixed'}
    else {
      if (prefix_L[1] == 'O') {
        receptorType='Ors only'
        if (length(clusterRecptors) == 1) {receptorTypeN = '1 Or'}
        else if (length(clusterRecptors) == 2) {receptorTypeN = '2 Or'}
        else {receptorTypeN = '3+ Or'}
      } else {
        receptorType='Irs only'
        if (length(clusterRecptors) == 1) {receptorTypeN = '1 Ir'}
        else if (length(clusterRecptors) == 2) {receptorTypeN = '2 Ir'}
        else {receptorTypeN = '3+ Ir'}
      }
      
    }
    
    clusterReceptor.df[as.numeric(clusterN) + 1, ] = 
      c(clusterN, paste(clusterRecptors, collapse = ','), receptorType, receptorTypeN)
    
  }
  else {
    clusterReceptor.df[as.numeric(clusterN) + 1, ] = c(clusterN, NA, NA, NA)
  }
  
}

clusterReceptor.df$ReceptorType <- 
  factor(clusterReceptor.df$ReceptorType, levels = c('Ors only', 'Irs only', 'Mixed'))
clusterReceptor.df$ReceptorTypeN <- 
  factor(clusterReceptor.df$ReceptorTypeN, levels = c('1 Or', '2 Or', '3+ Or', "1 Ir", "2 Ir", "3+ Ir", "Mixed"))

# clusterReceptor.df <- clusterReceptor.df %>%
#   group_by(ReceptorType) %>%
#   mutate(label_y = cumsum(Weight))

clusterReceptorSum.df <- clusterReceptor.df %>% 
  filter(! is.na(ReceptorType)) %>%
  group_by(ReceptorType, ReceptorTypeN) %>%
  summarise(n = n()) %>%
  mutate(label_y = cumsum(n),
         ReceptorTypeN = factor(ReceptorTypeN, levels = c('3+ Or','2 Or', '1 Or', "3+ Ir", "2 Ir", "1 Ir",  "Mixed")),
         ReceptorType = factor(ReceptorType, levels =  c('Ors only', 'Irs only', 'Mixed')))


clusterReceptorSum.df

p1 <- clusterReceptorSum.df %>%
  ggplot(aes(ReceptorType, y=n, fill=ReceptorTypeN)) + 
  geom_col(color="black") + 
  geom_text(aes(y = label_y, label=ReceptorTypeN), colour = "black", vjust = 1.5, size = 6) +
  theme(legend.position = "none", text = element_text(size = 14))+ ylab('Count')+
  scale_fill_manual(values=c('#deebf7','#9ecae1','#3182bd', '#e5f5e0','#a1d99b','#31a354', '#636363'))

pdf(file = file.path(picFolderPath, paste0('barGraph_receptorN_inEachCluster.pdf')), width=5, height=5)
print(p1)
dev.off()





### Illustration of cluster-specific genes
FeaturePlot(neuron.batchCorrect, features = c("Ir41k", "Or6"), blend = TRUE)
c("Ir41k", "Or132", 'Or6')

### Modify cluster names to receptor names

# Update name
new.cluster.ids <- c('Unknown','Unknown','Unknown','Or113','Unknown','Or87','Unknown',
                     'Unknown','Unknown','Or94','Or84','Unknown','Or36','Unknown','Or11',
                     'Ir41c','Unknown','Or63/Or64','Or38','Or23','Ir41a','Unknown','Or122',
                     'Ir75b','Or69','Ir41p','Or84/Or85','Ir75h','Ir75g','Ir93a','Or103',
                     'Unknown','Or103/Or115','Ir31a1','Or47/Or4','Or47/Or82/Or2','Ir41f',
                     'Or59','Ir87a1','Or79','Or130','Ir87a2','Or132','Or110','Or6','Ir64a',
                     'Or10','Ir75k','Or113/Or112','Or47/Or82','Or27/Or54','Or94/Gr77','Ir41k','Ir41l','Or44'
)

new.cluster.name.df <- read_csv(file = 'Analysis/integratedData_decontX/LogNormalize/7.exoressedReceptors/chemoreceptorClusterAnnotation.csv')
new.cluster.ids <- new.cluster.name.df$Receptor_S

length(new.cluster.ids)
length(levels(neuron.batchCorrect))

# cheange the cluster names
names(new.cluster.ids) <- new.cluster.name.df$ClusterN
new.cluster.ids
levels(neuron.batchCorrect)
neuron.batchCorrect <- RenameIdents(neuron.batchCorrect, new.cluster.ids)
neuron.batchCorrect@meta.data$chemoreceptorAnno <- Idents(neuron.batchCorrect)

print(
  (DimPlot(neuron.batchCorrect, reduction = "umap", label = TRUE, repel = TRUE) + theme(legend.position = "none")) +
        (DimPlot(neuron.batchCorrect, reduction = "tsne", label = TRUE, repel = TRUE) + theme(legend.position = "none"))
  )
neuron.batchCorrect@meta.data$seurat_clusters


pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_', 'chemoreceptorName.pdf')),
    width=15, height=7.5)
print(
  (DimPlot(neuron.batchCorrect, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters') + 
     theme(legend.position = "none", plot.title = element_blank())) +
    (DimPlot(neuron.batchCorrect, reduction = "umap", label = TRUE, repel = TRUE) + theme(legend.position = "none"))
)
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('tsne_', 'chemoreceptorName.pdf')),
    width=15, height=7.5)
print(
  (DimPlot(neuron.batchCorrect, reduction = "tsne", label = TRUE, repel = FALSE, group.by = 'seurat_clusters') + 
     theme(legend.position = "none", plot.title = element_blank())) +
    (DimPlot(neuron.batchCorrect, reduction = "tsne", label = TRUE, repel = FALSE) + theme(legend.position = "none"))
)
dev.off()




Idents(neuron.batchCorrect) <- neuron.batchCorrect@meta.data$seurat_clusters




### Generate feature plots for each receptor
eachReceptorFolder <-  file.path(picFolderPath, 'eachChemoreceptro')
dir.create(eachReceptorFolder, showWarnings = FALSE)

for (receptor in receptorExpSort_max35L) {
  pdf(file = file.path(eachReceptorFolder, 
                       paste0('FeaturePlot_', receptor, '.pdf')),
      width=5, height=5)
  print(FeaturePlot(neuron.batchCorrect, reduction = 'tsne', features = receptor, pt.size = 0.5))
  dev.off()
}



