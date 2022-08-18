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

rm(list = ls())   #remove all objects
theme_set(theme_cowplot())

setwd()
antenna.neuron <- readRDS("9_neurons_harmony_rmclusts_res4.rds")  #Uploaded on Zenodo as "SeuratObject2_Antenna_mergedBatches_Neurons.rds"

pipeline <- "harmony_rmclusts_res4"
pointsize = 1

DefaultAssay(antenna.neuron) <- "SCT"

dim(antenna.neuron@meta.data)[1]

######################################################################################
# tSNE
######################################################################################

# UMAP/tSNE
p1 <- DimPlot(antenna.neuron, reduction = "tsne", label = TRUE, pt.size = pointsize)
filename <- paste0("1a_tsne2.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p1)


p1 <- DimPlot(antenna.neuron, reduction = "umap", label = TRUE, pt.size = pointsize)
filename <- paste0("1b_umap.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 4000, height = 4000,p1)


######################################################################################
# GENERATE GENE LIST
######################################################################################

dat <- antenna.neuron
dat.assay <- "SCT"

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
receptorL <- append(coreceptors, GeneList$Gene) 

#receptorL <- GeneList$Gene

rm(ae, dat, GeneList, coreceptors)



######################################################################################
# DOTPLOT SET UP
######################################################################################

DefaultAssay(antenna.neuron) <- "SCT"

# Plot receptor 
receptor.dot <- DotPlot(antenna.neuron, features = receptorL) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)

#.....................................................................................
# RAW DOT PLOT
# ....................................................................................

filename <- paste0("1c1_receptordotplot_raw_",pipeline,".pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 16000, height = 4000,receptor.dot)

dotplot.features <- c("Or132","Or81","Ir41l","Ir41m","Or80","Or45","Or64","Or63",
                      "Or85","Or84","Or122","Or26","Or79","Or27","Or77","Or32",
                      "Or78","Or59","Ir87a1","Ir87a2","Gr77","Or94","Or100","Or97",
                      "Or110","Or16","Or38","Or121","Ir75g","Ir75e","Ir75f","Ir75i",
                      "Or103","Or114","Or115","Or72","Or52","Or71","Or6","Or82","Or3",
                      "Or47","Ir75h","Ir41p","Ir41a","Or11","Ir41f","Ir41c","Ir75d",
                      "Ir64a","Or36","Or88","Or87","Or4","Ir31a2","Ir31a1","Ir21a",
                      "Ir93a","Or23","Or69","Or112","Or113","Or105","Or111","Or104",
                      "Ir41b","Ir75a","Ir75b","Ir75c")


dotplot.features.rev <- append(rev(dotplot.features),c("Gr15","Or50","Or44","Or133","Or67","Or43"))
coreceptors <- c("LOC5575210","Orco","Ir25a","Ir76b","Ir8a")
dotplot.features.append <- append(coreceptors,dotplot.features.rev)

#dotplot.features.append <- append(coreceptors,dotplot.features.rev)

receptor.dot <- DotPlot(antenna.neuron, features = dotplot.features.append) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
receptor.dot

filename <- paste0("1c3_receptordotplot_",pipeline,".pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 10000, height = 4000,receptor.dot)

heatmap.genes5 <- c("LOC5575210","Orco"   ,"Ir25a"  ,"Ir76b"  ,
"Ir8a"      ,"Ir41k"  ,"Ir41p"  ,"Ir41a"  ,
"Ir41o"     ,"Ir41j"  ,"Or84"   ,"Or87"   ,
 "Gr77"     , "Ir64a" , "Or6"   , "Or47"  , 
 "Or132"    , "Or38"  , "Or25"  , "Or33"  , 
 "Or29"     , "Or27"  , "Or36"  , "Or88"  , 
 "Or97"     , "Or113" , "Or112" , "Or111" , 
 "Or104"    , "Or105" , "Or16"  , "Or41"  , 
 "Or42"     , "Or69"  , "Or63"  , "Or64"  , 
 "Or45"     , "Or122" , "Or121" , "Or91"  , 
 "Or130"    , "Or116" , "Or52"  , "Or72"  , 
 "Or10"     , "Or11"  , "Or94"  , "Or100" , 
 "Or103"    , "Or115" , "Or114" , "Or125" , 
 "Ir31a1"   , "Ir31a2", "Ir75k" , "Ir75b" , 
 "Ir75c"    , "Ir75a" , "Ir41c" , "Ir41b" , 
 "Ir41f"    , "Or4"   , "Or71"  , "Or82"  , 
 "Ir100a"   , "Or85"  , "Ir75h" , "Or58"  , 
 "Or23"     , "Or59"  , "Or31"  , "Or2"   , 
 "Ir93a"    , "Ir21a" , "Ir75g" , "Ir75f" , 
 "Ir75e"    , "Ir75i" , "Ir41e" , "Or3"   , 
 "Gr15"     , "Or50"  , "Ir87a1", "Ir87a2", 
 "Or110"    , "Or22"  , "Ir75d" , "Ir41l" , 
 "Or80"     , "Or81"  , "Ir41m" , "Or79"  , 
 "Or32"     , "Or78"  , "Or26"  , "Or44"  , 
 "Or133"    , "Or67"  , "Or43")

receptor.dot <- DotPlot(antenna.neuron, features = heatmap.genes5) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
receptor.dot

filename <- paste0("1c4_receptordotplot_heatmapgenes_",pipeline,".pdf")
ggsave(filename, path=outputFolder, limitsize = FALSE, units = "px", width = 10000, height = 4000,receptor.dot)



'''
my_levels <- clusterExpSortL_rmNoChemoL


receptor.dot <- DotPlot(dat, features = dotplot.features.append, assay = "SCT") + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
receptor.dot
filename <- paste0("1a3_neuronFigures_receptordotplot_",pipeline,".pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 10000, height = 4000,receptor.dot)
'''

dat <- antenna.neuron
levels(dat@active.ident)
my_levels <- c("32","19","37","10","41","12","24","35","33","22","4","30","17","31",
               "20","0","27","23","21","1","36","15","28","18","13","34","7","14",
               "16","9","5","29","38","39")
dat@active.ident <- factor(x = dat@active.ident, levels = my_levels)
levels(dat@active.ident)

receptor.dot <- DotPlot(dat, features = dotplot.features.append) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
receptor.dot


#.....................................................................................
# MORE SET UP
# ....................................................................................

receptor.dot.value <- receptor.dot$data %>% as.data.frame()

coReceptorL <- c('Orco', 'Ir25a', 'Ir76b', 'Ir8a','LOC5575210')

receptor.dot.value.noCoR <- receptor.dot.value %>% filter(! features.plot %in% coReceptorL)
receptor.dot.value.noCoR <- receptor.dot.value.noCoR %>% arrange(desc(avg.exp))

receptorExpSortL <- unique(receptor.dot.value.noCoR$features.plot)
receptorExpSortL <- c(as.character(receptorExpSortL), coReceptorL)
#head(receptorExpSortL)
#tail(receptorExpSortL)

clusterExpSortL <- unique(receptor.dot.value.noCoR$id) %>% as.character()
#clusterExpSortL

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
#head(mat.T)



######################################################################################
# SCATTER PLOTS (ORCO V IR25A)
######################################################################################

#.....................................................................................
# SET UP
# ....................................................................................

DefaultAssay(antenna.neuron) <- "SCT"

#yellow, green, teal, blue, purple
#scatterColors <- c('#fde725', '#5ec962', '#21918c', '#3b528b', '#440154')
#scatterColors <- c('#A06CB4', '#DF6C78', 'blue', '#fde725', '#5ec962', '#21918c', '#3b528b', '#440154')
scatterColors <- c('#A06CB4', '#DF6C78', '#911A2E', '#CD9139', '#B4B4B6', '#21918c', '#3b528b', '#440154')

#911A2E = Ir76b
#CD9139 = Ir8a
#B4B4B6 = Grey

### Plotting scatter plot: 1. single-cell level, 2. cluster level

### 1, single-cell level
## Orco vs. Ir25a

# head(antenna.neuron@assays$RNA@counts[c('Orco','Ir25a'),])
# orcoIl25a.mtx <- antenna.neuron@assays$RNA@counts[c('Orco','Ir25a'),]

antenna.neuron$Orco_UMIs <- antenna.neuron@assays$SCT@counts['Orco',]
antenna.neuron$Ir25a_UMIs <- antenna.neuron@assays$SCT@counts['Ir25a',]
antenna.neuron$Ir8a_UMIs <- antenna.neuron@assays$SCT@counts['Ir8a',]
antenna.neuron$Ir76b_UMIs <- antenna.neuron@assays$SCT@counts['Ir76b',]
antenna.neuron$LOC5575210_UMIs <- antenna.neuron@assays$SCT@counts['LOC5575210',]

antenna.neuron$Orco_Exp <- antenna.neuron@assays$SCT@data['Orco',]
antenna.neuron$Ir25a_Exp <- antenna.neuron@assays$SCT@data['Ir25a',]
antenna.neuron$Ir8a_Exp <- antenna.neuron@assays$SCT@data['Ir8a',]
antenna.neuron$Ir76b_Exp <- antenna.neuron@assays$SCT@data['Ir76b',]
antenna.neuron$LOC5575210_Exp <- antenna.neuron@assays$SCT@data['LOC5575210',]

antenna.neuron[[geneName]] <- antenna.neuron@assays$SCT@data[geneName,]

output <- "12a_neuronfigure_scatterplots"
dir.create(output)

#.....................................................................................
# ORCO V IR25A UMI
# ....................................................................................

## UMI: Orco vs. Ir25a
pmain <- antenna.neuron@meta.data %>%
  ggplot( aes(Orco_UMIs, Ir25a_UMIs) ) + 
  geom_point(size=0.3) + 
  scale_x_log10() + scale_y_log10()
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = antenna.neuron@meta.data, aes(x = Orco_UMIs), fill=scatterColors[1]) +scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = antenna.neuron@meta.data, aes(x = Ir25a_UMIs), fill=scatterColors[2]) + scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("2a1_neuronFigures_orco-v-ir25a_UMI_",pipeline,".pdf")
ggsave(filename, path=output, limitsize = FALSE, units = "px", width = 4000, height = 4000,p3)

#head(antenna.neuron@meta.data)

#.....................................................................................
# ORCO V IR25A NORMALIZED EXPRESSION (EXP>1)
# ....................................................................................

## normalized exp: Orco vs. Ir25a (Orco cutoff 1)
pmain <- antenna.neuron@meta.data %>%
  ggplot( aes(Orco_Exp, Ir25a_Exp) ) + 
  geom_point(size=0.3) +
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 1, color=scatterColors[3])
pmain <- pmain #+  geom_bin2d(bins = 20)

# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = antenna.neuron@meta.data, aes(x = Orco_Exp), fill=scatterColors[1]) +
  geom_vline(xintercept = 1, color=scatterColors[3])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = antenna.neuron@meta.data, aes(x = Ir25a_Exp), fill=scatterColors[2]) + 
  geom_vline(xintercept = 1, color=scatterColors[3])+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("2a2_neuronFigures_orco-v-ir25a_normexp1_",pipeline,".pdf")
ggsave(filename, path=output, limitsize = FALSE, units = "px", width = 4000, height = 4000,p3)

#.....................................................................................
# ORCO V IR25A NORMALIZED EXPRESSION (ORCO EXP>2)
# ....................................................................................

### normalized exp: Orco vs. Ir25a (Orco cutoff 2)
pmain <- antenna.neuron@meta.data %>%
  ggplot( aes(Orco_Exp, Ir25a_Exp) ) + 
  geom_point(size=0.3) +
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 2, color=scatterColors[3])
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = antenna.neuron@meta.data, aes(x = Orco_Exp), fill=scatterColors[1]) +
  geom_vline(xintercept = 2, color=scatterColors[3])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = antenna.neuron@meta.data, aes(x = Ir25a_Exp), fill=scatterColors[2]) + 
  geom_vline(xintercept = 1, color=scatterColors[3])+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("2a3_neuronFigures_orco-v-ir25a_normexp2_",pipeline,".pdf")
ggsave(filename, path=output, limitsize = FALSE, units = "px", width = 4000, height = 4000,p3)

#.....................................................................................
# ORCO V IR25A CLUSTER LEVEL
# ....................................................................................

plot_avgExp.clusters.df <- 
  receptor.dot.value %>%
  filter(features.plot %in% c(coReceptorL)) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 



### 2. Orco vs. Ir25a in cluster level
plot_avgExp.clusters.df <- 
  receptor.dot.value %>%
  filter(features.plot %in% c('Orco', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Orco, Ir25a, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Orco') + ylab('Average expression of Ir25a')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Orco), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Ir25a), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("2a4_neuronFigures_orco-v-ir25a_clusters_",pipeline,".pdf")
ggsave(filename, path=output, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)


### 2. Ir76b vs. Ir25a in cluster level
plot_avgExp.clusters.df <- 
  receptor.dot.value %>%
  filter(features.plot %in% c('Ir76b', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Ir76b, Ir25a, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Ir76b') + ylab('Average expression of Ir25a')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Ir76b), fill=scatterColors[3])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Ir25a), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p7 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p8 <- insert_yaxis_grob(p7, ydens, grid::unit(.2, "null"), position = "right")
p9 <- ggdraw(p8)
p9

filename <- paste0("2a6_neuronFigures_ir25a-v-ir8a_clusters_",pipeline,".pdf")
ggsave(filename,path=output, limitsize = FALSE, units = "px", width = 2000, height = 2000,p6)

### 2. Ir8a vs. Ir25a in cluster level
plot_avgExp.clusters.df <- 
  receptor.dot.value %>%
  filter(features.plot %in% c('Ir8a', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Ir8a, Ir25a, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Ir8a') + ylab('Average expression of Ir25a')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Ir8a), fill=scatterColors[4])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Ir25a), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p5 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p5 <- insert_yaxis_grob(p5, ydens, grid::unit(.2, "null"), position = "right")
p6 <- ggdraw(p5)
p6

filename <- paste0("2a5_neuronFigures_ir25a-v-ir8a_clusters_",pipeline,".pdf")
ggsave(filename, path=output,limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)


### 2. LOC5575210 vs. Ir25a in cluster level
plot_avgExp.clusters.df <- 
  receptor.dot.value %>%
  filter(features.plot %in% c('LOC5575210', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 


pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(LOC5575210, Ir25a, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of LOC5575210') + ylab('Average expression of Ir25a')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC5575210), fill=scatterColors[5])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Ir25a), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p10 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p11 <- insert_yaxis_grob(p10, ydens, grid::unit(.2, "null"), position = "right")
p12 <- ggdraw(p11)
p12

filename <- paste0("2a7_neuronFigures_ir25a-v-LOC5575210_clusters_",pipeline,".pdf")
ggsave(filename, path=output,limitsize = FALSE, units = "px", width = 2000, height = 2000,p9)


g = arrangeGrob(p3,p9,p6,p12, ncol = 4)
filename <- paste0("2a8_neuronFigures_coreceptorscatter_clusters_",pipeline,".pdf")
ggsave(filename, path=output, limitsize = FALSE, units = "px", width = 8000, height = 2000, g)




######################################################################################
# VENN DIAGRAM RAW COUNTS
######################################################################################

#.....................................................................................
# EXPRESSION LIST
# ....................................................................................

receptorSCT = antenna.neuron@assays$SCT@counts[receptorExpSortL,]
receptorSCT.data = antenna.neuron@assays$SCT@data[receptorExpSortL,]

for (coreceptor in coReceptorL) {
  print(coreceptor)
  positiveL <- names( receptorSCT.data[coreceptor,][receptorSCT.data[coreceptor,] >= 1] )
  positiveL2 <- names( receptorSCT.data[coreceptor,][receptorSCT.data[coreceptor,] >= 2] )
  # print(positiveL)
  
  antenna.neuron@meta.data[paste0(coreceptor, '_norExp1')] <- 
    as.numeric(
      lapply(rownames(antenna.neuron@meta.data), function(x){
        # print(x)
        if (as.character(x) %in% positiveL) {
          return(1)
        } else {return(0)}
      })
    )
  
  if (coreceptor == 'Orco') {
    antenna.neuron@meta.data[paste0(coreceptor, '_norExp2')] <- 
      as.numeric(
        lapply(rownames(antenna.neuron@meta.data), function(x){
          # print(x)
          if (as.character(x) %in% positiveL2) {
            return(1)
          } else {return(0)}
        })
      )
  }
}


write_lines( 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Orco_norExp1 == 1,]), 
             file = paste0(outputFolder, '/Orco_norExp1.list'))
write_lines( 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Orco_norExp2 == 1,]), 
             file = paste0(outputFolder, '/Orco_norExp2.list'))
write_lines( 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir25a_norExp1 == 1,]), 
             file = paste0(outputFolder, '/Ir25a_norExp1.list'))
write_lines( 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir76b_norExp1 == 1,]), 
             file = paste0(outputFolder, '/Ir76b_norExp1.list'))
write_lines( 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir8a_norExp1 == 1,]), 
             file = paste0(outputFolder, '/Ir8a_norExp1.list'))




library(VennDiagram)

dataset1 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Orco_norExp1 == 1,])
dataset2 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir25a_norExp1 == 1,])
dataset3 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir76b_norExp1 == 1,])
dataset4 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir8a_norExp1 == 1,])
  
name1 <- 'Orco , norm.exp > 1'
name2 <- 'Ir25a , norm.exp > 1'
name3 <- 'Ir76b , norm.exp > 1'
name4 <- 'Ir8a , norm.exp > 1'

venn.diagram(
  x = list(dataset1, dataset2, 
           dataset3, dataset4),
  category.names = c(name1, name2, name3, name4),
  filename = paste0('coreceptors_normexp1.png'),
  output=TRUE
)

dataset1 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Orco_norExp2 == 1,])
dataset2 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir25a_norExp1 == 1,])
dataset3 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir76b_norExp1 == 1,])
dataset4 <- 
  row.names(antenna.neuron@meta.data[antenna.neuron@meta.data$Ir8a_norExp1 == 1,])

name1 <- 'Orco , norm.exp > 2'
name2 <- 'Ir25a , norm.exp > 1'
name3 <- 'Ir76b , norm.exp > 1'
name4 <- 'Ir8a , norm.exp > 1'

venn.diagram(
  x = list(dataset1, dataset2, 
           dataset3, dataset4),
  category.names = c(name1, name2, name3, name4),
  filename = paste0('coreceptors_normexp1_orco2.png'),
  output=TRUE
)




######################################################################################
# CHORD PLOT (SCT VALS)
######################################################################################

#.....................................................................................
# SET UP
# ....................................................................................
library(chorddiag)
library(htmltools)

#rm(list = ls())   #remove all objects
theme_set(theme_cowplot())
#dev.off()
  
setwd("E:/GDrive Mirror/Lab/Lab Projects/2. scRNAseq/MOC/2022-03-14_pipelinerevisions/1_FINALPIPELINE")
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")

pipeline <- "harmony_rmclusts_res4"
filename <- '12_neuronfigureplots'
outputFolder = paste0(getwd(), '/', filename)
setwd(outputFolder)
#getwd()

DefaultAssay(antenna.neuron) <- "SCT"


filename <- paste0("1q_neuronFigures_chordplot_",pipeline)
dir.create(toString(filename))

setwd(paste(outputFolder,filename,sep="/"))
getwd()
dat.assay <- 'SCT'

#.....................................................................................
# GENERATE GENE LIST
# ....................................................................................

dat <- antenna.neuron
dat.assay <- "SCT"

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

#coreceptors <- c('LOC5575210','Orco','Ir25a','Ir76b','Ir8a')
#receptorL <- append(coreceptors, GeneList$Gene) 

receptorL <- GeneList$Gene

rm(ae, dat, GeneList, coreceptors)



##########################################################
### chord plot

receptorUMI.mtx <- antenna.neuron[ rownames(antenna.neuron[[dat.assay]]@counts) %in% receptorL, ]@assays$SCT@counts
dim(receptorUMI.mtx)
receptorNor.mtx <- antenna.neuron[ rownames(antenna.neuron[[dat.assay]]@data) %in% receptorL, ]@assays$SCT@data
dim(receptorNor.mtx)

# Filter non-expressed genes
receptorUMI.mtx.expCol=receptorUMI.mtx[rowSums(receptorUMI.mtx) != 0, ]
receptorNor.mtx.expCol=receptorNor.mtx[rowSums(receptorNor.mtx) != 0, ]
print(c(dim(receptorUMI.mtx.expCol), dim(receptorNor.mtx.expCol)))
# Filter non-expressed cells
receptorUMI.mtx.expCol=receptorUMI.mtx.expCol[, colSums(receptorUMI.mtx) != 0 ]
receptorNor.mtx.expCol=receptorNor.mtx.expCol[, colSums(receptorNor.mtx) != 0, ]
print(c(dim(receptorUMI.mtx.expCol), dim(receptorNor.mtx.expCol)))

print("Writing receptorUMI.csv ....")
write.csv(as.data.frame(receptorUMI.mtx.expCol), file = 'receptorUMI.csv')
print("Writing receptorNorExp.csv ....")
write.csv(as.data.frame(receptorNor.mtx.expCol), file = 'receptorNorExp.csv')

receptorN <- nrow(receptorUMI.mtx.expCol)

# receptorUMI.mtx.expCol %>% summarise_if(is.numeric, max)
# apply(receptorUMI.mtx.expCol, 1, max)
# length(apply(receptorUMI.mtx.expCol, 1, max))



### Use different UMI thresholds to check the coexpression

# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
colnames(coexp.df) <- rownames(receptorUMI.mtx.expCol)
rownames(coexp.df) <- rownames(receptorUMI.mtx.expCol)

coexp.df2 <- coexp.df
coexp.df3 <- coexp.df
coexp.df5 <- coexp.df

# coexp.df3 <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
# colnames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)
# rownames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)


# normalized expression

print("Running loop for calculating UMI thresholds ....")

for (i in 1:receptorN) {
  for (j in 1:receptorN) {
    # print(c('i = ', i))
    # print(c('j = ', j))
    
    if (i==j) {
      coexp.df[i, j] <- 0
      coexp.df2[i, j] <- 0
      coexp.df3[i, j] <- 0
      coexp.df5[i, j] <- 0
    }
    
    else {
      lapply(1:length(ncol(receptorUMI.mtx.expCol)), function(x){
        # print(c('i = ', i))
        # print(c('j = ', j))
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 1 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        coexp.df[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 2 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 2])
        coexp.df2[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 3 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 3])
        coexp.df3[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 5 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 5])
        coexp.df5[i,j]<<-jExpressedN
      })
    }
  }
}

print("Writing coexp_UMI1.csv ....")
write.csv(coexp.df, file = 'coexp_UMI1.csv')
print("Writing coexp_UMI2.csv ....")
write.csv(coexp.df2, file = 'coexp_UMI2.csv')
print("Writing coexp_UMI3.csv ....")
write.csv(coexp.df3, file = 'coexp_UMI3.csv')
print("Writing coexp_UMI5.csv ....")
write.csv(coexp.df5, file = 'coexp_UMI5.csv')


# names(hex_codes1) <- colnames(coexp.df)


hex_codes1 <- hue_pal()(receptorN)
print(chorddiag(as.matrix(coexp.df), groupColors = hex_codes1, groupnamePadding = 10, showTicks = FALSE))
chorddiag(as.matrix(coexp.df2), groupColors = hex_codes1, groupnamePadding = 20)
chorddiag(as.matrix(coexp.df3), groupColors = hex_codes1, groupnamePadding = 30, showTooltips = FALSE, showZeroTooltips = FALSE)


#.....................................................................................
# UMI CHORD PLOT
# ....................................................................................


### Use different cell# cutoff for plotting chord plot
### Remove coexpression <= 10 or 20

### coexp.df: UMI 1
# groupnamePadding: shift of group name
coexp.df.filter10 <- coexp.df
coexp.df.filter10[coexp.df <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter10)) {
  if (sum(coexp.df.filter10[,i]) == 0) {
    print(colnames(coexp.df.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter10 <- coexp.df.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df.filter10))
p <- chorddiag(as.matrix(coexp.df.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 40, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_UMI1.html", libdir="coexp_UMI1")
rm(p)

### coexp.df2: UMI 2
# groupnamePadding: shift of group name
coexp.df2.filter10 <- coexp.df2
coexp.df2.filter10[coexp.df2 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df2.filter10)) {
  if (sum(coexp.df2.filter10[,i]) == 0) {
    print(colnames(coexp.df2.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df2.filter10 <- coexp.df2.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df2.filter10))
p <- chorddiag(as.matrix(coexp.df2.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)


save_html(p, file = "coexp_UMI2.html", libdir="coexp_UMI2")
rm(p)

### coexp.df3: UMI 3
# groupnamePadding: shift of group name
coexp.df3.filter10 <- coexp.df3
coexp.df3.filter10[coexp.df3 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df3.filter10)) {
  if (sum(coexp.df3.filter10[,i]) == 0) {
    print(colnames(coexp.df3.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df3.filter10 <- coexp.df3.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df3.filter10))
p <- chorddiag(as.matrix(coexp.df3.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_UMI3.html", libdir="coexp_UMI3")
rm(p)

### coexp.df5: UMI 5
# groupnamePadding: shift of group name
coexp.df5.filter10 <- coexp.df5
coexp.df5.filter10[coexp.df5 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df5.filter10)) {
  if (sum(coexp.df5.filter10[,i]) == 0) {
    print(colnames(coexp.df5.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df5.filter10 <- coexp.df5.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df5.filter10))
p <- chorddiag(as.matrix(coexp.df5.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 20, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_UMI5.html", libdir="coexp_UMI5")
rm(p)

#.....................................................................................
# NORMALIZED EXP CHORD PLOT
# ....................................................................................


###
### Use different normalization expression thresholds to check the coexpression


# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
colnames(coexp.df) <- rownames(receptorNor.mtx.expCol)
rownames(coexp.df) <- rownames(receptorNor.mtx.expCol)

coexp.df05 <- coexp.df
coexp.df2 <- coexp.df
coexp.df3 <- coexp.df


# coexp.df3 <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
# colnames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)
# rownames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)


# normalized expression
for (i in 1:receptorN) {
  for (j in 1:receptorN) {
    # print(c('i = ', i))
    # print(c('j = ', j))
    
    if (i==j) {
      coexp.df[i, j] <- 0
      coexp.df2[i, j] <- 0
      coexp.df3[i, j] <- 0
      coexp.df05[i, j] <- 0
    }
    
    else {
      lapply(1:length(ncol(receptorNor.mtx.expCol)), function(x){
        # print(c('i = ', i))
        # print(c('j = ', j))
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 1 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        coexp.df[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 2 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 2])
        coexp.df2[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 3 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 3])
        coexp.df3[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 0.5 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 0.5])
        coexp.df05[i,j]<<-jExpressedN
      })
    }
  }
}

write.csv(coexp.df, file = 'coexp_norExp1.csv')
write.csv(coexp.df2, file = 'coexp_norExp2.csv')
write.csv(coexp.df3, file = 'coexp_norExp3.csv')
write.csv(coexp.df05, file = 'coexp_norExp0.5.csv')



### Use different cell# cutoff for plotting chord plot
### Remove coexpression <= 10 or 20

### coexp.df: norExp1
# groupnamePadding: shift of group name
coexp.df.filter10 <- coexp.df
coexp.df.filter10[coexp.df <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter10)) {
  if (sum(coexp.df.filter10[,i]) == 0) {
    print(colnames(coexp.df.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter10 <- coexp.df.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df.filter10))
p <- chorddiag(as.matrix(coexp.df.filter10), groupColors = hex_codes1, 
               groupnamePadding = 40, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp1.html", libdir="coexp_norExp1")
rm(p)


### coexp.df2: norExp 2
# groupnamePadding: shift of group name
coexp.df2.filter10 <- coexp.df2
coexp.df2.filter10[coexp.df2 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df2.filter10)) {
  if (sum(coexp.df2.filter10[,i]) == 0) {
    print(colnames(coexp.df2.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df2.filter10 <- coexp.df2.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df2.filter10))
p <- chorddiag(as.matrix(coexp.df2.filter10), groupColors = hex_codes1, 
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp2.html", libdir="coexp_norExp2")
rm(p)


### coexp.df3: norExp 3
# groupnamePadding: shift of group name
coexp.df3.filter10 <- coexp.df3
coexp.df3.filter10[coexp.df3 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df3.filter10)) {
  if (sum(coexp.df3.filter10[,i]) == 0) {
    print(colnames(coexp.df3.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df3.filter10 <- coexp.df3.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df3.filter10))
p <- chorddiag(as.matrix(coexp.df3.filter10), groupColors = hex_codes1, 
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp3.html", libdir="coexp_norExp3")
rm(p)


### coexp.df05: norExp 0.5
# groupnamePadding: shift of group name
coexp.df05.filter10 <- coexp.df05
coexp.df05.filter10[coexp.df05 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df05.filter10)) {
  if (sum(coexp.df05.filter10[,i]) == 0) {
    print(colnames(coexp.df05.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df05.filter10 <- coexp.df05.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df05.filter10))
p <- chorddiag(as.matrix(coexp.df05.filter10), groupColors = hex_codes1, 
               groupnamePadding = 40, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp0.5.html", libdir="coexp_norExp0.5")
rm(p)


######################################################################################
# CHORD PLOT (RNA VALS)
######################################################################################

#.....................................................................................
# SET UP
# ....................................................................................
library(chorddiag)
library(htmltools)



rm(list = ls())   #remove all objects
theme_set(theme_cowplot())
#dev.off()

setwd("E:/GDrive Mirror/Lab/Lab Projects/2. scRNAseq/MOC/2022-03-14_pipelinerevisions/1_FINALPIPELINE")
antenna.neuron <- readRDS("9_neurons_harmony_rm0_res4.rds")

pipeline <- "harmony_rmclusts_res4"
filename <- '12_neuronfigureplots'
outputFolder = paste0(getwd(), '/', filename)
setwd(outputFolder)
#getwd()

DefaultAssay(antenna.neuron) <- "RNA"


filename <- paste0("1r_neuronFigures_chordplot_RNA_",pipeline)
dir.create(toString(filename))

setwd(paste(outputFolder,filename,sep="/"))
getwd()
dat.assay <- 'RNA'

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

#coreceptors <- c('LOC5575210','Orco','Ir25a','Ir76b','Ir8a')
#receptorL <- append(coreceptors, GeneList$Gene) 

receptorL <- GeneList$Gene

rm(ae, dat, GeneList, coreceptors)



##########################################################
### chord plot

receptorUMI.mtx <- antenna.neuron[ rownames(antenna.neuron[[dat.assay]]@counts) %in% receptorL, ]@assays[[dat.assay]]@counts
dim(receptorUMI.mtx)
receptorNor.mtx <- antenna.neuron[ rownames(antenna.neuron[[dat.assay]]@data) %in% receptorL, ]@assays[[dat.assay]]@data
dim(receptorNor.mtx)

# Filter non-expressed genes
receptorUMI.mtx.expCol=receptorUMI.mtx[rowSums(receptorUMI.mtx) != 0, ]
receptorNor.mtx.expCol=receptorNor.mtx[rowSums(receptorNor.mtx) != 0, ]
print(c(dim(receptorUMI.mtx.expCol), dim(receptorNor.mtx.expCol)))
# Filter non-expressed cells
receptorUMI.mtx.expCol=receptorUMI.mtx.expCol[, colSums(receptorUMI.mtx) != 0 ]
receptorNor.mtx.expCol=receptorNor.mtx.expCol[, colSums(receptorNor.mtx) != 0, ]
print(c(dim(receptorUMI.mtx.expCol), dim(receptorNor.mtx.expCol)))

print("Writing receptorUMI.csv ....")
write.csv(as.data.frame(receptorUMI.mtx.expCol), file = 'receptorUMI.csv')
print("Writing receptorNorExp.csv ....")
write.csv(as.data.frame(receptorNor.mtx.expCol), file = 'receptorNorExp.csv')

receptorN <- nrow(receptorUMI.mtx.expCol)

# receptorUMI.mtx.expCol %>% summarise_if(is.numeric, max)
# apply(receptorUMI.mtx.expCol, 1, max)
# length(apply(receptorUMI.mtx.expCol, 1, max))



### Use different UMI thresholds to check the coexpression

# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
colnames(coexp.df) <- rownames(receptorUMI.mtx.expCol)
rownames(coexp.df) <- rownames(receptorUMI.mtx.expCol)

coexp.df2 <- coexp.df
coexp.df3 <- coexp.df
coexp.df5 <- coexp.df

# coexp.df3 <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
# colnames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)
# rownames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)


# normalized expression

print("Running loop for calculating UMI thresholds ....")

for (i in 1:receptorN) {
  for (j in 1:receptorN) {
    # print(c('i = ', i))
    # print(c('j = ', j))
    
    if (i==j) {
      coexp.df[i, j] <- 0
      coexp.df2[i, j] <- 0
      coexp.df3[i, j] <- 0
      coexp.df5[i, j] <- 0
    }
    
    else {
      lapply(1:length(ncol(receptorUMI.mtx.expCol)), function(x){
        # print(c('i = ', i))
        # print(c('j = ', j))
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 1 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        coexp.df[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 2 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 2])
        coexp.df2[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 3 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 3])
        coexp.df3[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorUMI.mtx.expCol[, receptorUMI.mtx.expCol[i,] >= 5 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 5])
        coexp.df5[i,j]<<-jExpressedN
      })
    }
  }
}

print("Writing coexp_UMI1.csv ....")
write.csv(coexp.df, file = 'coexp_UMI1.csv')
print("Writing coexp_UMI2.csv ....")
write.csv(coexp.df2, file = 'coexp_UMI2.csv')
print("Writing coexp_UMI3.csv ....")
write.csv(coexp.df3, file = 'coexp_UMI3.csv')
print("Writing coexp_UMI5.csv ....")
write.csv(coexp.df5, file = 'coexp_UMI5.csv')


# names(hex_codes1) <- colnames(coexp.df)


hex_codes1 <- hue_pal()(receptorN)
print(chorddiag(as.matrix(coexp.df), groupColors = hex_codes1, groupnamePadding = 10, showTicks = FALSE))
chorddiag(as.matrix(coexp.df2), groupColors = hex_codes1, groupnamePadding = 20)
chorddiag(as.matrix(coexp.df3), groupColors = hex_codes1, groupnamePadding = 30, showTooltips = FALSE, showZeroTooltips = FALSE)


#.....................................................................................
# UMI CHORD PLOT
# ....................................................................................


### Use different cell# cutoff for plotting chord plot
### Remove coexpression <= 10 or 20

### coexp.df: UMI 1
# groupnamePadding: shift of group name
coexp.df.filter10 <- coexp.df
coexp.df.filter10[coexp.df <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter10)) {
  if (sum(coexp.df.filter10[,i]) == 0) {
    print(colnames(coexp.df.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter10 <- coexp.df.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df.filter10))
p <- chorddiag(as.matrix(coexp.df.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 40, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_UMI1.html", libdir="coexp_UMI1")
rm(p)

### coexp.df2: UMI 2
# groupnamePadding: shift of group name
coexp.df2.filter10 <- coexp.df2
coexp.df2.filter10[coexp.df2 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df2.filter10)) {
  if (sum(coexp.df2.filter10[,i]) == 0) {
    print(colnames(coexp.df2.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df2.filter10 <- coexp.df2.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df2.filter10))
p <- chorddiag(as.matrix(coexp.df2.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)


save_html(p, file = "coexp_UMI2.html", libdir="coexp_UMI2")
rm(p)

### coexp.df3: UMI 3
# groupnamePadding: shift of group name
coexp.df3.filter10 <- coexp.df3
coexp.df3.filter10[coexp.df3 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df3.filter10)) {
  if (sum(coexp.df3.filter10[,i]) == 0) {
    print(colnames(coexp.df3.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df3.filter10 <- coexp.df3.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df3.filter10))
p <- chorddiag(as.matrix(coexp.df3.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_UMI3.html", libdir="coexp_UMI3")
rm(p)

### coexp.df5: UMI 5
# groupnamePadding: shift of group name
coexp.df5.filter10 <- coexp.df5
coexp.df5.filter10[coexp.df5 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df5.filter10)) {
  if (sum(coexp.df5.filter10[,i]) == 0) {
    print(colnames(coexp.df5.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df5.filter10 <- coexp.df5.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df5.filter10))
p <- chorddiag(as.matrix(coexp.df5.filter10), groupColors = hex_codes1, width = 1000, height = 1000,
               groupnamePadding = 20, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_UMI5.html", libdir="coexp_UMI5")
rm(p)

#.....................................................................................
# NORMALIZED EXP CHORD PLOT
# ....................................................................................


###
### Use different normalization expression thresholds to check the coexpression


# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
colnames(coexp.df) <- rownames(receptorNor.mtx.expCol)
rownames(coexp.df) <- rownames(receptorNor.mtx.expCol)

coexp.df05 <- coexp.df
coexp.df2 <- coexp.df
coexp.df3 <- coexp.df


# coexp.df3 <- data.frame(matrix(ncol = receptorN, nrow = receptorN))
# colnames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)
# rownames(coexp.df3) <- rownames(receptorUMI.mtx.expCol)


# normalized expression
for (i in 1:receptorN) {
  for (j in 1:receptorN) {
    # print(c('i = ', i))
    # print(c('j = ', j))
    
    if (i==j) {
      coexp.df[i, j] <- 0
      coexp.df2[i, j] <- 0
      coexp.df3[i, j] <- 0
      coexp.df05[i, j] <- 0
    }
    
    else {
      lapply(1:length(ncol(receptorNor.mtx.expCol)), function(x){
        # print(c('i = ', i))
        # print(c('j = ', j))
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 1 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        coexp.df[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 2 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 2])
        coexp.df2[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 3 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 3])
        coexp.df3[i,j]<<-jExpressedN
        
        iExpressedCells.mx = as.matrix(receptorNor.mtx.expCol[, receptorNor.mtx.expCol[i,] >= 0.5 ])
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 0.5])
        coexp.df05[i,j]<<-jExpressedN
      })
    }
  }
}

write.csv(coexp.df, file = 'coexp_norExp1.csv')
write.csv(coexp.df2, file = 'coexp_norExp2.csv')
write.csv(coexp.df3, file = 'coexp_norExp3.csv')
write.csv(coexp.df05, file = 'coexp_norExp0.5.csv')



### Use different cell# cutoff for plotting chord plot
### Remove coexpression <= 10 or 20

### coexp.df: norExp1
# groupnamePadding: shift of group name
coexp.df.filter10 <- coexp.df
coexp.df.filter10[coexp.df <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter10)) {
  if (sum(coexp.df.filter10[,i]) == 0) {
    print(colnames(coexp.df.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter10 <- coexp.df.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df.filter10))
p <- chorddiag(as.matrix(coexp.df.filter10), groupColors = hex_codes1, 
               groupnamePadding = 40, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp1.html", libdir="coexp_norExp1")
rm(p)


### coexp.df2: norExp 2
# groupnamePadding: shift of group name
coexp.df2.filter10 <- coexp.df2
coexp.df2.filter10[coexp.df2 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df2.filter10)) {
  if (sum(coexp.df2.filter10[,i]) == 0) {
    print(colnames(coexp.df2.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df2.filter10 <- coexp.df2.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df2.filter10))
p <- chorddiag(as.matrix(coexp.df2.filter10), groupColors = hex_codes1, 
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp2.html", libdir="coexp_norExp2")
rm(p)


### coexp.df3: norExp 3
# groupnamePadding: shift of group name
coexp.df3.filter10 <- coexp.df3
coexp.df3.filter10[coexp.df3 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df3.filter10)) {
  if (sum(coexp.df3.filter10[,i]) == 0) {
    print(colnames(coexp.df3.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df3.filter10 <- coexp.df3.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df3.filter10))
p <- chorddiag(as.matrix(coexp.df3.filter10), groupColors = hex_codes1, 
               groupnamePadding = 30, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp3.html", libdir="coexp_norExp3")
rm(p)


### coexp.df05: norExp 0.5
# groupnamePadding: shift of group name
coexp.df05.filter10 <- coexp.df05
coexp.df05.filter10[coexp.df05 <= 10] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df05.filter10)) {
  if (sum(coexp.df05.filter10[,i]) == 0) {
    print(colnames(coexp.df05.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df05.filter10 <- coexp.df05.filter10[-zeroGenes_V, -zeroGenes_V]
hex_codes1 <- hue_pal()(ncol(coexp.df05.filter10))
p <- chorddiag(as.matrix(coexp.df05.filter10), groupColors = hex_codes1, 
               groupnamePadding = 40, showTicks = TRUE, groupPadding=1, showZeroTooltips=FALSE)

save_html(p, file = "coexp_norExp0.5.html", libdir="coexp_norExp0.5")
rm(p)

}
