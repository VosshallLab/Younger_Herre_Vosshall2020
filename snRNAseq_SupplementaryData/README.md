# Table of Contents
Updated August 18, 2022 - for questions, contact ogoldman@rockefeller.edu  

## Repositories
All supplementary snRNA-seq analysis is in the following repository:
https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData

Antenna snRNAseq R data files, scripts, and expression matricies are on Zenodo:
https://doi.org/10.5281/zenodo.5818542

Maxillary Palp snRNAseq R data files, scripts, and expression matricies are on Zenodo:
https://doi.org/10.5281/zenodo.5818951

For raw sequencing files for both antenna and maxillary palp [NCBI BioProject: PRJNA794050](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=794050)

## *Antenna: List of Data Items*

| LOCATION 	| ITEM DESCRIPTION 	| DATA FILE USED (download [here](https://doi.org/10.5281/zenodo.5818542)) 	| SCRIPT USED TO GENERATE 	|
|---	|---	|---	|---	|
| Figure 5B 	| All Cells Marker Heatmap 	| SeuratObject1_Antenna_mergedBatches_AllCells.rds 	| Antenna_Figures2_AllCells.R 	|
| Figure 5C 	| Neurons tSNE 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures1_Neurons.R 	|
| Figure 5D 	| Chemosensor Dot plot 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures1_Neurons.R 	|
| Figure 5E 	| Chord Plot (SCT values)	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures1_Neurons.R 	|
| Figure 5F 	| _Or82_+ cells heatmap 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| Figure S5A 	| Orco Feature Plot (DecontX illustration) 	| cellranger_and_decontX_outputs 	| Antenna_Figures2_AllCells.R 	|
| Figure S5B 	| Batch Correction Feature Plot 	| SeuratObject1_Antenna_mergedBatches_AllCells.rds 	| Antenna_Pipeline.R 	|
| Figure S5C-E 	| Genes, UMIs, Mitochondrial Genes Violin Plots 	| SeuratObject1_Antenna_mergedBatches_AllCells.rds 	| Antenna_Figures2_AllCells.R 	|
| Figure S5F 	| Chemosensor Pseudobulk Scatter Plot	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures4_NeuronPseudobulk.R 	|
| Figure S5G 	| All Cells tSNE 	| SeuratObject1_Antenna_mergedBatches_AllCells.rds 	| Antenna_Pipeline.R 	|
| Figure S5H 	| Cell Marker Feature Plot 	| SeuratObject1_Antenna_mergedBatches_AllCells.rds 	| Antenna_Pipeline.R 	|
| Figure S5I 	| Neural Marker Dot Plot 	| SeuratObject1_Antenna_mergedBatches_AllCells.rds 	| Antenna_Pipeline.R 	|
| Figure S5J 	| Chord Plot (RNA values)	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures1_Neurons.R 	|
| Figure S5K 	| Co-receptor Venn Diagram | SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures1_Neurons.R 	|
| Figure S5L 	| Cluster 3 heatmap 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| Figure S5M 	| Cluster 0 heatmap 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| Figure S5N 	| _Or36_+ cells heatmap 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| Figure S5O 	| _Or47_+ cells heatmap 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| Figure S5P 	| _Or84_+ cells heatmap 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| Figure S6A 	| All neurons heatmap 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| Figure S6A 	| Co-receptor violin plots 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds 	| Antenna_Figures3_NeuronExampleHeatmaps.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_Figures4_NeuronPseudobulk_Output)  	| Chemosensor PseudoBulk Table	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| Antenna_Figures4_NeuronPseudobulk.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_SuppData1_Heatmaps_ViolinPlots/10b_byCluster_heatmaps_avggenelist_reclustered_res4) 	| Heatmaps (filtered by cluster) 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| Antenna_SuppData1_Heatmaps.R 	|
Antenna_SuppData1_Heatmaps.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_SuppData1_Heatmaps_ViolinPlots/10c_byCluster_violinplot_topgenes_reclustered_res4)	 	| Violin Plots (on chemosensor-filtered populations) 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| Antenna_SuppData1_Heatmaps.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_SuppData1_Heatmaps_ViolinPlots/10d_byCluster_batchcheck_dimensionreduction) 	 	| tSNE batch illustration (each cluster) 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| Antenna_SuppData1_Heatmaps.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_SuppData1_Heatmaps_ViolinPlots/11a_byGene_heatmaps)	| Heatmaps (filtered by chemosensor expression) 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| Antenna_SuppData1_Heatmaps.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_SuppData1_Heatmaps_ViolinPlots/11b_byGene_violinplot) 	 	| Violin Plots (for each chemosensor) 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| 
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_SuppData2_FeaturePlots)	 	| Feature Plots (for each chemosensor) 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| Antenna_SuppData2_FeaturePlots.R 	|
| Github 	| Antenna_receptorList.rds (at least 1 UMI) 	| SeuratObject2_Antenna_mergedBatches_Neurons.rds	| Antenna_SuppData1_Heatmaps.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna/Antenna_SuppData3_PlotsOnAmbiguousNeuralClusters)	 	| Further analysis on cluster removal 	| SeuratObject1_Antenna_mergedBatches_AllCells.rds	| Antenna_Pipeline.R 	|




| TISSUE 	| LOCATION 	| ITEM DESCRIPTION 	| DATA FILE USED (all on [Zenodo](https://zenodo.org/record/5818543)) 	| SCRIPT FILE (all on [Zenodo](https://zenodo.org/record/5818543))	|
|---	|---	|---	|---	|---	|
| Antenna 	| Figure 5B 	| All Cells Marker Heatmap 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure 5C 	| Neurons tSNE 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure 5D 	| Dot plot (manually edited, for full version see Figure S10A) 	| N/A 	| N/A 	|
| Antenna 	| Figure 5E 	| Chord Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure 5F-I 	| Simplified Coexpression Heatmap (Ir64a, Or4, Or82, I41k) 	| SeuratObject_AntennaNeurons.RData 	| SimplifiedCoexpressionHeatmaps.R 	|
| Antenna 	| Figure S9B 	| Doublet Removal Scatter Plot 	| antenna_2.DoubletFinder.08192021.rds antenna_2.DoubletFinder.10062021.rds 	| antenna_2.doubletfinder.R 	|
| Antenna 	| Figure S9C-E 	| Genes, UMIs, Mitochondrial Genes Violin Plots 	| antenna_3.merge.filter.rds 	| antenna_3.integrate.R 	|
| Antenna 	| Figure S9F 	| Batch Correction Feature Plot 	| antenna_3.merge.filter.rds antenna_4.seurat.combined.rds 	| antenna_3.integrate.R antenna_4.anchor.R 	|
| Antenna 	| Figure S9G 	| All Cells tSNE 	| antenna_4.seurat.combined.rds 	| antenna_4.anchor.R 	|
| Antenna 	| Figure S9H 	| Cell Marker Feature Plot 	| SeuratObject_AntennaNeurons.RData 	| FeaturePlots.R 	|
| Antenna 	| Figure S9I 	| Neural Marker Dot Plot 	| SeuratObject_AntennaNeurons.RData 	| FeaturePlots.R 	|
| Antenna 	| Figure S10A 	| Chemosensor Dot Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10B 	| Dot Plot Criteria Scatter Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10C 	| Cluster Coexpression Heatmap 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| Figure S10D 	| Gene Expression Violin Plot (Orco, Ir25a) 	| SeuratObject_AntennaNeurons.RData 	| ViolinPlots.R 	|
| Antenna 	| Figure S10E 	| Cluster Gene Expression Scatter Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10F 	| Cell Gene Expression Scatter Plot 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S10G 	| Gene Expression Venn Diagram 	| antenna_6.neuralClusters.Rdata 	| antenna_7.expressedReceptors.R 	|
| Antenna 	| Figure S11A 	| Gene Expression Violin Plot (Ir41a, Or82, Or4, Ir64a, Gr77, Ir8a, Or47,   Ir41k) 	| SeuratObject_AntennaNeurons.RData 	| ViolinPlots.R 	|
| Antenna 	| Figure S11B 	| Coexpression Heatmap Ir41k 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| Figure S11C-F 	| Simplified Coexpression Heatmap (Ir41a, Ir8a, Gr77, Or47) 	| SeuratObject_AntennaNeurons.RData 	| SimplifiedCoexpressionHeatmaps.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/PseudobulkTables) 	| PseudoBulk Tables 	| SeuratObject_AntennaNeurons.RData 	| PseudoBulkTables.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna_FeaturePlots) 	| Gene Expression Feature Plots 	| SeuratObject_AntennaNeurons.RData 	| FeaturePlots.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna_ViolinPlots) 	| Gene Expression Violin Plots 	| SeuratObject_AntennaNeurons.RData 	| ViolinPlots.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antennae_BroadCoexpressionHeatmaps_byOriginalClusters) 	| Broad Coexpression Heatmaps by Neuron Cluster 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/Antenna_BroadCoexpressionHeatmaps_byChemoreceptor) 	| Broad Coexpression Heatmaps by Chemoreceptor 	| SeuratObject_AntennaNeurons.RData 	| CoexpressionHeatmaps.R 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Raw Expression Matrix (Batch 1    & 2) 	| antenna1_raw_feature_bc_matrix.h5 antenna2_raw_feature_bc_matrix.h5 	| N/A 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Filtered Expression Matrix (Batch 1    & 2) 	| antenna1_filtered_feature_bc_matrix.h5 antenna2_filtered_feature_bc_matrix.h5 	| N/A 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Seurat Object: All Antenna Cells (Batch 1 & 2 combined) 	| SeuratObject_AntennaAllCells.rds 	| N/A 	|
| Antenna 	| [Zenodo](https://zenodo.org/record/5818543) 	| Seurat Object: All Antenna Neurons (Batch 1 & 2 combined) 	| SeuratObject_AntennaNeurons.rds 	| N/A 	|


## *Maxillary Palp: List of Data Items*
|  LOCATION 	| ITEM DESCRIPTION 	| DATA FILE (all on [Zenodo](https://doi.org/10.5281/zenodo.5818951)) 	| SCRIPT FILE (all on [Zenodo](https://doi.org/10.5281/zenodo.5818951))	|
|---	|---	|---	|---	|
| Figure 6B 	| Cell Marker Heatmap  	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Figure 6C 	| tSNE 	| maxpalp_3.neuron.Rdata 	| maxpalp_3.check_receptor.R 	|
| Figure 6D 	| Neuron Gene Heatmap 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Figure 6E 	| Chord Plot 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Figure 6F-M 	| Feature Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Figure S7A 	| Orco Feature Plot (DecontX illustration) 	| maxpalp_1.Decontx.08192021.Rdata 	| maxpalp_1.decontX.R 	|
| Figure S7B-D 	| Genes, UMIs, Mitochondrial Genes Violin Plots 	| SeuratObject_MaxPalpAllCells.RData 	| ViolinPlots.R 	|
| Figure S7E 	| Cell Marker Feature Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Figure S7F 	| Dot Plot Criteria Scatter Plot 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Figure S7G 	| Neural Marker Dot Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Figure S7H-J 	| Genes, UMIs, Mitochondrial Genes Violin Plots (Neuron Clusters) 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| Figure S7K 	| Gene Expression Feature Plot 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| Figure S7L 	| Gene Expression Violin Plot by Cluster (Gr1, Gr2, Gr3, Ir25a, Orco, Or8, Or49, Ir76b) 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| Figure S14M-N 	| Cell Gene Expression Scatter Plot 	| maxpalp_4.neuralClusters.RData 	| maxpalp_5.expressedReceptors.R 	|
| Figure S7O 	| Neuron Cluster Marker Dot Plot 	| SeuratObject_MaxPalpNeurons.RData 	| FeaturePlots.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/PseudobulkTables) 	| PseudoBulk Tables 	| SeuratObject_MaxPalpAllCells.RData 	| PseudoBulkTables.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/MaxPalp_FeaturePlots) 	| Gene Expression Feature Plots 	| SeuratObject_MaxPalpAllCells.RData 	| FeaturePlots.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/MaxPalp_ViolinPlots) 	| Gene Expression Violin Plots 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| [Github](https://github.com/VosshallLab/Younger_Herre_Vosshall2020/tree/main/snRNAseq_SupplementaryData/MaxPalp_ViolinPlots_byCluster) 	| Gene Expression Violin Plots by Cluster 	| SeuratObject_MaxPalpNeurons.RData 	| ViolinPlots.R 	|
| [Zenodo](https://zenodo.org/record/581895) 	| Raw Expression Matrix 	| maxpalp_raw_feature_bc_matrix.h5 	| N/A 	|
| [Zenodo](https://zenodo.org/record/5818952) 	| Filtered Expression Matrix 	| maxpalp_filtered_feature_bc_matrix.h5 	| N/A 	|
| [Zenodo](https://zenodo.org/record/5818952) 	| Seurat   Object: All Maxillary Palp Cells 	| SeuratObject_MaxPalpAllCells.RData 	| N/A 	|
| [Zenodo](https://zenodo.org/record/5818952) 	| Seurat   Object: All Maxillary Palp Neurons 	| SeuratObject_MaxPalpNeurons.RData 	| N/A 	|
