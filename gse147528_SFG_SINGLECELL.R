library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SingleR)
library(celldex)
library(scRNAseq)
library(ggpubr)
library(hdf5r)
library(SCP)
#get dir
h5_files <- list.files(pattern = "*.h5")
h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, CreateSeuratObject)
h5_seurat  <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], 
                    add.cell.ids = c("SFG2", "SFG1", "SFG3","SFG4","SFG6","SFG7","SFG5","SFG9","SFG8","SFG10"), project = "project")
view(h5_seurat@meta.data)
# create a sample column
h5_seurat$sample <- rownames(h5_seurat@meta.data)
# split sample column
h5_seurat@meta.data <- separate(h5_seurat@meta.data, col = 'sample', into = c('Patient', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
h5_seurat$mitoPercent <- PercentageFeatureSet(h5_seurat, pattern='^MT-')

saveRDS(h5_seurat,file = 'gse147528__SFG_raw.rds')
# filtering
merged_seurat_filtered <- subset(h5_seurat, subset = nCount_RNA > 500 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 20)

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered,nfeatures=2000)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered,resolution=1)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- RunTSNE(object = merged_seurat_filtered, dims = 1:20)

# Create a new column named 'group' in the metadata
merged_seurat_filtered@meta.data$group <- ifelse(merged_seurat_filtered@meta.data$Patient %in% c('SFG1', 'SFG2', 'SFG3'), 'braak 0',
                                                 ifelse(merged_seurat_filtered@meta.data$Patient %in% c('SFG4', 'SFG5', 'SFG6', 'SFG7'), 'braak 2',
                                                        'braak 6'))
saveRDS(merged_seurat_filtered,file = 'gse147528_sfg_filtered.rds')



#Read RDS object
merged_seurat_filtered<-gse147528_sfg_filtered
# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'group',
              cols = c('#13d18c','#ffcc33','#df6f75'))

DimPlot(merged_seurat_filtered, reduction = 'tsne', group.by = 'seurat_clusters',label = TRUE,
        label.size = 4, label.color = "black")

FeaturePlot(merged_seurat_filtered, reduction = 'umap', features = c('MARCH1','AP2S1','ATOX1','ATP1A1','ATPIF1','ATRX','BAG6','BLOC1S1','BLVRB','CDC40','CGGBP1','CPA3','CRYBA2','DHRS7B','DOCK3','DPCD','EDEM3','GDE1','GHITM','GLS','GPI','H2AFY2','HMGCR','ICA1','IL6ST','KCNE4','KIAA0754','LHX2','MAP7D2','MOCS2','MRPS22','MRPS28','MTSS1','NAPA','NDUFAF3','NOTCH2NL','ORC4','ORC5','PEF1','PRR3','PSMB3','SKI','SLIRP','SNAPC3','SNCA','SNX1','SPHK2','SRGAP1','TRIM58','TSTA3','TUBA4A','UBLCP1','WBSCR17','YWHAZ','ZFP36L1','ZNF275','ZNF721'))
FeaturePlot(merged_seurat_filtered, reduction = 'umap', features = c('ELAVL4'))
FeaturePlot(merged_seurat_filtered, reduction = 'umap', features = c('PIP4K2A','BAZ1A','ELAVL4','FGD4','MAP7D2','NOTCH2NL','SNAP23','ZFP36L1'))

b<-DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'cluster')
b

OPC = c(8)
Astrocyte = c(7,17,19,21,31)
Microglia = c(9)
Oligo=c(0,1)

Endothelial = c(24)

Excit = c(2,5,6,10,11,12,13,15,16,20,22,23,25,26,27,28,29)
Inhibit = c(3,4,14,18,30)


cl.mat <- cbind(c( rep("OPC", length(OPC)), rep("Astrocyte", length(Astrocyte)), rep("Microglia", length(Microglia)),rep("Oligo", length(Oligo)),
                   rep("Endothelial", length(Endothelial)), rep("Excitatory Neuron", length(Excit)), rep("Inhibitory Neuron", length(Inhibit))),
                c(OPC, Astrocyte, Microglia,Oligo, Endothelial, Excit, Inhibit))

cl.vec <- merged_seurat_filtered$RNA_snn_res.1
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
merged_seurat_filtered$cluster <- ct.vec

merged_seurat_filtered$cluster <- factor(merged_seurat_filtered$cluster, levels = c('OPC', 'Astrocyte','Microglia','Oligo', 'Endothelial', 
                                                                                    'Excitatory Neuron', 'Inhibitory Neuron'))
DimPlot(merged_seurat_filtered, reduction = "umap", label = TRUE, pt.size = 0.2,group.by ='cluster',cols='Set2',label.size = 2,label.box = T,repel = T) 
comparisons <- list(c('braak 0','braak 2'),c('braak 2','braak 6'),c('braak 0','braak 6'))
Idents(object = merged_seurat_filtered)<-merged_seurat_filtered$cluster

VlnPlot(merged_seurat_filtered, features = c("ELAVL4"),pt.size = 0,
        idents = c('Inhibitory Neuron'),group.by = 'group', split.by = "group",cols = c('#13d18c','#ffcc33','#df6f75')) +
  stat_compare_means(comparisons = comparisons,label = "p.signif")+ geom_boxplot(width=.1)+
  ylim(0,8)+color_palette('Set1')

######volcano
markers <- FindAllMarkers(merged_seurat_filtered , min.pct = 0.25,test.use = "wilcox")
Oligo.MAker<-markers[markers$cluster=='Oligo',]
# add a column of NAs
Oligo.MAker$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Oligo.MAker$diffexpressed[Oligo.MAker$avg_log2FC > 0.5 & Oligo.MAker$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Oligo.MAker$diffexpressed[Oligo.MAker$avg_log2FC < -0.5 & Oligo.MAker$p_val_adj < 0.05] <- "DOWN"

Oligo.MAker$delabel <- ifelse(Oligo.MAker$gene %in% c('PIP4K2A'), Oligo.MAker$gene, NA)



ggplot(data = data.frame(Oligo.MAker), aes(x = avg_log2FC, y = -log10(p_val_adj),col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + 
  geom_point(size = 0.3, alpha = 0.8) +
  guides(color = guide_legend(override.aes = list(size = 2)))+ theme_bw()+
  scale_color_manual(values = c("darkseagreen", "grey70", "lightcoral"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'group', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +  
  ggtitle('SFG-GSE147528')+
  # Plot title 
  geom_text_repel(max.overlaps = 20,size = 2.5,colour='gray25')+ # To show all labels 
  theme(panel.border = element_blank())


