library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(SingleR)
library(celldex)
library(scRNAseq)
library(ggpubr)
library(hdf5r)



# Initialize the first Seurat object with the first sample
first_sample_folder <- list.dirs(recursive = FALSE)[1]
count_matrix_file <- list.files(first_sample_folder, pattern = "\\.txt$", full.names = TRUE)
count_matrix <- read.table(count_matrix_file, header = TRUE, row.names = 1)
seurat_obj <- CreateSeuratObject(counts = count_matrix)
cell_id <- tools::file_path_sans_ext(basename(count_matrix_file))
# Rename cells with the file name
seurat_obj <- RenameCells(object = seurat_obj, add.cell.id = cell_id)
# Loop through each sample folder
sample_folders <- list.dirs(recursive = FALSE)[-1]

for (sample_folder in sample_folders) {
  # Get the file paths for the count matrix and barcode files
  count_matrix_file <- list.files(sample_folder, pattern = "\\.txt$", full.names = TRUE)
  # Read the count matrix and barcode files
  count_matrix <- read.table(count_matrix_file, header = TRUE, row.names = 1)
  # Create a new Seurat object for the current sample
  sample_seurat <- CreateSeuratObject(counts = count_matrix)
  cell_id <- tools::file_path_sans_ext(basename(count_matrix_file))
  # Rename cells with the file name
  sample_seurat <- RenameCells(object = sample_seurat, add.cell.id = cell_id)
  # Integrate the current sample into the main Seurat object
  seurat_obj <- merge(seurat_obj, y = sample_seurat)
}


seurat_obj$sample <- rownames(seurat_obj@meta.data)
# split sample column
seurat_obj@meta.data <- separate(seurat_obj@meta.data, col = 'sample', into = c('GSE', 'Sample','library','barcode'), 
                                sep = '_')

# calculate mitochondrial percentage
seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern='^MT-')

saveRDS(seurat_obj,file = 'gse198323_raw.rds')
# filtering
merged_seurat_filtered <- subset(seurat_obj, subset = nCount_RNA > 500 &
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
merged_seurat_filtered@meta.data$group <- ifelse(merged_seurat_filtered@meta.data$Sample %in% c('Sample23', 'Sample25', 'Sample27','Sample32','Sample33','Sample35','Sample36','Sample37'), 
                                                 'AD','control')
saveRDS(merged_seurat_filtered,file = 'gse198323_hippocampus_filtered.rds')


#Read RDS object
merged_seurat_filtered<-gse198323_hippocampus_filtered

# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Sample')
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'group',
              cols = c('#f6c943','#9c7bd6'))

DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE,
        label.size = 4, label.color = "black")

FeaturePlot(merged_seurat_filtered, reduction = 'umap', features = c('PIP4K2A','BAZ1A','ELAVL4','FGD4','MAP7D2','NOTCH2NL','SNAP23','ZFP36L1'))

b<-DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'cluster')
b

Choroid=c(24)
OPC = c(8,37)
Astrocyte = c(6,13,28,31,34,35,36)
Microglia = c(18)
Oligo=c(4,25,29)
Ependymal=c(23)
Endothelial = c(26)

Glut.N = c(0,1,2,3,5,7,9,11,12,14,17,20,21,22,30,32,33)
GABA.N = c(10,15,16,19,27)

cl.mat <- cbind(c( rep("OPC", length(OPC)), rep("Astrocyte", length(Astrocyte)), rep("Microglia", length(Microglia)),rep("Oligo", length(Oligo)),
                   rep("Endothelial", length(Endothelial)), rep("Glut.N", length(Glut.N)), rep("GABA.N", length(GABA.N)),
                rep("Ependymal", length(Ependymal)), rep("Choroid", length(Choroid))),
                c(OPC, Astrocyte, Microglia,Oligo, Endothelial, Glut.N, GABA.N,Ependymal,Choroid))

cl.vec <- merged_seurat_filtered$RNA_snn_res.1
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
merged_seurat_filtered$cluster <- ct.vec

merged_seurat_filtered$cluster <- factor(merged_seurat_filtered$cluster, levels = c('OPC', 'Astrocyte','Microglia','Oligo', 'Endothelial', 
                                                                                    'Glut.N', 'GABA.N','Ependymal','Choroid'))
DimPlot(merged_seurat_filtered, reduction = "umap", label = TRUE, pt.size = 0.5,group.by ='cluster',cols='Set2',label.size = 2,label.box = T,repel = T) 
comparisons <- list(c('control','AD'))
merged_seurat_filtered$group <- factor(merged_seurat_filtered$group, levels = c("control", "AD"))


Idents(object = merged_seurat_filtered)<-merged_seurat_filtered$cluster

VlnPlot(merged_seurat_filtered, features = c("ELAVL4"),pt.size = 0,
        idents = c('GABA.N'), group.by = 'group', split.by = "group",cols = c('#f6c943','#9c7bd6')) +
  stat_compare_means(comparisons = comparisons,label = "p.signif")+ geom_boxplot(width=.1)+
  ylim(0,7)



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
  ggtitle('hippocampus-GSE198323')+
  # Plot title 
  geom_text_repel(max.overlaps = 20,size = 2.5,colour='gray25')+ # To show all labels 
  theme(panel.border = element_blank())


