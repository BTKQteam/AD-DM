library(affy)
library(limma)
library(sva)
library(GEOquery)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(dplyr)
library(ggrepel)
library(IOBR)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(dplyr)
library(tidyr)
# load series and platform data from GEO

gset <- getGEO("GSE95849", GSEMatrix =TRUE, getGPL=T)
if (length(gset) > 1) idx <- grep("GPL22448", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

##clear data matrix
matrix<-as.data.frame(exprs(gset))

#create matrix with gene Gene_symbol and expression value
features_new <- fData(gset)
features_new <- select(features_new,Gene_symbol)
matrix <- cbind(features_new,matrix)


#remove duplicate in columns "gene_Gene_symbol" and NA value
matrix<-matrix[matrix$Gene_symbol!="",]
matrix<-matrix[matrix$Gene_symbol!="previous version conserved probe",]

matrix$Gene_symbol <- sub("\\|.*", "", matrix$Gene_symbol)

matrix<-matrix[!duplicated(matrix$Gene_symbol),]
#set "GENE_Gene_symbol" as index
rownames(matrix) <- matrix$Gene_symbol
matrix[, "Gene_symbol"] <- NULL
# log2 transformation
ex <- matrix
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
matrix <- log2(ex) }
groups <- make.names(c("DM","control"))
# box-and-whisker plot
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE95849 (DM vs control)")
boxplot(mat_DM_final, boxwex=0.7, notch=T, main=title, outline=F, las=2, col=DM$color)
legend(, groups, fill=palette(), bty="n")
##Extract sample Info
sampleInfo <- pData(gset)
sampleInfo<-select(sampleInfo,geo_accession,source_name_ch1)
##rename sampleInfo
sampleInfo<-sampleInfo%>%dplyr::rename(ID=geo_accession,phenotype=source_name_ch1)
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("diabetes mellitus","DM",phenotype,fixed = T))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("healthy","control",phenotype,fixed = T))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub('Diabetic peripheral neuropathy','DM_neuro',phenotype,fixed = T))


#subset 
DM <-sampleInfo[(sampleInfo$phenotype=='DM'|sampleInfo$phenotype=='control'),]
DM$color<-ifelse(DM$phenotype=='DM','#619CFF','#F8766D')
#create matrix of that pair
mat_DM <- matrix[,DM$ID]


#BEFORE EXCLUDE OUTLIER, RUN DEG FIRST:

# Convert the phenotype data to a design matrix
design <- model.matrix(~0+ phenotype,data=DM)
# Create a linear model and fit it to the gene expression data
fit <- lmFit(mat_DM, design)
# Create contrasts for the groups of interest
contrast.matrix <- makeContrasts(DM_vs_control = phenotypeDM - phenotypecontrol, levels=colnames(coef(fit)))

fit2 <- contrasts.fit(fit, contrast.matrix)

# Compute moderated t-statistics for each contrast
fit2<- eBayes(fit2,0.01)
full_results <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
logFC <- full_results$logFC
P.Value <- full_results$P.Value
ggplot(data = data.frame(full_results), aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(P.Value >= 0.05, "non-significant",
                                ifelse(logFC < -0.5 & P.Value < 0.05, "significant-control",
                                       ifelse(logFC > 0.5 & P.Value < 0.05, "significant-alzheimer", "grey50")))),
             size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("non-significant" = "grey80",
                                "significant-control" = "darkseagreen",
                                "significant-alzheimer" = "lightcoral")) +guides(color = FALSE) +
  labs(x = "Log2 Fold Change (DM vs Control)", y = "-log10(P-Value)",
       title = "Volcano Plot for DM vs Control (GSE95849)") +theme_bw()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+ theme(panel.border = element_blank())

full_results<-full_results[full_results$P.Value<0.05 & abs(full_results$logFC)>0.5,]
full_results$group<-ifelse(full_results$logFC>0.5,"up","down")
full_results$Gene<-rownames(full_results)
## PCA analysis
#MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(mat_DM))

## Join the PCs to the sample information
cbind(DM, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=phenotype,label=paste("", ID))) + geom_point(size=0.1) +
  geom_text_repel(aes(label = ID), size = 2.5, nudge_x = 1, nudge_y = 1)+
theme_bw()+stat_ellipse()+geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+ theme(panel.border = element_blank())
write.csv(full_results,file = '95849-DM.csv')

outliers<- 'GSM2527027'
#

#subset the matrix excluding outliers

mat_DM_final<-mat_DM[,!colnames(mat_DM)%in%outliers]
DM<-DM[!rownames(DM)%in%outliers,]
# Convert the phenotype data to a design matrix
design2 <- model.matrix(~0+ phenotype,data=DM)
# Create a linear model and fit it to the gene expression data
fit_new <- lmFit(mat_DM_final, design2)

# Create contrasts for the groups of interest
contrast.matrix2 <- makeContrasts(DM_vs_control = phenotypeDM - phenotypecontrol, levels=colnames(coef(fit_new)))

fit2_new <- contrasts.fit(fit_new, contrast.matrix2)
# Compute moderated t-statistics for each contrast
fit2_new <- eBayes(fit2_new,0.01)
full_results_new <- topTable(fit2_new, adjust="fdr", sort.by="B", number=Inf)
full_results_new$Gene<-rownames(full_results_new)

logFC <- full_results_new$logFC
P.Value <- full_results_new$P.Value
ggplot(data = data.frame(fit2_new), aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(P.Value >= 0.05, "non-significant",
                                ifelse(logFC < -0.5 & P.Value < 0.05, "significant-control",
                                       ifelse(logFC > 0.5 & P.Value < 0.05, "significant-alzheimer", "grey50")))),
             size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("non-significant" = "grey80",
                                "significant-control" = "darkseagreen",
                                "significant-alzheimer" = "lightcoral")) +guides(color = FALSE) +
  labs(x = "Log2 Fold Change (DM vs Control)", y = "-log10(P-Value)",
       title = "Volcano Plot for DM vs Control (GSE95849)") +theme_bw()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+ theme(panel.border = element_blank())



full_results_new<-full_results_new[full_results_new$P.Value<0.05 & abs(full_results_new$logFC)>0.5,]
full_results_new$group<-ifelse(full_results_new$logFC>0.5,"up","down")
full_results_new$Gene<-rownames(full_results_new)
write.csv(full_results_new,file = 'gse95849_outlier_full_pval.csv')
## PCA analysis
#MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca2 <- prcomp(t(mat_DM_final))

## Join the PCs to the sample information
cbind(DM, pca2$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=phenotype,label=paste("", ID))) + geom_point(size=0.1) +
  geom_text_repel(aes(label = ID), size = 2.5, nudge_x = 1, nudge_y = 1)+
  theme_bw()+stat_ellipse()+geom_hline(yintercept = 0, linetype = "dashed",lwd=0.3)+
  geom_vline(xintercept = 0, linetype = "dashed",lwd=0.3)+ theme(panel.border = element_blank())+
  scale_color_manual(values = c("#3498DB", "#F39C12"),labels = c("control", "DM"))
#write.csv(full_results_new,file = '95849-DM_outlier.csv')


topN <- 40
annot<-data.frame(group=DM$phenotype,row.names = rownames(DM))

##
gene_names<-full_results_new[order(full_results_new$logFC,decreasing = T),]
up_gene<-gene_names$Gene[1:20]
down_gene<-gene_names$Gene[(nrow(gene_names)-19):nrow(gene_names)]
gene_list<-c(up_gene,down_gene)
gene_matrix <- mat_DM_final[gene_list,]
pheatmap(gene_matrix,
         labels_row = rownames(gene_matrix),
         scale="row",annotation_col = annot,
         show_colnames = FALSE)
####WGCNA analysis
library('WGCNA')
library('gridExtra')
library("CorLevelPlot")

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers, and transpose the matrix of genes as columns
power <- c(c(1:20), seq(from = 22, to = 30, by = 2))
W_matrix <-t(mat_DM_final)
# Call the network topology analysis function
sft <- pickSoftThreshold(W_matrix,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +geom_vline(xintercept = 11, color = 'red')+
  theme_classic()


print(a1+a2)
# convert matrix to numeric
W_matrix[] <- sapply(W_matrix, as.numeric)

soft_power <- 7
temp_cor <- cor
cor <- WGCNA::cor
# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(W_matrix), method = "average");
# plot sample tree

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(W_matrix,
                          maxBlockSize = 25000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# Correlate traits --------------------------------------------------------
# create traits file - binarize categorical variables
traits <- DM %>% 
  mutate(DM = ifelse(grepl('DM', phenotype), 1, 0)) %>% select(3) #(8) means the column "AD" is the 8th columns of 'hippocampus'

# binarize categorical variables


control <-DM %>% 
  mutate(control = ifelse(grepl('control', phenotype), 1, 0)) %>% select(3)
traits <- cbind(traits, control)


# Define numbers of genes and samples
nSamples <- nrow(W_matrix)
nGenes <- ncol(W_matrix)  


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

names(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[43:44],
             y = names(heatmap.data)[1:42],
             col = c("blue1", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)
turquoise<-module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()
write.csv(turquoise,file = 'turquoise-DM-95849.csv',row.names = F)



#####################################
#### IMMUNE INFILTRATION ############
#####################################

signature_score_calculation_methods
names(signature_tme)[1:20]
names(signature_metabolism)[1:20]

mat_DM_final[1:5, 1:10]



#Analyze all collected signature scores (integrating three methods: PCA, ssGSEA and z-score).

sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = mat_DM_final,
                             signature       = signature_collection,
                             method          = "pca",
                             mini_gene_count = 2)
##The signature gene sets derived from GO, KEGG, HALLMARK and REACTOME datasets.
sig_hallmark<-calculate_sig_score(pdata           = NULL,
                                  eset            = mat_DM_final,
                                  signature       = hallmark,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)

sig_kegg<-calculate_sig_score(pdata           = NULL,
                              eset            = mat_DM_final,
                              signature       = kegg,
                              method          = "ssgsea",
                              mini_gene_count = 2)
sig_go<-calculate_sig_score(pdata           = NULL,
                            eset            = mat_DM_final,
                            signature       = go_bp,
                            method          = "ssgsea",
                            mini_gene_count = 2)
###METHOD 1: CIBERSORT

cibersort<-deconvo_tme(eset = matrix, method = "cibersort_abs", arrays = FALSE, perm = 200 )
res<-cell_bar_plot(input = cibersort[1:18,], title = "CIBERSORT Cell Fraction")


##METHOD 2: ESTIMATE
estimate<-deconvo_tme(eset = mat_DM_final, method = "estimate")
res_es<-cell_bar_plot(input = estimate[1:11,], title = "ESTIMATE Cell Fraction")

##METHOD 3: MCP
mcp<-deconvo_tme(eset = matrix, method = "mcpcounter")

##METHOD 4: xCELL
xcell<-deconvo_tme(eset = matrix, method = "xcell",arrays = FALSE)
##METHOD 5: EPIC
epic<-deconvo_tme(eset = matrix, method = "epic", arrays = FALSE)
####5.2 Phenotype Module
DM$ID<-NULL
DM <- tibble::rownames_to_column(DM, var = "ID")

#unique(sampleInfo$group)

#match<-sam
##id1 and id2 must match!
res<-iobr_cor_plot(pdata_group           = km_data_all,
                   id1                   = "ID",
                   feature_data          = xcell,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "risk",
                   is_target_continuous  = FALSE,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   signature_group       = (sig_group)[23],
                   ProjectID             = "DM-xcell_PIP4K2A",
                   
                   palette_box           = "aaas",
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 4,
                   feature_limit         = 26,
                   character_limit       = 30,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE)
?iobr_cor_plot





#############################
#############################
####### METABOLISM ##########
#############################
km_data_all<-rownames_to_column(km_data_all,var='ID')


sig_meta<-calculate_sig_score(pdata           = NULL,
                              eset            = matrix,
                              signature       = signature_metabolism,
                              method          = "ssgsea",
                              mini_gene_count = 2)

res_meta<-iobr_cor_plot(pdata_group           = km_data_all,
                   id1                   = "ID",
                   feature_data          = sig_meta,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "risk",
                   is_target_continuous  = F,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   signature_group       = (sig_group)[c(28,30,31)],
                   ProjectID             = "DM-metabolism",
                   
                   palette_box           = "aaas",
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 2,
                   feature_limit         = 200,
                   character_limit       = 50,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE)
write.csv(res_meta,file = 'DM_meta_list_PIP4K2A.csv')
write.csv(sig_meta,file = 'DM_meta_score.csv')
###test
pdata_group<-rownames_to_column(as.data.frame(t(mat_DM_final)),var = "ID")
pdata_group<-as.data.frame(pdata_group[,c("ID","LHX2")])
head(pdata_group)

res<-iobr_cor_plot(pdata_group           = pdata_group,
                   id1                   = "ID",
                   feature_data          = sig_meta,
                   id2                   = "ID",
                   target                = "LHX2",
                   group                 = "group3",
                   is_target_continuous  = T,
                   padj_cutoff           = 1,
                   index                 = 3,
                   category              = "signature",
                   signature_group       = (sig_group)[c(28,30,31)],
                   ProjectID             = "DM-metabolism_test",
                   
                   palette_box           = "jco",
                   palette_corplot       = "normal",
                   palette_heatmap       = 2,
                   feature_limit         = 30,
                   character_limit       = 30,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE)


######lay toan bo samples:
selected_genes<-c('PIP4K2A')
mat_all_transpose<-t(mat_DM_final)
mat_all<-as.data.frame(mat_all_transpose[,selected_genes])
km_data_all<-data.frame(mat_all)
colnames(km_data_all)<-'PIP4K2A'
median_risk_all<-median(km_data_all$PIP4K2A)
km_data_all$risk <- ifelse(km_data_all$PIP4K2A <= median_risk_all, "low", "high")
km_data_all$group<-DM$phenotype
go_subset<-select(sig_go,GOBP_REGULATION_OF_INSULIN_RECEPTOR_SIGNALING_PATHWAY,
                  GOBP_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,
                  GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS,
                  GOBP_NEGATIVE_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY,
                  GOBP_NEGATIVE_REGULATION_OF_LIPID_KINASE_ACTIVITY,
                  GOBP_PHOSPHATIDYLINOSITOL_PHOSPHATE_BIOSYNTHETIC_PROCESS,
                  GOBP_VESICLE_MEDIATED_TRANSPORT,
                  GOBP_MACROAUTOPHAGY,
                  GOBP_CYTOSOLIC_TRANSPORT,
                  GOBP_REGULATION_OF_VESICLE_FUSION,
                  GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGOSOME_ASSEMBLY,
                  GOBP_POSITIVE_REGULATION_OF_VACUOLE_ORGANIZATION,
                  GOBP_POSITIVE_REGULATION_OF_ORGANELLE_ASSEMBLY,
                  GOBP_REGULATION_OF_AUTOPHAGY,
                  GOBP_POSITIVE_REGULATION_OF_MACROAUTOPHAGY,
                  GOBP_AUTOPHAGOSOME_LYSOSOME_FUSION)
km_data_all<-cbind(km_data_all,go_subset)

### boxplot
my_comparisons<- c("control", "DM")
km_data_all$group <- factor(km_data_all$group, levels = c("control", "DM"))

ggboxplot(km_data_all, x = "group", y = "GOBP_NEGATIVE_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
          color = "black", fill = 'group', palette = c("#00AFBB",'#E7B800', "#FC4E07"),
          add = "jitter", shape = "group",add.params = list(fill = "white"))+ stat_compare_means(method = 't.test')

