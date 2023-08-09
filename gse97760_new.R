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

gset <- getGEO("GSE97760", GSEMatrix =TRUE, getGPL=T)
if (length(gset) > 1) idx <- grep("GPL16699", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

##clear data matrix
matrix<-as.data.frame(exprs(gset))

#create matrix with gene GENE_SYMBOL and expression value
features_new <- fData(gset)
features_new <- select(features_new,GENE_SYMBOL)
matrix <- cbind(features_new,matrix)


#remove duplicate in columns "gene_GENE_SYMBOL" and NA value
matrix<-matrix[matrix$GENE_SYMBOL!="",]
matrix$GENE_SYMBOL <- gsub(" ///.*", "", matrix$GENE_SYMBOL)

matrix<-matrix[!duplicated(matrix$GENE_SYMBOL),]
#set "GENE_GENE_SYMBOL" as index
rownames(matrix) <- matrix$GENE_SYMBOL
matrix[, "GENE_SYMBOL"] <- NULL


# log2 transformation
ex <- matrix
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
matrix <- log2(ex) }


##Extract sample Info
sampleInfo <- pData(gset)

##Extract sample Info
sampleInfo <- pData(gset)
sampleInfo<-select(sampleInfo,geo_accession,characteristics_ch1.2)
##rename sampleInfo
sampleInfo<-sampleInfo%>%dplyr::rename(ID=geo_accession,phenotype=characteristics_ch1.2)
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("disease: healthy","control",phenotype,fixed = T))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("disease: advanced Alzheimer's disease","AD",phenotype,fixed = T))

#BEFORE EXCLUDE OUTLIER, RUN DEG FIRST:

# Convert the phenotype data to a design matrix
design <- model.matrix(~0+ phenotype,data=sampleInfo)
# Create a linear model and fit it to the gene expression data
fit <- lmFit(matrix, design)
# Create contrasts for the groups of interest
contrast.matrix <- makeContrasts(AD_vs_control = phenotypeAD - phenotypecontrol, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)

# Compute moderated t-statistics for each contrast
fit2<- eBayes(fit2,0.01)
full_results <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
full_results$Gene<-rownames(full_results)

logFC <- full_results$logFC
Pval <- full_results$P.Value
ggplot(data = data.frame(full_results), aes(x = logFC, y = -log10(Pval))) +
  geom_point(aes(color = ifelse(Pval >= 0.05, "non-significant",
                                ifelse(logFC < -0.5 & Pval < 0.05, "significant-control",
                                       ifelse(logFC > 0.5 & Pval < 0.05, "significant-alzheimer", "grey50")))),
             size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("non-significant" = "grey80",
                                "significant-control" = "darkseagreen",
                                "significant-alzheimer" = "lightcoral")) +guides(color = FALSE) +
  labs(x = "Log2 Fold Change (AD vs Control)", y = "-log10(P-Value)",
       title = "Volcano Plot for AD vs Control (GSE97760)") +theme_bw()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+ theme(panel.border = element_blank())

full_results<-full_results[full_results$P.Value<0.05 & abs(full_results$logFC)>0.5,]
full_results$group<-ifelse(full_results$logFC>0.5,"up","down")


topN <- 40
##
ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Gene)
## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- as.matrix(matrix[ids_of_interest,])
annot<-data.frame(group=sampleInfo$phenotype,row.names = rownames(sampleInfo))
pheatmap(gene_matrix,
         labels_row = rownames(gene_matrix),
         scale="row",annotation_col = annot,
         show_colnames = FALSE)


## PCA analysis
#MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(matrix))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=phenotype,label=paste("", ID))) + geom_point(size=0.1) +
  geom_text_repel(aes(label = ID), size = 2.5, nudge_x = 1, nudge_y = 1)+
  theme_bw()+stat_ellipse()+geom_hline(yintercept = 0, linetype = "dashed",lwd=0.3)+
  geom_vline(xintercept = 0, linetype = "dashed",lwd=0.3)+ theme(panel.border = element_blank())+
  scale_color_manual(values = c("#F39C12","#3498DB"))
#write.csv(full_results,file = '97760_Pval.csv')



#####################################
#### IMMUNE INFILTRATION ############
#####################################

signature_score_calculation_methods
names(signature_tme)[1:20]
names(signature_metabolism)[1:20]

matrix[1:5, 1:10]



#Analyze all collected signature scores (integrating three methods: PCA, ssGSEA and z-score).

sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = matrix,
                             signature       = signature_collection,
                             method          = "pca",
                             mini_gene_count = 2)
##The signature gene sets derived from GO, KEGG, HALLMARK and REACTOME datasets.
sig_hallmark<-calculate_sig_score(pdata           = NULL,
                                  eset            = matrix,
                                  signature       = hallmark,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)

sig_meta<-calculate_sig_score(pdata           = NULL,
                              eset            = matrix,
                              signature       = signature_metabolism,
                              method          = "zscore",
                              mini_gene_count = 2)
sig_kegg<-calculate_sig_score(pdata           = NULL,
                              eset            = matrix,
                              signature       = kegg,
                              method          = "ssgsea",
                              mini_gene_count = 2)
sig_go<-calculate_sig_score(pdata           = NULL,
                              eset            = matrix,
                              signature       = go_bp,
                              method          = "ssgsea",
                              mini_gene_count = 2)

x###METHOD 1: CIBERSORT

cibersort<-deconvo_tme(eset = matrix, method = "cibersort", arrays = FALSE, perm = 200 )
res<-cell_bar_plot(input = cibersort[1:19,], title = "CIBERSORT Cell Fraction of AD")
?cell_bar_plot

##METHOD 2: ESTIMATE
estimate<-deconvo_tme(eset = matrix, method = "estimate")
res_es<-cell_bar_plot(input = estimate[1:19,], title = "ESTIMATE Cell Fraction")
##METHOD 3: MCP
mcp<-deconvo_tme(eset = matrix, method = "mcpcounter")

##METHOD 4: xCELL
xcell<-deconvo_tme(eset = matrix, method = "xcell",arrays = FALSE)
##METHOD 5: EPIC
epic<-deconvo_tme(eset = matrix, method = "epic", arrays = FALSE)
km_data_all<-rownames_to_column(km_data_all,var='ID')
####5.2 Phenotype Module
##id1 and id2 must match!
res_meta<-iobr_cor_plot(pdata_group     = km_data_all,
                   id1                   = "ID",
                   feature_data          = sig_meta,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "risk",
                   is_target_continuous  = FALSE,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   signature_group       = (sig_group)[c(28,30,31)],
                   ProjectID             = "AD-PIP4K2A-meta",
                   
                   palette_box           = "jco",
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 2,
                   feature_limit         = 200,
                   character_limit       = 50,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE)
write.csv(res_meta,file = '97760_meta_list_PIP4K2A.csv')
write.csv(sig_go,file = '97760_go_score_PIP4K2A.csv')
##### AUC ###########
#####################
library(pROC)
risk_list<-read.csv('intersection_AD_DM_new.csv',header = F)
risk_list<-unlist(risk_list)
mat_AUC<-t(matrix[risk_list,])
mat_AUC<-as.data.frame(mat_AUC)
mat_AUC <- mat_AUC[, !colSums(is.na(mat_AUC)) > 0]

label<-as.factor(ifelse(sampleInfo$phenotype=='AD',1,0))


auc_scores <- data.frame(Gene = colnames(mat_AUC), AUC = NA)  # Create an empty dataframe
# Create an empty plot
#gene_to_plot<-c(1,5,16,19,20,35,49,73,77,84,94,97) #('STOML1','ASB16','HDAC7','SCAMP2','OS9','ZFP36L1','F5','PARP9','FAM20B','BMP2K','MLEC','BACE2)
# Create an empty plot
plot(NULL, xlim = c(1, 0), ylim = c(0, 1), xlab = "Specificity", ylab = "Sensitivity", main = "BACE2")

for (i in 1:ncol(mat_AUC)) {
  gene <- mat_AUC[, i]
  roc_obj <- roc(label, gene)
  auc_scores$AUC[i] <- auc(roc_obj)
  # Calculate confidence interval
  ci <- ci.auc(roc_obj, conf.level = 0.95)
  auc_scores$CI[i] <- paste(ci[1], ci[2], sep = " , ")
}
write.csv(auc_scores,file = '97760_auc_new_july.csv')
##plot specific genes
gene <- mat_AUC[, 'PIP4K2A']
roc_obj <- roc(label, gene)
plot.roc(roc_obj, add = TRUE,print.auc=T,plot=T,smooth=T,col='red4')
lines(c(0, 1), c(1, 0), col = "gray", lty = "dashed",add=T)


mat_AD<-cbind(mat_AUC,label)
AD<-sampleInfo[sampleInfo$phenotype=='AD',]
mat_AD_meta<-matrix[,rownames(AD)]
#######################################################
############  construct the DMAD signature ############
#######################################################

AD<-sampleInfo[sampleInfo$phenotype=='AD',]
mat_AD<-matrix[,rownames(AD)]
mat_AD_transpose<-t(mat_AD)
selected_genes<-c('PIP4K2A')

mat_AD<-as.data.frame(mat_AD_transpose[,selected_genes])
km_data<-data.frame(mat_AD)
colnames(km_data)<-'PIP4K2A'


median_risk<-median(km_data$PIP4K2A)


km_data$risk <- ifelse(km_data$PIP4K2A <= median_risk, "low", "high")

######lay toan bo samples:
selected_genes<-c('PIP4K2A')
mat_all_transpose<-t(matrix)
mat_all<-as.data.frame(mat_all_transpose[,selected_genes])
km_data_all<-data.frame(mat_all)
colnames(km_data_all)<-'PIP4K2A'
median_risk_all<-median(km_data_all$PIP4K2A)
km_data_all$risk <- ifelse(km_data_all$PIP4K2A <= median_risk_all, "low", "high")
km_data_all$group<-sampleInfo$phenotype
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
my_comparisons<- c("control", "AD")
km_data_all$group <- factor(km_data_all$group, levels = c("control", "AD"))

ggboxplot(km_data_all, x = "group", y = "GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS",
          color = "black", fill = 'group', palette = c("#00AFBB", "#FC4E07","#E7B800"),
          add = "jitter", shape = "group",add.params = list(fill = "white"))+ stat_compare_means(method = 't.test')

