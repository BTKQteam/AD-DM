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

gset <- getGEO("GSE48350", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

##clear data matrix
matrix<-as.data.frame(exprs(gset))

#create matrix with gene Gene.Symbol and expression value
features_new <- fData(gset)
features_new <- select(features_new,Gene.symbol)
matrix <- cbind(features_new,matrix)


#remove duplicate in columns "gene_Gene.Symbol" and NA value
matrix<-matrix[matrix$Gene.symbol!="",]
matrix$Gene.symbol <- sub("///.*", "", matrix$Gene.symbol)

matrix<-matrix[!duplicated(matrix$Gene.symbol),]
#set "GENE_Gene.Symbol" as index
rownames(matrix) <- matrix$Gene.symbol
matrix[, "Gene.symbol"] <- NULL
# log2 transformation
ex <- matrix
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
matrix <- log2(ex) }

##Extract sample Info
sampleInfo <- pData(gset)
sampleInfo<-select(sampleInfo,geo_accession,source_name_ch1)
##rename sampleInfo
sampleInfo<-sampleInfo%>%dplyr::rename(ID=geo_accession,tissue=source_name_ch1)
sampleInfo <- sampleInfo %>%
  mutate(group = case_when(
    grepl("hippocampus_", tissue) ~ "AD_Hippocampus",
    grepl("post-central gyrus_", tissue) ~ "AD_Post_central_gyrus",
    grepl("superior frontal gyrus_", tissue) ~ "AD_superior_frontal_gyrus",
    grepl("entorhinal cortex_", tissue) ~ "AD_Entorhinal_cortex",
    grepl("brain, postcentral gyrus", tissue) ~ "Control_Post_central_gyrus",
    grepl("brain, superior frontal gyrus", tissue) ~ "Control_superior_frontal_gyrus",
    grepl("brain, hippocampus", tissue) ~ "Control_Hippocampus",
    grepl("brain, entorhinal cortex", tissue) ~ "Control_Entorhinal_cortex"
  ))

#subset 
sfg <-sampleInfo[(sampleInfo$group=='AD_superior_frontal_gyrus'|sampleInfo$group=='Control_superior_frontal_gyrus'),]
#create matrix of that pair
mat_sfg <- matrix[,sfg$ID]

#BEFORE EXCLUDE OUTLIER, RUN DEG FIRST:

# Convert the phenotype data to a design matrix
design <- model.matrix(~0+ group,data=sfg)
# Create a linear model and fit it to the gene expression data
fit <- lmFit(mat_sfg, design)
# Create contrasts for the groups of interest
contrast.matrix <- makeContrasts(sfg_vs_control = groupAD_superior_frontal_gyrus - groupControl_superior_frontal_gyrus, levels=colnames(coef(fit)))

fit2 <- contrasts.fit(fit, contrast.matrix)

# Compute moderated t-statistics for each contrast
fit2<- eBayes(fit2,0.01)
full_results <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
full_results$Gene<-rownames(full_results)

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
  labs(x = "Log2 Fold Change (AD-hippo vs Control)", y = "-log10(P-Value)",
       title = "Volcano Plot for AD-sfg vs Control (GSE48350)") +theme_bw()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+ theme(panel.border = element_blank())

full_results<-full_results[full_results$P.Value<0.05 & abs(full_results$logFC)>0.5,]
full_results$group<-ifelse(full_results$logFC>0.5,"up","down")


## PCA analysis
#MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(mat_sfg))

## Join the PCs to the sample information
cbind(sfg, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("", ID))) + geom_point(size=0.1) + geom_text_repel(aes(label = ID), size = 2, nudge_x = 1, nudge_y = 1)+theme_bw()+stat_ellipse()+geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+ theme(panel.border = element_blank())


outliers<- c('GSM300250','GSM300191','GSM300213','GSM300251','GSM300288','GSM300311',
             'GSM318840','GSM300284','GSM300201','GSM300188','GSM300199','GSM300191',
             'GSM300229','GSM300195','GSM300180','GSM300217')
#subset the matrix excluding outliers

mat_sfg_final<-mat_sfg[,!colnames(mat_sfg)%in%outliers]
sfg_new<-sfg[!rownames(sfg)%in%outliers,]
# Convert the phenotype data to a design matrix
design2 <- model.matrix(~0+ group,data=sfg_new)
# Create a linear model and fit it to the gene expression data
fit_new <- lmFit(mat_sfg_final, design2)

# Create contrasts for the groups of interest
contrast.matrix2 <- makeContrasts(sfg_vs_control = groupAD_superior_frontal_gyrus - groupControl_superior_frontal_gyrus, levels=colnames(coef(fit_new)))

fit2_new <- contrasts.fit(fit_new, contrast.matrix2)
# Compute moderated t-statistics for each contrast
fit2_new <- eBayes(fit2_new,0.01)
full_results_new <- topTable(fit2_new, adjust="fdr", sort.by="B", number=Inf)
logFC <- full_results_new$logFC
P.value <- full_results_new$P.Value
ggplot(data = data.frame(fit2_new), aes(x = logFC, y = -log10(P.value))) +
  geom_point(aes(color = ifelse(P.value >= 0.05, "non-significant",
                                ifelse(logFC < -0.5 & P.value < 0.05, "significant-control",
                                       ifelse(logFC > 0.5 & P.value < 0.05, "significant-alzheimer", "grey50")))),
             size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("non-significant" = "grey80",
                                "significant-control" = "darkseagreen",
                                "significant-alzheimer" = "lightcoral")) +guides(color = FALSE) +
  labs(x = "Log2 Fold Change (sfg vs Control)", y = "-log10(P-Value)",
       title = "Volcano Plot for sfg vs Control (GSE15932)") +theme_bw()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed")+ theme(panel.border = element_blank())


full_results_new<-full_results_new[full_results_new$P.Value<0.05 & abs(full_results_new$logFC)>0.5,]
full_results_new$group<-ifelse(full_results_new$logFC>0.5,"up","down")
full_results_new$Gene<-rownames(full_results_new)
write.csv(full_results_new,file = '48350_SGF_pval.csv')

pca <- prcomp(t(mat_sfg_final))

## Join the PCs to the sample information
cbind(sfg_new, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("", ID))) + geom_point(size=0.1) + geom_text_repel(aes(label = ID), size = 2, nudge_x = 1, nudge_y = 1)+theme_bw()+stat_ellipse()+geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+ theme(panel.border = element_blank())

##### AUC ###########
#####################
library(pROC)
risk_list<-read.csv('intersection_AD_DM_new.csv',header = F)
risk_list<-unlist(risk_list)
mat_AUC<-t(mat_sfg_final[risk_list,])

mat_AUC<-as.data.frame(mat_AUC)
mat_AUC <- mat_AUC[, !colSums(is.na(mat_AUC)) > 0]

label<-as.factor(ifelse(sfg_new$group=='AD_superior_frontal_gyrus',1,0))


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
write.csv(auc_scores,file = '48350_auc_new.csv')




######lay toan bo samples:
selected_genes<-c('PIP4K2A')
mat_all_transpose<-t(mat_sfg_final)
mat_all<-as.data.frame(mat_all_transpose[,selected_genes])
km_data_all<-data.frame(mat_all)
colnames(km_data_all)<-'PIP4K2A'
median_risk_all<-median(km_data_all$PIP4K2A)
km_data_all$risk <- ifelse(km_data_all$PIP4K2A <= median_risk_all, "low", "high")
km_data_all$group<-sfg_new$group
km_data_all <- km_data_all %>%
  mutate(group_new = case_when(
    grepl("Control", group) ~ "Control",
    grepl("AD", group) ~ "AD"
  ))
### boxplot
my_comparisons<- c("Control", "AD")
km_data_all$group_new <- factor(km_data_all$group_new, levels = c("Control", "AD"))

ggboxplot(km_data_all, x = "group_new", y = "PIP4K2A",
          color = "group_new", palette = c("#00AFBB", "#FC4E07","#E7B800"),
          add = "jitter", shape = "group_new",add.params = list(fill = "white"))+ stat_compare_means(method = 't.test')

