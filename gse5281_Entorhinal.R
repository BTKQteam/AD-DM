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

gset <- getGEO("GSE5281", GSEMatrix =TRUE, getGPL=T)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

##clear data matrix
matrix<-as.data.frame(exprs(gset))

#create matrix with gene Gene.Symbol and expression value
features_new <- fData(gset)
features_new <- select(features_new,Gene.Symbol)
matrix <- cbind(features_new,matrix)


#remove duplicate in columns "gene_Gene.Symbol" and NA value
matrix<-matrix[matrix$Gene.Symbol!="",]
matrix$Gene.Symbol <- gsub(" ///.*", "", matrix$Gene.Symbol)

matrix<-matrix[!duplicated(matrix$Gene.Symbol),]
#set "GENE_Gene.Symbol" as index
rownames(matrix) <- matrix$Gene.Symbol
matrix[, "Gene.Symbol"] <- NULL
# log2 transformation
ex <- matrix
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
matrix <- log2(ex) }


##Extract sample Info
sampleInfo <- pData(gset)
sampleInfo<-select(sampleInfo,geo_accession,source_name_ch1,characteristics_ch1.8,characteristics_ch1.9,characteristics_ch1.11,characteristics_ch1.11)
##rename sampleInfo
sampleInfo<-sampleInfo%>%dplyr::rename(ID=geo_accession,tissue=source_name_ch1,age=characteristics_ch1.11,sex=characteristics_ch1.9,phenotype=characteristics_ch1.8)
sampleInfo<-sampleInfo%>%mutate(sex=gsub("sex: ","",sex,fixed = T))
sampleInfo<-sampleInfo%>%mutate(sex=gsub("Sex: ","",sex,fixed = T))
sampleInfo<-sampleInfo%>%mutate(age=gsub("age: ","",age,fixed = T))
sampleInfo<-sampleInfo%>%mutate(age=gsub("Age: ","",age,fixed = T))
sampleInfo<-sampleInfo%>%mutate(age=gsub(" years","",age,fixed = T))
sampleInfo<-sampleInfo%>%mutate(age=gsub(" days","",age,fixed = T))

sampleInfo<-sampleInfo%>%mutate(tissue=gsub("brain, ","",tissue,fixed = T))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("Disease State: ","",phenotype,fixed = T))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("disease state: ","",phenotype,fixed = T))
sampleInfo$sex<-substr(sampleInfo$sex, 1, nchar(sampleInfo$sex)-2)
sampleInfo$phenotype<-substr(sampleInfo$phenotype, 1, nchar(sampleInfo$phenotype)-2)
sampleInfo[] <- lapply(sampleInfo, function(x) gsub(" ", "_", x))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("normal","control",phenotype,fixed = T))
sampleInfo<-sampleInfo%>%mutate(phenotype=gsub("Alzheimer's_Disease","AD",phenotype,fixed = T))
sampleInfo$group<- paste(sampleInfo$phenotype, sampleInfo$tissue, sep="_")
sampleInfo$genderbase<-paste(sampleInfo$group, sampleInfo$sex, sep="_")

##########################
##########################
##########################

#subset 
Entorhinal <-sampleInfo[(sampleInfo$group=='AD_Entorhinal_Cortex'|sampleInfo$group=='control_Entorhinal_Cortex'),]
#create matrix of that pair
mat_Entorhinal <- matrix[,Entorhinal$ID]

#BEFORE EXCLUDE OUTLIER, RUN DEG FIRST:

# Convert the phenotype data to a design matrix
design <- model.matrix(~0+ phenotype,data=Entorhinal)
# Create a linear model and fit it to the gene expression data
fit <- lmFit(mat_Entorhinal, design)
# Create contrasts for the groups of interest
contrast.matrix <- makeContrasts(Entorhinal_vs_control = phenotypeAD - phenotypecontrol, levels=colnames(coef(fit)))

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
       title = "Volcano Plot for AD-Entorhinal vs Control (GSE5281)") +theme_bw()+
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
gene_matrix <- as.matrix(mat_Entorhinal[ids_of_interest,])
annot<-data.frame(group=Entorhinal$group,row.names = rownames(Entorhinal))
pheatmap(gene_matrix,
         labels_row = rownames(gene_matrix),
         scale="row",annotation_col = annot,
         show_colnames = FALSE)
## PCA analysis
#MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(mat_Entorhinal))

## Join the PCs to the sample information
cbind(Entorhinal, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=phenotype,label=paste("", ID))) + geom_point(size=0.1) + geom_text_repel(aes(label = ID), size = 2, nudge_x = 1, nudge_y = 1)+theme_bw()+stat_ellipse()+geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+ theme(panel.border = element_blank())
#write.csv(full_results,file = '5281_Entorhinal_pval.csv')

####WGCNA analysis
library('WGCNA')
library('gridExtra')
library("CorLevelPlot")

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers, and transpose the matrix of genes as columns
power <- c(c(1:20), seq(from = 22, to = 30, by = 2))
W_matrix <-t(mat_hippocampus)
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
  geom_text(nudge_y = 0.1) +  geom_vline(xintercept = 9, color = 'red') +

  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


print(a1+a2)
# convert matrix to numeric
W_matrix[] <- sapply(W_matrix, as.numeric)

soft_power <- 9
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
traits <- hippocampus %>% 
  mutate(AD = ifelse(grepl('AD_Hippocampus', group), 1, 0)) %>% select(8) #(8) means the column "AD" is the 8th columns of 'hippocampus'

# binarize categorical variables


control <-hippocampus %>% 
  mutate(control = ifelse(grepl('control_hippocampus', group), 1, 0)) %>% select(8)
traits <- cbind(traits, control)


# Define numbers of genes and samples
nSamples <- nrow(W_matrix)
nGenes <- ncol(W_matrix)  


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[44:45],
             y = names(heatmap.data)[1:43],
             col = c("blue1", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)
brown<-module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()
write.csv(brown,file = 'brown-hippo-5281.csv',row.names = F)


#####################################
#### IMMUNE INFILTRATION ############
#####################################

signature_score_calculation_methods
names(signature_tme)[1:20]
names(signature_metabolism)[1:20]

mat_hippocampus[1:5, 1:10]



#Analyze all collected signature scores (integrating three methods: PCA, ssGSEA and z-score).

sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = mat_hippocampus,
                             signature       = signature_collection,
                             method          = "ssgsea",
                             mini_gene_count = 2)
##The signature gene sets derived from GO, KEGG, HALLMARK and REACTOME datasets.
sig_hallmark<-calculate_sig_score(pdata           = NULL,
                                  eset            = mat_hippocampus,
                                  signature       = hallmark,
                                  method          = "ssgsea",
                                  mini_gene_count = 2)


###METHOD 1: CIBERSORT

cibersort<-deconvo_tme(eset = mat_hippocampus, method = "cibersort", arrays = FALSE, perm = 200 )
res<-cell_bar_plot(input = cibersort[1:11,], title = "CIBERSORT Cell Fraction of AD-hippocampus")
?cell_bar_plot

##METHOD 2: ESTIMATE
estimate<-deconvo_tme(eset = mat_hippocampus, method = "estimate")
res_es<-cell_bar_plot(input = estimate[1:11,], title = "ESTIMATE Cell Fraction")


####5.2 Phenotype Module
hippocampus$ID<-NULL
hippocampus <- tibble::rownames_to_column(hippocampus, var = "ID")

#unique(sampleInfo$group)

#match<-sam
##id1 and id2 must match!
res<-iobr_cor_plot(pdata_group           = hippocampus,
                   id1                   = "ID",
                   feature_data          = estimate,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "group",
                   is_target_continuous  = FALSE,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   signature_group       = (sig_group)[25],
                   ProjectID             = "hippo-cibersort",
                   
                   palette_box           = "set2",
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

###testdata
pdata_group<-rownames_to_column(as.data.frame(t(mat_hippocampus)),var = "ID")
pdata_group<-as.data.frame(pdata_group[,c("ID","AGFG1","LHX2","AMER2","INPP5D","PMEPA1","SEC14L1","ZNF516")])
head(pdata_group)
sig_meta<-calculate_sig_score(pdata           = NULL,
                              eset            = mat_hippocampus,
                              signature       = signature_metabolism,
                              method          = "pca",
                              mini_gene_count = 2)

res_meta<-iobr_cor_plot(pdata_group           = hippocampus,
                   id1                   = "ID",
                   feature_data          = sig_meta,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "group",
                   is_target_continuous  = F,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   signature_group       = (sig_group)[c(28,30,31)],
                   ProjectID             = "AD-metabolism",
                   
                   palette_box           = "jco",
                   palette_corplot      = "pheatmap",
                   palette_heatmap       = 2,
                   feature_limit         = 200,
                   character_limit       = 50,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE)
write.csv(sig_meta,file = 'res_meta_hippo.csv')
write.csv(res_meta,file = 'sig_meta_score_hippo.csv')

##### AUC ###########
#####################
library(pROC)
risk_list<-read.csv('intersection_AD_DM_new.csv',header = F)
risk_list<-unlist(risk_list)
mat_AUC<-t(mat_Entorhinal[risk_list,])

mat_AUC<-as.data.frame(mat_AUC)
mat_AUC <- mat_AUC[, !colSums(is.na(mat_AUC)) > 0]

label<-as.factor(ifelse(Entorhinal$phenotype=='AD',1,0))


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
write.csv(auc_scores,file = '5281_auc_entorhinal_new.csv')

##plot specific genes
gene <- mat_AUC[, 97]
roc_obj <- roc(label, gene)
plot(roc_obj, add = TRUE,print.auc=T,plot=T,smooth=T,col='red4')
lines(c(0, 1), c(1, 0), col = "gray", lty = "dashed",add=T)


write.csv(auc_scores,file = 'auc_5281_validate_378gene.csv')

mat_AD<-cbind(mat_AUC,label)
#######################################################
############  construct the DMAD signature ############
#######################################################
mat_final_transpose<-t(mat_hippocampus)
selected_genes<-c('STOML1','ASB16','HDAC7','SCAMP2','OS9','ZFP36L1','F5','PARP9','FAM20B','BMP2K','MLEC','BACE2')
km_data_AD <- mat_AD[, c(selected_genes, "label")]
final_AD_matrix<-mat_AUC[,selected_genes]
# Step 1: Calculate combined expression value for each patient
patient_scores <- rowSums(final_AD_matrix)
# Step 2: Calculate average expression for each gene in the signature
average_expression <- colMeans(final_AD_matrix)

# Step 3: Calculate sum of average expressions
sum_average_expression <- sum(average_expression)
# Step 4: Calculate normalized score for each patient
normalized_scores <- patient_scores / sum_average_expression
median_risk<-median(normalized_scores)

km_data_AD$risk <- normalized_scores

km_data_AD$risk <- ifelse(km_data_AD$risk <= median_risk, "low score", "high score")
#write.csv(km_data_AD,file = 'kmdata_AD.csv')


##########################
### TEST METABOLISM ######
km_data_AD<-rownames_to_column(km_data_AD,var = "ID")

final_AD_matrix<-t(final_AD_matrix)
sig_meta<-calculate_sig_score(pdata           = NULL,
                              eset            = mat_hippocampus,
                              signature       = signature_metabolism,
                              method          = "pca",
                              mini_gene_count = 2)

sig_res<-calculate_sig_score(pdata           = NULL,
                              eset            = mat_hippocampus,
                              signature       = signature_collection,
                              method          = "pca",
                              mini_gene_count = 2)
res_meta<-iobr_cor_plot(pdata_group           = km_data_AD,
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
                        ProjectID             = "DMAD-metabolism",
                        
                        palette_box           = "jco",
                        palette_corplot      = "pheatmap",
                        palette_heatmap       = 2,
                        feature_limit         = 200,
                        character_limit       = 50,
                        show_heatmap_col_name = FALSE,
                        show_col              = FALSE,
                        show_plot             = TRUE)



#write.csv(sig_meta,file = 'res_meta_hippo.csv')
write.csv(res_meta,file = 'res_meta_score_panel.csv')

############################
############################
### ML model of 12 gene ####
############################
library(glmnet)
library(pROC)
library(caret)
ML_data<-km_data_AD
ML_data$ID<-NULL
ML_data$risk<-NULL
ML_data$label<-NULL
# Step 2: Split the data into training and testing sets
set.seed(123)
train_indices <- sample(1:nrow(ML_data), 0.7 * nrow(ML_data))
train_data <- ML_data[train_indices, ]
train_labels <- label[train_indices]
test_data <- ML_data[-train_indices, ]
test_labels <- label[-train_indices]


# Step 3: Perform LASSO regression
lasso_model <- cv.glmnet(as.matrix(train_data), as.factor(train_labels), family = "binomial",alpha=1)
optimal_lambda <- lasso_model$lambda.min
lasso_coefficients <- coef(lasso_model, s = optimal_lambda)
plot(lasso_model)
# Extract the log lambda values from the LASSO model
log_lambda <- log(lasso_model$lambda)
# Plot the coefficients
plot(lasso_model$glmnet.fit, 
     "lambda", label=FALSE)

# Step 4: Calculate diagnostic scores
diagnostic_scores <- rowSums(train_data * lasso_coefficients[-1]) + lasso_coefficients[1]

# Step 5: Calculate ROC curve and AUC
roc_obj <- roc(train_labels, diagnostic_scores)
auc_value <- auc(roc_obj)

# Step 6: Calculate confusion matrix and performance metrics
predictions <- ifelse(diagnostic_scores >= 0.5, "AD", "control")
# Convert predictions and test_labels to factors
predictions <- factor(predictions)
train_labels <- factor(train_labels)

# Align factor levels
levels(predictions) <- levels(train_labels)

# Compute confusion matrix
confusion_matrix <- confusionMatrix(predictions, train_labels)
specificity <- confusion_matrix$byClass["Specificity"]
accuracy <- confusion_matrix$overall["Accuracy"]
recall <- confusion_matrix$byClass["Sensitivity"]

# Print the results
print(paste("Specificity:", specificity))
print(paste("Accuracy:", accuracy))
print(paste("Recall:", recall))
print(paste("AUC:", auc_value))

plot(roc_obj,,print.auc=T,plot=T,smooth=T,col='blue4')
lines(c(0, 1), c(1, 0), col = "gray", lty = "dashed",add=T)



######lay toan bo samples:
selected_genes<-c('PIP4K2A')
mat_all_transpose<-t(mat_Entorhinal)
mat_all<-as.data.frame(mat_all_transpose[,selected_genes])
km_data_all<-data.frame(mat_all)
colnames(km_data_all)<-'PIP4K2A'
median_risk_all<-median(km_data_all$PIP4K2A)
km_data_all$risk <- ifelse(km_data_all$PIP4K2A <= median_risk_all, "low", "high")
km_data_all$group<-Entorhinal$phenotype

### boxplot
my_comparisons<- c("control", "AD")
km_data_all$group <- factor(km_data_all$group, levels = c("control", "AD"))

ggboxplot(km_data_all, x = "group", y = "PIP4K2A",
          color = "group", palette = c("#00AFBB", "#FC4E07","#E7B800"),
          add = "jitter", shape = "group",add.params = list(fill = "white"))+ stat_compare_means(method = 't.test')
