set.seed(50)
require('KnowSeq')
require(caret)
require(dplyr)
require(ROCR)
require(limma)
require("gplots")
require("RColorBrewer")

# READING DATA

# GSE156063 California San Francisco
# NP/OP swab
# 234 MUESTRAS ( no virus 100 / other virus 41 / SC2 93)
# INDEX outlier 50 81 113 153 157 171 223 ( 1 SC2 / 4 VARI / 2 NVARI )
# 15959 genes  227 muestras ( no virus 98 / other virus 37 / SC2 92)
# 135 NEG / 92 POS
# EDAD : 20-40 (25%) / 40-50 (25%) / 50-65 (25%) / 65-89 (25%)
# SEXO : F (54%) - M (46%)

GSE156063 <- read.csv('GSE156063_summary.csv' ,header=TRUE)
rownames(GSE156063) <- GSE156063[,1]
GSE156063[,1] <- NULL
GSE156063_m <- t(GSE156063[,1:15959])
GSE156063_lab <- GSE156063$GSE156063_lab
GSE156063_age <- GSE156063$GSE156063_age
GSE156063_gender <- GSE156063$GSE156063_gender
GSE156063_pcr <- GSE156063$GSE156063_pcr
GSE156063_rpm <- GSE156063$GSE156063_rpm
for (i in 1:length(GSE156063_lab)){
  if (GSE156063_lab[i]=='no virus'){
    GSE156063_lab[i] <- 'NON.VIRAL.ARI'
  } else if (GSE156063_lab[i]=='other virus'){
    GSE156063_lab[i] <- 'VIRAL.ARI'
  }
}



# GSE149273 The University of Arizona 
# NP/OP swab
# 90 muestras ( 30 Control / 30 RVA / 30 RVC )
# INDEX outlier 50 51 80 88 (1 Control / 2 RVA / 1 RVC)
# 17605 genes  86 muestras ( 29 Control / 28 RVA / 29 RVC )
# 30 PACIENTES / 3 MUESTAS POR PACIENTE
# SEXO: F (51%) - M(49%)
# AGE: 33-36 
GSE149273 <- read.csv('GSE149273_summary.csv' ,header=TRUE)
rownames(GSE149273) <- GSE149273[,1]
GSE149273[,1] <- NULL
GSE149273_m <- t(GSE149273[,1:17605])
GSE149273_lab <- GSE149273$GSE149273_lab
GSE149273_age <- GSE149273$GSE149273_age
GSE149273_gender <- GSE149273$GSE149273_gender
for (i in 1:length(GSE149273_lab)){
  if (GSE149273_lab[i]=='RVA' | GSE149273_lab[i]=='RVC' ){
    GSE149273_lab[i] <- 'VIRAL.ARI'
  }
}


# GSE152075 University of Washington
# NP swab
# 484 muestras ( Control 54 / SC2 430) 
# INDEX outlier 7  26  63  72  75  82  98 105 118 123 127 132 133 141 159 164 176 182 193 195 198 219 226 228 237 243 249 285 291 297 302 312
# 315 317 332 345 346 347 355 357 361 365 384 395 398 403 408 422 424 430 432 458 471 477 (4 Control / 50 SC2)
# 35486 genes 430 ( Control 50 / SC2 380) 
# AGE: UNKNOWN (3%) / 0-20 (3%) / 20-40 (20%) / 40-60 (35%) / 60-80 (28%) / +80 (11%)
# GENDER : UNKNOWN (11%) / M (41%) / F (48%)
GSE152075 <- read.csv('GSE152075_summary.csv' ,header=TRUE)
rownames(GSE152075) <- GSE152075[,1]
GSE152075[,1] <- NULL
GSE152075_m <- t(GSE152075[,1:35486])
GSE152075_lab <- GSE152075$GSE152075_lab
GSE152075_age <- GSE152075$GSE152075_age
GSE152075_gender <- GSE152075$GSE152075_gender
GSE152075_viral_load <- GSE152075$GSE152075_viral_load


# GSE163151 University of California, San Francisco
# NP swab 
# 351 muestras ( 138 covid / 11 control / NON VIRAL ARI 82 / VIRAL ARI 120)
# INDEX outlier  21  29  35  41  51  61  64  74 312 316 321 334 343 344 (1 Control / 13 VARI)
# 21909 genes 337 (138 covid / 10 control / NON VIRAL ARI 82 / VIRAL ARI 107 ) 
# Viral ARI : Adenovirus / Coronavirus 229E /  Coronavirus HKU1 / Coronavirus NL63 / Coronavirus OC43 / Human metapneumovirus / Influenza / Influenza A /  Influenza / Parainfluenza / Respiratory syncytial virus / Rhinovirus
GSE163151 <- read.csv('GSE163151_summary.csv' ,header=TRUE)
rownames(GSE163151) <- GSE163151[,1]
GSE163151[,1] <- NULL
GSE163151_m <- t(GSE163151[,1:21909])
GSE163151_lab <- GSE163151$GSE163151_lab
GSE163151_patogeno <- GSE163151$GSE163151_patogeno
for (i in 1:length(GSE163151_lab)){
  if (GSE163151_lab[i]=='COVID-19'){
    GSE163151_lab[i] <- 'SC2'
  } else if (GSE163151_lab[i]=='Donor control'){
    GSE163151_lab[i] <- 'Control'
  } else if (GSE163151_lab[i]=='Non-viral acute respiratory illness'){
    GSE163151_lab[i] <- 'NON.VIRAL.ARI'
  } else if (GSE163151_lab[i]=='Viral acute respiratory illness'){
    GSE163151_lab[i] <- 'VIRAL.ARI'
  }
}



# GSE162835 Al jalila Children's Specialty Hospital DUBAI
# NP SWAB 
# 50 muestras SC2
# INDEX outlier
# 20601 genes
# Severity Asymptomatic 37(74%) / Moderate 10(20%) / Severe 3 (6%)
# Sexo: F 18 (36%) / M 32 (64%)
# Age : 0-20 (6%) / 20-40 (50%) / 40-60 (28%) / 60-80 (16%)
GSE162835 <- read.csv('GSE162835_summary.csv' ,header=TRUE)
rownames(GSE162835) <- GSE162835[,1]
GSE162835[,1] <- NULL
GSE162835_m <- t(GSE162835[,1:20601])
GSE162835_lab <- GSE162835$GSE162835_labels
GSE162835_age <- GSE162835$GSE162835_age
GSE162835_gender <- GSE162835$GSE162835_gender
GSE162835_severity <- GSE162835$GSE162835_severity



# Superserie 5 bacth
I1 <- intersect(rownames(GSE156063_m), rownames(GSE149273_m))
I2 <- intersect(I1, rownames(GSE152075_m))
I3 <- intersect(I2,  rownames(GSE163151_m))
I4 <- intersect(I3, rownames(GSE162835_m))
GSE156063_I <- GSE156063_m[which(rownames(GSE156063_m) %in%  I4),]
GSE149273_I <- GSE149273_m[which(rownames(GSE149273_m) %in%  I4),]
GSE152075_I <- GSE152075_m[which(rownames(GSE152075_m) %in%  I4),]
GSE163151_I <- GSE163151_m[which(rownames(GSE163151_m) %in%  I4),]
GSE162835_I <- GSE162835_m[which(rownames(GSE162835_m) %in%  I4),]
GSE156063_I <- GSE156063_I[order(rownames(GSE156063_I)),] 
GSE149273_I <- GSE149273_I[order(rownames(GSE149273_I)),] 
GSE152075_I <- GSE152075_I[order(rownames(GSE152075_I)),] 
GSE163151_I <- GSE163151_I[order(rownames(GSE163151_I)),] 
labels <- c(GSE156063_lab,GSE149273_lab,GSE152075_lab,GSE163151_lab,GSE162835_lab)
expression_matrix <- cbind(GSE156063_I,GSE149273_I,GSE152075_I,GSE163151_I,GSE162835_I)


# Normalization
expression_matrix_norm_scale <- normalizeBetweenArrays(expression_matrix, method = 'scale')
# efecto batch 
expression_matrix_norm_scale_fix <- batchEffectRemoval(expression_matrix_norm_scale,labels, method = 'sva')
# outliers 
outliers_scale <- RNAseqQA(expression_matrix_norm_scale_fix,toRemoval = TRUE, toPNG = FALSE, toPDF = FALSE) #5
# elimino outliers
expression_matrix_norm_scale_fix_out <- expression_matrix_norm_scale_fix[,-which(colnames(expression_matrix_norm_scale_fix) %in% outliers_scale$outliers)]
labels_scale <- labels[-which(colnames(expression_matrix_norm_scale_fix) %in% outliers_scale$outliers)]

# TRAIN-TEST
Index_train_test_scale <- createDataPartition(labels_scale, p = .80, 
                                              list = FALSE, 
                                              times = 1)
train_labels_scale <- labels_scale[Index_train_test_scale]
test_labels_scale <- labels_scale[-Index_train_test_scale]
train_matrix_scale <- expression_matrix_norm_scale_fix_out[,Index_train_test_scale]
test_matrix_scale <- expression_matrix_norm_scale_fix_out[,-Index_train_test_scale]

#limma
DEGs_scale2 <- DEGsExtraction(train_matrix_scale, as.factor(train_labels_scale), lfc=1, cov=2, pvalue = 0.01, number = Inf, CV=TRUE)
#FS algorithms
gene_mrmr_scale2 <- featureSelection(t(train_matrix_scale),train_labels_scale,DEGs_scale2$Common_DEGs, mode ='mrmr')
gene_rf_scale2 <- featureSelection(t(train_matrix_scale),train_labels_scale,DEGs_scale2$Common_DEGs, mode ='rf')
gene_da_scale2 <- featureSelection(t(train_matrix_scale),train_labels_scale,DEGs_scale2$Common_DEGs, mode ='da',disease = 'COVID-19')


#knn
knn_train_mrmr_scale2 <- knn_trn(t(train_matrix_scale), as.factor(train_labels_scale), names(gene_mrmr_scale2), 5) 
knn_test_mrmr_scale2 <- knn_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), names(gene_mrmr_scale2), bestK =  knn_train_mrmr_scale2$bestK)

knn_train_rf_scale2 <- knn_trn(t(train_matrix_scale), as.factor(train_labels_scale), gene_rf_scale2, 5) 
knn_test_rf_scale2 <- knn_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), gene_rf_scale2, bestK =  knn_train_rf_scale2$bestK)

knn_train_da_scale2 <- knn_trn(t(train_matrix_scale), as.factor(train_labels_scale), names(gene_da_scale2), 5) 
knn_test_da_scale2 <- knn_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), names(gene_da_scale2), bestK =  knn_train_da_scale2$bestK)

knn_train_limma_scale2 <- knn_trn(t(train_matrix_scale), as.factor(train_labels_scale), gene_limma_scale2, 5) 
knn_test_limma_scale2 <- knn_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), gene_limma_scale2, bestK =  knn_train_limma_scale2$bestK)


#Validation plot
plot(x,knn_train_mrmr_scale2$accuracyInfo$meanAccuracy, type = 'l', col= 'black', ylab='Metric Performance', xlab='Genes', lwd=2, ylim = c(0.8,1), panel.first = grid(col='gray45'))
lines(x,knn_train_mrmr_scale2$sensitivityInfo$meanSensitivity, col='blue', lwd=2, lty=2)
lines(x,knn_train_mrmr_scale2$specificityInfo$meanSpecificity, col='#FF8B00', lwd=2, lty=4)
lines(x,knn_train_mrmr_scale2$F1Info$meanF1, col='red', lwd=2, lty=4)
legend(x=14.305 ,y =0.8345, c('Accuracy', 'Sensitivity','Specificity','F1-Score'), lty = c(1,2,4,5), col = c('black','blue','#FF8B00','red'))

#svm
svm_train_mrmr_scale2 <- svm_trn(t(train_matrix_scale), as.factor(train_labels_scale), names(gene_mrmr_scale2), 5) 
svm_test_mrmr_scale2 <- svm_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), names(gene_mrmr_scale2), bestParameters =  svm_train_mrmr_scale2$bestParameters)

svm_train_rf_scale2 <- svm_trn(t(train_matrix_scale), as.factor(train_labels_scale), gene_rf_scale2, 5) 
svm_test_rf_scale2 <- svm_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), gene_rf_scale2, bestParameters =  svm_train_rf_scale2$bestParameters)

svm_train_da_scale2 <- svm_trn(t(train_matrix_scale), as.factor(train_labels_scale), names(gene_da_scale2), 5) 
svm_test_da_scale2 <- svm_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), names(gene_da_scale2), bestParameters =  svm_train_da_scale2$bestParameters)

svm_train_limma_scale2 <- svm_trn(t(train_matrix_scale), as.factor(train_labels_scale), gene_limma_scale2, 5) 
svm_test_limma_scale2 <- svm_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), gene_limma_scale2, bestParameters =  svm_train_limma_scale2$bestParameters)

#rf
rf_train_mrmr_scale2 <- rf_trn(t(train_matrix_scale), as.factor(train_labels_scale), names(gene_mrmr_scale2), 5) 
rf_test_mrmr_scale2 <- rf_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), names(gene_mrmr_scale2), bestParameters =  rf_train_mrmr_scale2$bestParameters)

rf_train_rf_scale2 <- rf_trn(t(train_matrix_scale), as.factor(train_labels_scale), gene_rf_scale2, 5) 
rf_test_rf_scale2 <- rf_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), gene_rf_scale2, bestParameters =  rf_train_rf_scale2$bestParameters)

rf_train_da_scale2 <- rf_trn(t(train_matrix_scale), as.factor(train_labels_scale), names(gene_da_scale2), 5) 
rf_test_da_scale2 <- rf_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), names(gene_da_scale2), bestParameters =  rf_train_da_scale2$bestParameters)

rf_train_limma_scale2 <- rf_trn(t(train_matrix_scale), as.factor(train_labels_scale), gene_limma_scale2, 5) 
rf_test_limma_scale2 <- rf_test(t(train_matrix_scale), as.factor(train_labels_scale), t(test_matrix_scale), as.factor(test_labels_scale), gene_limma_scale2, bestParameters =  rf_train_limma_scale2$bestParameters)

#confusion matrix plot
dataPlot(svm_test_mrmr_scale2$cfMats[[15]]$table, labels = test_labels_scale  ,mode = "confusionMatrix",toPNG = FALSE, toPDF = FALSE)

# AUC plot
response <- as.factor(test_labels_scale)
aucs <- rep(NA, length(levels(response))) # store AUCs
legendLabels <- as.character()
colours <- c('red','blue','green','black')
plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
     ylab="Sensitivity",
     xlab="1 - Specificity",
     bty='n')

for (i in seq_along(levels(response))) {
  cur.class <- levels(response)[i]
  binaryTraining.labels <- as.factor(train_labels_scale == cur.class)
  binaryTest.labels <- as.factor(test_labels_scale == cur.class)
  
  binary_svm_cv_mrmr_results <- svm_trn(t(train_matrix_scale), binaryTraining.labels, names(gene_mrmr_scale2)[1:15],5)
  
  binary_svm_test_mrmr_results <- svm_test(t(train_matrix_scale), binaryTraining.labels, t(test_matrix_scale), binaryTest.labels, names(gene_mrmr_scale2)[1:15], bestParameters = svm_train_mrmr_scale2$bestParameters)
  
  score <- binary_svm_test_mrmr_results$predictions[[15]]
  score <- as.vector(score)
  score[score=='FALSE'] <- 0
  score[score=='TRUE'] <- 1
  binaryTest.labels <- as.vector(binaryTest.labels)
  binaryTest.labels[binaryTest.labels=='FALSE'] <- 0
  binaryTest.labels[binaryTest.labels=='TRUE'] <- 1
  pred <- prediction(as.numeric(score), as.numeric(binaryTest.labels))
  perf <- performance(pred, "tpr", "fpr")
  roc.x <- unlist(perf@x.values)
  roc.y <- unlist(perf@y.values)
  lines(roc.y ~ roc.x, col = colours[i], lwd = 2)
  # store AUC
  auc <- performance(pred, "auc")
  auc <- unlist(slot(auc, "y.values"))
  aucs[i] <- auc
  legendLabels[i] <- paste(levels(response)[i], " AUC: ",format(round(aucs[i], 4), nsmall = 3),sep = "")
}

print(paste0("Mean AUC under the precision-recall curve is: ", round(mean(aucs), 2)))

lines(x=c(0,1), c(0,1))
legend(x=0.5 ,y =0.3, legendLabels, lty=1, ncol= 1,inset = c(0,0),  col = colours, cex = 1.3,lwd=3)




# Clustering & heatmap
condition_colors <- unlist(lapply(train_labels_scale,function(x){
  if(grepl("NON.VIRAL.ARI",x)) 'blue' #pink
  else if (grepl('Control',x)) 'red' #grey
  else if (grepl('SC2',x)) 'green'
  else if (grepl('VIRAL.ARI',x)) 'black'
}))

# I like to have a line just to assign input in one step
input <- as.matrix(train_matrix_scale[which(rownames(train_matrix_scale) %in% names(gene_mrmr_scale2)),])


par(oma = c(5, 1, 0, 1))
heatmap.2(input, trace="none", density="none", col=bluered(20), cexRow=1, cexCol=0.2,
          ColSideColors=condition_colors, scale="row",
          hclust=function(x) hclust(x))
par(fig = c(0, 1, 0, 1), oma = c(20, 0, 0, 40), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
labelegend <- c('Non viral ARI','Control','SC2','Viral ARI')
legend("bottom", labelegend, ncol= 1,inset = .02, fill = c('blue','red','green','black'), cex = 1.1)
          
#t-SNE clustering
#BiocManager::install("M3C")
require(M3C)
tsne(expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out)%in%names(gene_mrmr_scale2)[1:15]),],labels=as.factor(labels_scale),controlscale=TRUE, scale=3, colvec = c('red','blue','green','black'))


# boxplot 
order_genes <- expression_matrix_norm_scale_fix_out[which(rownames(expression_matrix_norm_scale_fix_out) %in% names(gene_mrmr_scale2[1:15])),]
dataPlot(order_genes[c(7,9,8,2,15,10,11,4,13,5,12,6,1,3,14),],labels_scale,mode = "genesBoxplot",toPNG = FALSE, colours = c("darkred", "forestgreen", "darkorchid3", "dodgerblue3"),colnumber = 3)















