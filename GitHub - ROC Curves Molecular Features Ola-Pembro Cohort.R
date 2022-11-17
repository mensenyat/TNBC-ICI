load("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/Both_Cohorts/Environment Both Cohorts 1000RF v2 - after both cohorts.RData")
library(randomForest)
library(pROC)
library(verification)

# Compare Training vs Validation in ORR and PD #
molOla = read.table("Clinical_Data.txt", sep = "\t", header = T, row.names = 1)

data = read.table("Both SPY corrected.txt", sep = "\t", header = T)
datat = as.data.frame(t(data))
datat = as.data.frame(scale(datat))

molPem = read.table("Biomarker_signatures.txt", sep = "\t", header = T)
molPem$Patient = paste("X", molPem$Patient, sep = "")

Pat = read.table("Meta_both_SPY_cohorts.txt", sep = "\t", header = T)
molPem = molPem[molPem$Patient %in% Pat$Patient,]

molAll = merge(t(molPem), t(molOla), by = 0)
rownames(molAll) = molAll$Row.names
molAll = molAll[,-1]
molAll = as.data.frame(t(molAll))
rownames(molAll) = molAll$Patient
molAll = molAll[,-6]
molAll = molAll[order(match(rownames(molAll), rownames(datat))),]
molAll = as.data.frame(cbind(molAll, "PDL1" = datat$CD274, "PD1" = datat$PDCD1))
Response = molAll$Response
molAll = molAll[,-6]
molAll = as.data.frame(cbind(molAll, "Response" = Response))


selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = T, direction = "<")
plot(ROCtrain, col = "blue", lwd = 2)
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]
Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = T, direction = "<")
plot(ROCval, col = "blue", lwd = 2)
auc(ROCval)
ci.auc(ROCval)
obs = metaVal$Response
obs[obs == "PD"] = 1
obs[obs == "ORR"] = 0
roc.area(obs = as.numeric(obs), pred = Valresult[,2]) # Method to obtain the p-value

ROCwhole = roc(c(metaTrain$Response, metaVal$Response), c(fit.def$votes[,2], Valresult[,2]), plot = T, direction = "<")
plot(ROCwhole, col = "red", lwd = 2, main = "Whole cohort, AUC = 0.94")
ci.auc(ROCwhole)
obs = c(metaTrain$Response, metaVal$Response)
obs[obs == "PD"] = 1
obs[obs == "ORR"] = 0
roc.area(obs = as.numeric(obs), pred = c(fit.def$votes[,2], Valresult[,2])) 

# With one clinical feature #
setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/Both_Cohorts")

for(i in 1:(ncol(molAll)-1)){
  pdf(paste("Clinical ROCs/", colnames(molAll)[i], " ROC.pdf", sep = ""))
  plot(ROCwhole, col = "blue", lwd = 2)
  par(new=T)

  ROCFeat = roc(molAll$Response, as.numeric(molAll[,i]), plot = T, direction = "<")

  if(auc(ROCFeat) > 0.5){
    plot(ROCFeat, col = "red", lwd = 2, main = paste(colnames(molAll)[i], " AUC=", auc(ROCFeat)))
  } else {
    ROCFeat = roc(molAll$Response, as.numeric(molAll[,i]), plot = T, direction = ">")
    plot(ROCFeat, col = "red", lwd = 2, main = paste(colnames(molAll)[i], " AUC=",   perfAUC_BestClin@y.values))
  }
  dev.off()
}

#### ROC Curves adding molecular signatures ####
molClean = molAll[-ncol(molAll)]

dataMol = merge(datat, molClean, by = 0)
rownames(dataMol) = dataMol$Row.names
dataMol = dataMol[,-1]
dataMol = dataMol[order(match(rownames(dataMol), rownames(meta))),]
dataMolt = as.data.frame(as.matrix(t(dataMol)))
dataMolTrain = dataMolt[,colnames(dataMolt) %in% rownames(metaTrain)]
dataMolVal = dataMolt[,colnames(dataMolt) %in% rownames(metaVal)]

selected_genes = c(Genes1$Gene[as.numeric(Genes1$Freq) >= 585], colnames(molClean))
filt_rev = dataMolTrain[which(rownames(dataMolTrain) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt$Label = metaTrain$Response
fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)

#pdf("ROC Training Signature + MolSign.pdf")
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = T, direction = "<")
dev.off()
auc(ROCtrain)

dataVal = as.data.frame(dataMol[which(rownames(dataMol) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]
Valresult = predict(fit.def, filtered_val, type = "prob")

#pdf("ROC Validation Signature + MolSign.pdf")
ROCval = roc(metaVal$Response, Valresult[,2], plot = T, direction = "<")
dev.off()
auc(ROCval)
obs = metaVal$Response
obs[obs == "PD"] = 1
obs[obs == "ORR"] = 0
roc.area(obs = as.numeric(obs), pred = Valresult[,2])

#pdf("ROC Whole Cohort Signature + MolSign.pdf")
ROCwhole = roc(c(metaTrain$Response, metaVal$Response), c(fit.def$votes[,2], Valresult[,2]), plot = T, direction = "<")
dev.off()
auc(ROCwhole)
obs = c(metaTrain$Response, metaVal$Response)
obs[obs == "PD"] = 1
obs[obs == "ORR"] = 0
roc.area(obs = as.numeric(obs), pred = c(fit.def$votes[,2], Valresult[,2]))

#### Using only Clinical Features ####
dataMol = merge(datat, molClean, by = 0)
rownames(dataMol) = dataMol$Row.names
dataMol = dataMol[,-1]
dataMol = dataMol[order(match(rownames(dataMol), rownames(meta))),]
dataMolt = as.data.frame(as.matrix(t(dataMol)))
dataMolTrain = dataMolt[,colnames(dataMolt) %in% rownames(metaTrain)]
dataMolVal = dataMolt[,colnames(dataMolt) %in% rownames(metaVal)]


selected_genes = colnames(molClean)
filt_rev = dataMolTrain[which(rownames(dataMolTrain) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)

#pdf("ROC Training MolSign.pdf")
ROCtrainMol = roc(filt_revt$Label, fit.def$votes[,2], plot = T, direction = "<")
dev.off()
auc(ROCtrainMol)

dataVal = as.data.frame(dataMol[which(rownames(dataMol) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]
Valresult = predict(fit.def, filtered_val, type = "prob")
#pdf("ROC Validation MolSign.pdf")
ROCvalMol = roc(metaVal$Response, Valresult[,2], plot = T, direction = "<")
dev.off()
auc(ROCvalMol)
obs = metaVal$Response
obs[obs == "PD"] = 1
obs[obs == "ORR"] = 0
roc.area(obs = as.numeric(obs), pred = Valresult[,2])

#pdf("ROC Whole Cohort MolSign.pdf")
ROCwholeMol = roc(c(metaTrain$Response, metaVal$Response), c(fit.def$votes[,2], Valresult[,2]), plot = T, direction = "<")
dev.off()
auc(ROCwholeMol)
obs = c(metaTrain$Response, metaVal$Response)
obs[obs == "PD"] = 1
obs[obs == "ORR"] = 0
roc.area(obs = as.numeric(obs), pred = c(fit.def$votes[,2], Valresult[,2]))
