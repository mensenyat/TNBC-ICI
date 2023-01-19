library(dplyr)
library(M3C)
library(varSelRF)
library(tidyr)
library(pROC)
library(ggplot2)
library(verification)

load("~/Environment after ML.RData")

pdf("Pan-cancer ROC.pdf")

listName = c("GSE111636", "GSE165252", "GSE181815", "GSE183924", "GSE67501", "GSE78220", "IMvigor210", "E-MTAB-3218", "E-MTAB-4030")

# We have to repeat the RF each time because not all genes might be present in all datasets #

#### GSE111636 ####
dataGSE = read.table("GSE111636/GSE111636 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("GSE111636/GSE111636 metadata.csv", sep = ",", header = T)

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCR.simple_Best1 = list(predictions=fit.def$votes[,2], labels=filt_revt$Label)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[1], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

#### GSE165252 ####
dataGSE = read.table("GSE165252/GSE165252 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("GSE165252/GSE165252 metadata.csv", sep = ",", header = T)
metaGSE = metaGSE[metaGSE$Time == "baseline",]
metaGSE$Response[metaGSE$Response == "responder"] = "ORR"
metaGSE$Response[metaGSE$Response == "non-responder"] = "PD"

dataGSE = dataGSE[,colnames(dataGSE) %in% metaGSE$Sample_description]
metaGSE = metaGSE[order(match(metaGSE$Sample_description, colnames(dataGSE))),]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCR.simple_Best1 = list(predictions=fit.def$votes[,2], labels=filt_revt$Label)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[2], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

levels(metaGSE$Response) <- c(0, 1)
metaGSE$Response[metaGSE$Response == "ORR"] = 0
metaGSE$Response[metaGSE$Response == "PD"] = 1
ra <- roc.area(as.numeric(as.vector(metaGSE$Response)), GSEresult[,2])
ra$p.value

#### GSE181815 ####
dataGSE = read.table("GSE181815/GSE181815 data.csv", sep = ",", header = T)
rownames(dataGSE) = make.names(dataGSE$GeneSymbol, unique = T)
dataGSE = dataGSE[,-1]

colnames(dataGSE) = substr(colnames(dataGSE), 1, 7)
metaGSE = read.table("GSE181815/GSE181815 metadata.txt", sep = "\t", header = T)
metaGSE$characteristics..Response.to.pembrolizumab[metaGSE$characteristics..Response.to.pembrolizumab %in% 
                                                     c("Progressive Disease")] = "PD"
metaGSE$characteristics..Response.to.pembrolizumab[metaGSE$characteristics..Response.to.pembrolizumab %in% 
                                                     c("Partial Response", "Complete Response", "Stable Disease")] = "ORR"
metaGSE$Sample.name = gsub("-", ".", metaGSE$Sample.name)
metaGSE = metaGSE[metaGSE$Sample.name %in% colnames(dataGSE),]
metaGSE = metaGSE[order(match(metaGSE$Sample.name, colnames(dataGSE))),]

dataGSE = dataGSE[,colnames(dataGSE) %in% metaGSE$Sample.name]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCR.simple_Best1 = list(predictions=fit.def$votes[,2], labels=filt_revt$Label)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
ROCGSE = roc(metaGSE$characteristics..Response.to.pembrolizumab, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[3], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

#### GSE183924 ####
dataGSE = read.table("GSE183924/GSE183924 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("GSE183924/GSE183924 metadata.txt", sep = "\t", header = T)

metaGSE$Response = NA
metaGSE$Response[metaGSE$Relapse == "Yes"] = "PD"
metaGSE$Response[metaGSE$Relapse == "No"] = "ORR"
metaGSE = metaGSE[!is.na(metaGSE$Response),]

dataGSE = dataGSE[,colnames(dataGSE) %in% metaGSE$RNA.Seq.ID]
metaGSE = metaGSE[order(match(metaGSE$RNA.Seq.ID, colnames(dataGSE))),]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCR.simple_Best1 = list(predictions=fit.def$votes[,2], labels=filt_revt$Label)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
rownames(GSEresult) == metaGSE$RNA.Seq.ID
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[4], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

levels(metaGSE$Response) <- c(0, 1)
metaGSE$Response[metaGSE$Response == "ORR"] = 0
metaGSE$Response[metaGSE$Response == "PD"] = 1
ra <- roc.area(as.numeric(as.vector(metaGSE$Response)), GSEresult[,2])
ra$p.value

#### GSE67501 ####
dataGSE = read.table("GSE67501/GSE67501 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("GSE67501/GSE67501 metadata.csv", sep = ",", header = T)

dataGSE = dataGSE[,colnames(dataGSE) %in% metaGSE$Sample_geo_accession]
metaGSE = metaGSE[order(match(metaGSE$Sample_geo_accession, colnames(dataGSE))),]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
rownames(GSEresult) == metaGSE$Sample_geo_accession
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[5], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

#### GSE78220 ####
dataGSE = read.table("GSE78220/GSE78220 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("GSE78220/GSE78220 metadata.csv", sep = ",", header = T)
metaGSE$Response[metaGSE$Response == "Progressive Disease"] = "PD"
metaGSE$Response[metaGSE$Response %in% c("Partial Response", "Complete Response", "Stable Disease")] = "ORR"

dataGSE = dataGSE[,colnames(dataGSE) %in% metaGSE$Patient]
metaGSE = metaGSE[order(match(metaGSE$Patient, colnames(dataGSE))),]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
rownames(GSEresult) == metaGSE$Patient
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[6], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

levels(metaGSE$Response) <- c(0, 1)
metaGSE$Response[metaGSE$Response == "ORR"] = 0
metaGSE$Response[metaGSE$Response == "PD"] = 1
ra <- roc.area(as.numeric(as.vector(metaGSE$Response)), GSEresult[,2])
ra$p.value

#### IMvigor210 ####
dataGSE = read.table("IMvigor210/IMvigor210 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("IMvigor210/IMvigor210 metadata.txt", sep = "\t", header = T)
metaGSE = metaGSE[!metaGSE$Best.Confirmed.Overall.Response == "NE",]
metaGSE$Best.Confirmed.Overall.Response[metaGSE$Best.Confirmed.Overall.Response %in% c("PR", "CR", "SD")] = "ORR"
metaGSE$Response = metaGSE$Best.Confirmed.Overall.Response

dataGSE = dataGSE[,colnames(dataGSE) %in% rownames(metaGSE)]
metaGSE = metaGSE[order(match(rownames(metaGSE), colnames(dataGSE))),]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
rownames(GSEresult) == rownames(metaGSE)
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[7], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

levels(metaGSE$Response) <- c(0, 1)
metaGSE$Response[metaGSE$Response == "ORR"] = 0
metaGSE$Response[metaGSE$Response == "PD"] = 1
ra <- roc.area(as.numeric(as.vector(metaGSE$Response)), GSEresult[,2])
ra$p.value

#### E-MTAB-3218 ####
dataGSE = read.table("E-MTAB-3218/E-MTAB-3218 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("E-MTAB-3218/E-MTAB-3218 metadata.txt", sep = "\t", header = T)
metaGSE = metaGSE[metaGSE$Characteristics.biopsy.timepoint. == "Screen",]
metaGSE$Source.Name = paste("X", metaGSE$Source.Name, sep = "")
metaGSE = metaGSE[!metaGSE$Characteristics.BOR. == "NE",]

metaGSE$Response = NA
metaGSE$Response[metaGSE$Characteristics.BOR3. %in% c("CRPR", "SD")] = "ORR"
metaGSE$Response[metaGSE$Characteristics.BOR3. == "PD"] = "PD"

dataGSE = dataGSE[,colnames(dataGSE) %in% metaGSE$Source.Name]
metaGSE = metaGSE[order(match(metaGSE$Source.Name, colnames(dataGSE))),]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
rownames(GSEresult) == metaGSE$Source.Name
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[8], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

levels(metaGSE$Response) <- c(0, 1)
metaGSE$Response[metaGSE$Response == "ORR"] = 0
metaGSE$Response[metaGSE$Response == "PD"] = 1
ra <- roc.area(as.numeric(as.vector(metaGSE$Response)), GSEresult[,2])
ra$p.value

#### E-MTAB-4030 ####
dataGSE = read.table("E-MTAB-4030/E-MTAB-4030 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE = read.table("E-MTAB-4030/E-MTAB-4030 metadata.txt", sep = "\t", header = T)
metaGSE = metaGSE[metaGSE$Characteristics.sampling.time.point. == "screen",]
metaGSE$Source.Name = paste("X", metaGSE$Source.Name, sep = "")
metaGSE = metaGSE[!metaGSE$Characteristics.BestOverallResponse.RECIST. == "NonEvaluable",]

metaGSE$Response = NA
metaGSE$Response[metaGSE$Characteristics.BestOverallResponse.RECIST. %in% 
                   c("Complete/Partial Response", "Stable Disease")] = "ORR"
metaGSE$Response[metaGSE$Characteristics.BestOverallResponse.RECIST. == "Progressive Disease"] = "PD"

dataGSE = dataGSE[,colnames(dataGSE) %in% metaGSE$Source.Name]
metaGSE = metaGSE[order(match(metaGSE$Source.Name, colnames(dataGSE))),]

selected_genes = Genes1$Gene[as.numeric(Genes1$Freq) >= 585]
selected_genes = selected_genes[selected_genes %in% rownames(dataGSE)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCtrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<")
auc(ROCtrain)
ci.auc(ROCtrain)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCval = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<")
auc(ROCval)
ci.auc(ROCval)

dataGSEt = as.data.frame(as.matrix(t(dataGSE)))
dataGSEt = as.data.frame(scale(dataGSEt))
filtered_GSE = dataGSEt[,which(colnames(dataGSEt) %in% selected_genes)]
filtered_GSE = filtered_GSE[,order(match(colnames(filtered_GSE), colnames(filt_revt)))]
colnames(filtered_GSE) == colnames(filt_revt)

GSEresult = predict(fit.def, filtered_GSE, type = "prob")
rownames(GSEresult) == metaGSE$Source.Name
ROCGSE = roc(metaGSE$Response, GSEresult[,2], plot = F, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .001), boot.n=500, conf.level = 0.95)
plot(ROCGSE, main = paste(listName[9], " AUC = ", round(ci.auc(ROCGSE)[2],2), 
                          " (", round(ci.auc(ROCGSE)[1],2), "-", round(ci.auc(ROCGSE)[3],2), ")", sep = "")) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))

levels(metaGSE$Response) <- c(0, 1)
metaGSE$Response[metaGSE$Response == "ORR"] = 0
metaGSE$Response[metaGSE$Response == "PD"] = 1
ra <- roc.area(as.numeric(as.vector(metaGSE$Response)), GSEresult[,2])
ra$p.value


dev.off()
