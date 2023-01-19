library(randomForest)

# Validation cohort #
fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCTrain = roc(filt_revt$Label, fit.def$votes[,2], plot = T)

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCVal = roc(metaVal$Response, Valresult[,2], plot = T)
auc(ROCVal)
ci.auc(ROCVal)

# Whole cohort #
ROCWhole = roc(c(metaVal$Response, metaTrain$Response), c(Valresult[,2], fit.def$votes[,2]), plot = T)
auc(ROCWhole)
ci.auc(ROCWhole)

# Signature from Iwase, 2021 #
selected_genes = c("APOD", "ASPN", "CCL5", "CD52", "COL2A1", "CXCL11", "CXCL13", "DUSP5", "FOXC1",
                   "GZMB", "HTRA1", "IDO1", "IL23A", "ITM2A", "KMO", "KRT16", "KYNU", "MIA",
                   "PSMB9", "PTGDS", "RARRES3", "RTP4", "S100A8", "SFRP1", "SPTLC2", "TNFAIP8", "TNFSF10")

filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.Iwase <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Iwaseresult = predict(fit.Iwase, filtered_val, type = "prob")
ROCvalIwase = roc(metaVal$Response, Iwaseresult[,2], plot = T, direction = "<")
plot(ROCvalIwase, col = "blue", lwd = 2)
auc(ROCvalIwase)

ROCwholeIwase = roc(c(metaTrain$Response, metaVal$Response), c(fit.Iwase$votes[,2], Iwaseresult[,2]), plot = T, direction = "<")
plot(ROCwholeIwase, col = "red", lwd = 2, main = "Whole cohort")
obs = c(metaTrain$Response, metaVal$Response)
obs[obs == "PD"] = 1
obs[obs == "ORR"] = 0
roc.area(obs = as.numeric(obs), pred = c(fit.Iwase$votes[,2], Iwaseresult[,2]))

roc.test(ROCwhole, ROCwholeIwase)
