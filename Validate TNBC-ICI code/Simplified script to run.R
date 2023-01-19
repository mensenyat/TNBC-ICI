# Installation steps #
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("M3C")

install.packages("pROC")
install.packages("splitTools")
install.packages("tidyverse")
install.packages("varSelRF")

if (!require("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("emilelatour/lamisc")

# Start of the script #
library(lamisc)
library(varSelRF)
library(tidyverse)
library(pROC)
library(M3C)
library(splitTools)

#setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/Both_Cohorts")
#setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer")
#score = read.table("Score immunotherapy pan-cancer.txt", sep = "\t")
#score = score[order(-abs(score$finalScore)),]
#scoreHigh = score[1:500,]
#
#setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/")
#load("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/Both_Cohorts/Environment Both Cohorts 1000RF v2 - after both cohorts.RData")
#rm(list=setdiff(ls(), c("metaTrain", "metaVal")))
#
#data = read.table("Both SPY corrected.txt", sep = "\t", header = T)
#
#meta = read.table("Meta_both_SPY_cohorts.txt", sep = "\t", header = T)
#meta$Patient = meta$Patient
#
#data = data[rownames(data) %in% rownames(scoreHigh),]
#datat = as.data.frame(t(data))

#### FROM HERE ####
load("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/GitHub/Simplified to run/Initial data.RData")

scoreHigh = score[1:500,] # Previously calculated pan-cancer score. Select the top 500 genes with the highest score.
# To see how the pan-cancer score was calculated, go to "Identify pan-cancer ICI score.R" 

data = data[rownames(data) %in% rownames(scoreHigh),]
datat = as.data.frame(t(data))
datat = scale(datat)
datat = as.data.frame(as.matrix(t(na.omit(as.matrix(t(datat)))))) # Remove genes with NA

dataTrain = as.data.frame(datat[which(rownames(datat) %in% rownames(metaTrain)),])
dataTraint = as.data.frame(as.matrix(t(dataTrain)))
dataTraint = as.data.frame(na.omit(dataTraint))
dataTrain = as.data.frame(as.matrix(t(dataTraint)))

metaORR = meta[which(meta$Response == "ORR"),]
metaPD = meta[which(meta$Response == "PD"),]

dataORR = dataTraint[which(colnames(dataTraint) %in% metaORR$Patient)]
dataPD = dataTraint[which(colnames(dataTraint) %in% metaPD$Patient)]

# First, we will select only the genes that follow the same trend that we observed pan-cancer #
meanORR = apply(dataORR, 1, mean)
meanPD = apply(dataPD, 1, mean)

fold = (meanORR-meanPD)
foldScore = merge(score, fold, by=0)
foldScore = foldScore[foldScore$finalScore < 0 & foldScore$y < 0 | foldScore$finalScore > 0 & foldScore$y > 0,]

# Repeat only with genes following the trend #
data = data[rownames(data) %in% foldScore$Row.names,]
datat = as.data.frame(t(data))
datat = datat
datat = scale(datat)
datat = as.data.frame(as.matrix(t(na.omit(as.matrix(t(datat)))))) # Remove genes with NA
dataTrain = as.data.frame(datat[which(rownames(datat) %in% rownames(metaTrain)),])

dataTraint = as.data.frame(as.matrix(t(dataTrain)))
dataTraint = as.data.frame(na.omit(dataTraint))
dataTrain = as.data.frame(as.matrix(t(dataTraint)))

# Create dataframe with Subtype data #
tumorMetadata = as.data.frame(cbind("barcode" = as.character(metaTrain$Patient), "Response" = as.character(metaTrain$Response)))
tumorLabels = as.character(tumorMetadata$Response)
plottingColors = c("Red", "cadetblue")
names(plottingColors) = unique(tumorLabels)

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
dataValt = as.data.frame(as.matrix(t(dataVal)))


setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/Both_Cohorts")
N = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
AUC = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
AUCval = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
AUCvalSum = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
SelGenes = list(NA)

# THIS PART IS RESOURCE AND TIME CONSUMING - DON'T DO IT IF YOU ONLY WANT TO CHECK THE CODE, OR CHANGE THE NUMBER OF ITERATIONS #
Cl_Genes = for(l in 1:1000) { # This result may vary between different iterations of the same code due to random forest being RANDOM
  # Random Forest #
  RF_Cl = varSelRF(dataTrain, as.factor(tumorLabels), c.sd = 2, mtryFactor = 1, ntree = 1000,
                   ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2,
                   whole.range = TRUE, recompute.var.imp = TRUE, verbose = FALSE,
                   returnFirstForest = FALSE, fitted.rf = NULL, keep.forest = FALSE)
   
  sel_hist = RF_Cl$selec.history
  sel_hist2=data.frame(x = sel_hist$Vars.in.Forest)
  All_variables = sel_hist2 %>% separate(x, as.character(1:ncol(dataTrain)))
  SelGenes[[l]] = All_variables  
  
  for(i in 1:(nrow(All_variables)-6)){ # Only signatures with >10 genes
    selected_genes = All_variables[i,]
    filt_rev = dataTraint[which(rownames(dataTraint) %in% selected_genes),]
    filt_revt = as.data.frame(as.matrix(t(filt_rev)))
    filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
    filt_revt$Label = metaTrain$Response
    
    fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
    ROCTrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<", quiet = T)
    auc(ROCTrain)
    
    N[l,i] = length(selected_genes[which(!is.na(selected_genes))])
    AUC[l,i] = auc(ROCTrain)
    
    dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
    dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
    filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]
    
    Valresult = predict(fit.def, filtered_val, type = "prob")
    ROCVal = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<", quiet = T)
    
    AUCval[l,i] = auc(ROCVal)
  }
  print(l)
}

#AUCmean = apply(AUC, 2, mean, na.rm=TRUE)
#AUCsd = apply(AUC, 2, sd, na.rm=TRUE)
#AUCvalmean = apply(AUCval, 2, mean, na.rm=TRUE)
#AUCvalsd = apply(AUCval, 2, sd, na.rm=TRUE)
#
#
#CL = as.data.frame(cbind("AUC" = unlist(AUC), "N" = unlist(N), 
#                         "AUCsd" = rep(AUCsd, each = 3), "AUCmean" = rep(AUCmean, each = 3),
#                         "AUCval" = unlist(AUCval), "AUCvalmean" = rep(AUCvalmean, each = 3), 
#                         "AUCValsd" = rep(AUCvalsd, each = 3)))
#CL = na.omit(CL)
#
#p = ggplot(CL, aes(x=N, y=AUCvalmean)) + 
#  geom_line(data = CL, aes(x = N, y = AUCvalmean), color = "blue") +
#  geom_line(data = CL, aes(x = N, y = AUCmean), color = "red") +
#  geom_point(data = CL, aes(x = N, y = AUCvalmean), color = "blue")+ 
#  geom_point(data = CL, aes(x = N, y = AUCmean), color = "red")+ 
#  geom_errorbar(aes(ymin=AUCvalmean-AUCvalsd, ymax=AUCvalmean+AUCvalsd), width=.5, color = "blue") +
#  geom_errorbar(aes(ymin=AUCmean-AUCsd, ymax=AUCmean+AUCsd), width=.5, color = "red") +
#  xlab('N') +
#  ylab('AUC') +
#  xlim(10,150) +
#  ylim(0.5,1) +
#  scale_color_brewer(palette="Paired")+theme_classic() 
#p 

Sign = as.data.frame(matrix(NA, ncol = ncol(SelGenes[[1]]), nrow = nrow(AUCval)))

for(i in 1:nrow(AUCvalSum)){
  AUCvalSum[i,] = AUCval[i,]# +AUC[i,]
  Sign[i,] = as.character(SelGenes[[i]][which.max(AUCvalSum[i,]),])
}

rep1 = count_duplicates(Sign)
rep1 = rep1[order(-rep1$dupe_count),]

#### Best combination per number of repetitions ####
AllGenes = na.omit(unlist(rep1[,-ncol(rep1)]))

Genes1 = as.data.frame(table(AllGenes))
Genes1$Freq = as.numeric(Genes1$Freq)
Genes1 = Genes1[order(-Genes1$Freq),]

# Compute AUC with the selected combinations with different lengths #
N2 = as.data.frame(matrix(NA, ncol = length(unique(Genes1$Freq)), nrow = 30))
AUC2 = as.data.frame(matrix(NA, ncol = length(unique(Genes1$Freq)), nrow = 30))
AUCval2 = as.data.frame(matrix(NA, ncol = length(unique(Genes1$Freq)), nrow = 30))

# THIS PART IS RESOURCE AND TIME CONSUMING - DON'T DO IT IF YOU ONLY WANT TO CHECK THE CODE, OR CHANGE THE NUMBER OF ITERATIONS #
for(l in 1:30){
  for(i in 1:length(unique(Genes1$Freq))){ # We perform the analysis removing "one freq" at a time
    selected_genes = Genes1$AllGenes[as.numeric(Genes1$Freq) >= unique(Genes1$Freq)[i]]
    filt_rev = dataTraint[which(rownames(dataTraint) %in% selected_genes),]
    filt_revt = as.data.frame(as.matrix(t(filt_rev)))
    filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
    filt_revt$Label = metaTrain$Response
    
    fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
    ROCTrain = roc(filt_revt$Label, fit.def$votes[,2], plot = F, direction = "<", quiet = T)
    auc(ROCTrain)
    
    N2[l,i] = length(selected_genes[which(!is.na(selected_genes))])
    AUC2[l,i] = auc(ROCTrain)
    
    dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
    dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
    filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]
    
    Valresult = predict(fit.def, filtered_val, type = "prob")
    ROCR.simple_Best1 = list(predictions=Valresult[,2], labels=metaVal$Response)
    ROCVal = roc(metaVal$Response, Valresult[,2], plot = F, direction = "<", quiet = T)
    auc(ROCVal)
    
    AUCval2[l,i] = auc(ROCVal)
  }
  print(l)
}

AUCmean = apply(AUC2, 2, mean, na.rm=TRUE)
AUCsd = apply(AUC2, 2, sd, na.rm=TRUE)
AUCval2mean = apply(AUCval2, 2, mean, na.rm=TRUE)
AUCval2sd = apply(AUCval2, 2, sd, na.rm=TRUE)

CL = as.data.frame(cbind("N" = unlist(N2), "AUCsd" = rep(AUCsd, each = 30), "AUCmean" = rep(AUCmean, each = 30),
                         "AUCval2mean" = rep(AUCval2mean, each = 30), "AUCval2sd" = rep(AUCval2sd, each = 30)))
CL = na.omit(CL)

p = ggplot(CL, aes(x=N, y=CL$AUCval2mean)) + 
  geom_line(data = CL, aes(x = N, y = AUCval2mean), color = "blue") +
  geom_line(data = CL, aes(x = N, y = AUCmean), color = "red") +
  geom_point(data = CL, aes(x = N, y = AUCval2mean), color = "blue") + 
  geom_point(data = CL, aes(x = N, y = AUCmean), color = "red") + 
  geom_errorbar(aes(ymin=AUCval2mean-AUCval2sd, ymax=AUCval2mean+AUCval2sd), width=.1, color = "blue") +
  geom_errorbar(aes(ymin=AUCmean-AUCsd, ymax=AUCmean+AUCsd), width=.1, color = "red") +
  xlab('N') +
  ylab('AUC') +
  #xlim(11,130) +
  ylim(0.5,1) +
  scale_color_brewer(palette="Paired") + theme_classic() 
p 

#### Plot ROC Curves ####
selected_genes = Genes1$AllGenes[as.numeric(Genes1$Freq) >= 585] # This is our selected cutoff. It may diverge in each case
selected_genes = c("ADAMTS12", "BATF2", "CBWD5", "C6orf89", "CD8A", "CDH6", "COL5A3", "CXCL13", "CXCL9", "FGD3", "GALNT5",
                   "GBP1", "GBP4", "GUCY1A2", "GZMB", "IGDCC4", "ITGA2", "ITGB3", "KLHDC7B", "LDB2", "MED9", "MGLL", "NLRC5",
                   "PHLDB1", "PRSS1", "PTPN2", "RABEP1", "RBBP8", "RPA1", "RPA4", "RTN3", "SALL1", "SCN4B", "THEG",
                   "TMEM92", "TRIM37", "KIAA0100") # If you want to replicate the plots from the study, use this list of genes
                                                   # However, slightly different results are expected due to the random component of the algorithm

filt_rev = dataTraint[which(rownames(dataTraint) %in% selected_genes),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCTrain = roc(filt_revt$Label, fit.def$votes[,2], plot = T, direction = "<", quiet = T)
auc(ROCTrain)

pdf("ROC Curve Training.pdf")
ci.sp.obj <- ci.sp(ROCTrain, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCTrain) 
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))
dev.off()

dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCVal = roc(metaVal$Response, Valresult[,2], plot = T)
auc(ROCVal)
ci.auc(ROCVal)

pdf("ROC Curve Validation.pdf")
ci.sp.obj <- ci.sp(ROCVal, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCVal) 
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))
dev.off()

# Whole cohort #
ROCWhole = roc(c(metaVal$Response, metaTrain$Response), c(Valresult[,2], fit.def$votes[,2]), plot = T)
auc(ROCWhole)
ci.auc(ROCWhole)

pdf("ROC Curve Whole cohort.pdf")
ci.sp.obj <- ci.sp(ROCWhole, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCWhole) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))
dev.off()






#### Analyze in SPY non-TNBC and plot ROC curves ####
OlaHR = read.table("GSE173839_data.txt", sep = "\t", header = T, row.names = 1)
genes = read.table("Probes_to_genes_Updated.txt", sep = "\t", header = T)
OlaHR = merge(genes, OlaHR, by.x = 1, by.y = 0)
rownames(OlaHR) <- make.names(OlaHR$GeneName, unique = TRUE)
OlaHR = OlaHR[,-c(1:7)]

metaOlaHR = read.table("SPY/Metadata.txt", sep = "\t", header = T)
metaOlaHR = metaOlaHR[which(metaOlaHR$Hormone == "HR"),]
metaOlaHR = metaOlaHR[which(metaOlaHR$X.Sample_characteristics_ch1.4 == "arm: durvalumab/olaparib"),]

OlaHR = OlaHR[,which(colnames(OlaHR) %in% metaOlaHR$Patient)]
OlaHRt = as.data.frame(as.matrix(t(OlaHR)))
OlaHRt = data.frame(apply(OlaHRt, 2, function(x) as.numeric(as.character(x))), row.names = rownames(OlaHRt))
OlaHRt = as.data.frame(scale(OlaHRt))

metaOlaHR = metaOlaHR[metaOlaHR$Patient %in% colnames(OlaHR),]
metaOlaHR = metaOlaHR[order(match(metaOlaHR$Patient, rownames(OlaHRt))),]
rownames(metaOlaHR) = metaOlaHR$Patient

# Pembro cohort #
dataPemHR = read.table("GeneExp_SPY2_Neoadjuvant_GPL20078_Corrected.txt", sep = "\t", header = T, row.names = 1)
metaPemHR = read.table("Clinical_Data_GPL20078.txt", sep = "\t", header = T)
metaPemHR = metaPemHR[which(metaPemHR$HR == 1 & metaPemHR$HER2 == 0),]
metaPemHR = metaPemHR[which(metaPemHR$Treatment == "Paclitaxel+Pembrolizumab"),]
rownames(metaPemHR) = paste("X",metaPemHR$Patient.ID, sep = "")

dataPemHR = dataPemHR[,which(colnames(dataPemHR) %in% rownames(metaPemHR))]
dataPemHRt = as.data.frame(as.matrix(t(dataPemHR)))
dataPemHRt = as.data.frame(scale(dataPemHRt))
metaPemHR = metaPemHR[order(match(rownames(metaPemHR), rownames(dataPemHRt))),]

# Validation #
dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]
filtered_val = dataVal[,which(colnames(dataVal) %in% selected_genes)]

Valresult = predict(fit.def, filtered_val, type = "prob")
ROCVal = roc(metaVal$Response, Valresult[,2], plot = T)
auc(ROCVal)
ci.auc(ROCVal)

pdf("ROC Curve Validation.pdf")
ci.sp.obj <- ci.sp(ROCVal, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCVal) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))
dev.off()

# Whole cohort #
ROCWhole = roc(c(metaVal$Response, metaTrain$Response), c(Valresult[,2], fit.def$votes[,2]), plot = T)
auc(ROCWhole)
ci.auc(ROCWhole)

pdf("ROC Curve Whole cohort.pdf")
ci.sp.obj <- ci.sp(ROCWhole, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCWhole) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))
dev.off()

# Ola #
filtered_val = OlaHRt[,which(colnames(OlaHRt) %in% selected_genes)]
Olaresult = predict(fit.def, filtered_val, type = "prob")
ROCOla = roc(metaOlaHR$Response, Olaresult[,2], plot = T, direction = "<")
auc(ROCOla)
ci.auc(ROCOla)

pdf("ROC Curve Olaparib HR.pdf")
ci.sp.obj <- ci.sp(ROCOla, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCOla) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))
dev.off()

# Pembro #
filtered_val = dataPemHRt[,which(colnames(dataPemHRt) %in% selected_genes)]
GSEresult = predict(fit.def, filtered_val, type = "prob")
ROCGSE = roc(metaPemHR$Response, GSEresult[,2], plot = T, direction = "<")
auc(ROCGSE)
ci.auc(ROCGSE)

pdf("ROC Curve Pembrolizumab HR.pdf")
ci.sp.obj <- ci.sp(ROCGSE, sensitivities=seq(0, 1, .01), boot.n=500, conf.level = 0.95)
plot(ROCGSE) # restart a new plot
plot(ci.sp.obj, type="shape", col=alpha("grey", 0.1))
dev.off()


#### UMAP - Volcano - Heatmap ####
filt_cl = datat[,colnames(datat) %in% selected_genes]

#pdf("UMAP 37 genes Whole Cohort.pdf")
umap(as.matrix(t(filt_cl)),labels=as.factor(meta$Response),controlscale=TRUE,scale=3)
dev.off()

# Volcano Plot
#pdf("Volcano Plot ORR vs PD.pdf")
plot(Zratio, -log10(pvalue), main = "ORR vs PD - Volcano", xlim=c(-3, 3), ylim=c(0,3), xlab="fold")
points (Zratio[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "red")
points (Zratio[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "blue")

abline(v = fold_cutoff, col = "red", lwd = 3)
abline(v = -fold_cutoff, col = "blue", lwd = 3)
abline(h = -log10(qvalue_cutoff), col = "grey", lwd = 3)
dev.off()


#### Identify DEG using all genes ####
setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro")
data = read.table("Both SPY corrected.txt", sep = "\t", header = T)
datat = as.data.frame(t(data))
#spy = data.frame(apply(spy, 2, function(x) as.numeric(as.character(x))), row.names = rownames(spy))

meta = read.table("Meta_both_SPY_cohorts.txt", sep = "\t", header = T)
meta$Patient = meta$Patient
setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/SPY-Pembro/Both_Cohorts")

metaTrain = meta[part$train,]
metaVal = meta[part$valid,]
table(metaVal$Batch)
table(metaTrain$Batch)

datat = scale(datat)
datat = as.data.frame(as.matrix(t(na.omit(as.matrix(t(datat)))))) # Remove genes with NA

dataTrain = as.data.frame(datat[which(rownames(datat) %in% rownames(metaTrain)),])
dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])

dataTraint = as.data.frame(as.matrix(t(dataTrain)))
dataTraint = as.data.frame(na.omit(dataTraint))

dataValt = as.data.frame(as.matrix(t(dataVal)))

dataTrain = as.data.frame(as.matrix(t(dataTraint)))

metaORR = meta[which(meta$Response == "ORR"),]
metaPD = meta[which(meta$Response == "PD"),]

dataORR = dataTraint[which(colnames(dataTraint) %in% metaORR$Patient)]
dataPD = dataTraint[which(colnames(dataTraint) %in% metaPD$Patient)]

meanORR = apply(dataORR, 1, mean)
meanPD = apply(dataPD, 1, mean)

fold = (meanORR-meanPD)
Zratio = (fold/sd(fold))

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(dataORR)) {
  x = dataORR[i,]
  y = dataPD[i,]
  
  t = t.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}

qvalue = p.adjust(pvalue, method = "fdr", n = length(pvalue)) # Correction by False Discovery Rate

combined  = cbind(qvalue, fold, abs(fold), dataTraint)
combined_t = as.data.frame(as.matrix(t(combined)))

# Filter per Z-Ratio and qvalue #
fold_cutoff = 1.5
qvalue_cutoff = 0.05

# fold-change filter for "biological" significance
filter_by_fold = abs(Zratio) >= fold_cutoff
#dim(dataTrain[filter_by_fold, ]) 

# P-value filter for "statistical" significance
filter_by_qvalue = pvalue <= qvalue_cutoff
#dim(dataTrain[filter_by_qvalue, ]) 

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_qvalue 

filtered = dataTraint[filter_combined,]
dim(filtered)
filtered_combined = merge(combined, filtered, by=0)
Probes_with_changes = filtered_combined[,1:4]
UpProbes = Probes_with_changes$Row.names[which(Probes_with_changes$fold>0)]
DownProbes = Probes_with_changes$Row.names[which(Probes_with_changes$fold<0)]

# Create dataframe with Subtype data #
tumorMetadata = cbind("barcode" = as.character(meta$Patient), "Response" = as.character(meta$Response))
filtered_t = as.data.frame(as.matrix(t(filtered)))
filtered_metadata = merge(tumorMetadata, filtered_t, by.x=1, by.y=0)
rownames(filtered_metadata) = filtered_metadata[,1]
names(filtered_metadata) = c(as.character(colnames(tumorMetadata)), as.character(colnames(filtered_t)))
tumorLabels = as.character(filtered_metadata$Response)
plottingColors = c("Red", "cadetblue")
names(plottingColors) = unique(tumorLabels)

filtered_metadata_t = as.data.frame(as.matrix(t(filtered_metadata)))
filtered_Clean = filtered_metadata[,-c(1,2)]
filtered_Clean_t = as.data.frame(as.matrix(t(filtered_Clean)))
filtered_RF = filtered_Clean
filtered_RF$Label = filtered_metadata$Response

filtered_val = dataVal

#pdf("UMAP All DEG Training Cohort.pdf")
umap(filtered_Clean_t,labels=as.factor(metaTrain$Response),controlscale=TRUE,scale=3)
dev.off()


#pdf("UMAP All DEG Validation Cohort.pdf")
umap(dataValt,labels=as.factor(metaVal$Response),controlscale=TRUE,scale=3)
dev.off()

# Volcano Plot
#pdf("Volcano Plot ORR vs PD All Genes.pdf")
plot(Zratio, -log10(pvalue), main = "ORR vs PD - Volcano", xlim=c(-4, 4), ylim=c(0,4), xlab="fold")
points (Zratio[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "red")
points (Zratio[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "blue")

abline(v = fold_cutoff, col = "red", lwd = 3)
abline(v = -fold_cutoff, col = "blue", lwd = 3)
abline(h = -log10(qvalue_cutoff), col = "grey", lwd = 3)
dev.off()

#### Heatmap Training - Validation - Whole Cohort ####
# Training #
tumorMetadata = cbind("barcode" = as.character(meta$Patient), "Response" = as.character(meta$Response))
filtered_t = as.data.frame(as.matrix(t(filtered)))
filtered_metadata = merge(tumorMetadata, filtered_t, by.x=1, by.y=0)
rownames(filtered_metadata) = filtered_metadata[,1]
names(filtered_metadata) = c(as.character(colnames(tumorMetadata)), as.character(colnames(filtered_t)))
tumorLabels = as.character(filtered_metadata$Response)
plottingColors = c("Red", "cadetblue")
names(plottingColors) = unique(tumorLabels)

filtered_metadata_t = as.data.frame(as.matrix(t(filtered_metadata)))
filtered_Clean = filtered_metadata[,-c(1,2)]
filtered_Clean_t = as.data.frame(as.matrix(t(filtered_Clean)))
filtered_RF = filtered_Clean
filtered_RF$Label = filtered_metadata$Response

filtered_metadata$Color[filtered_metadata$Response=="PD"]="red"
filtered_metadata$Color[filtered_metadata$Response=="ORR"] = "yellow"
colv = as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(filtered_Clean))))))
rowv = as.dendrogram(hclust(as.dist(1-cor(as.matrix(filtered_Clean)))))

library(gplots)
pdf("Heatmap_All_DEG_Training.pdf")
heatmap.2(as.matrix(t(filtered_Clean)), cexCol=0.7,
          col = rev(redblue(256)), Rowv=rowv, Colv=colv, scale = "row", trace=NULL, tracecol=NULL, ColSideColors = filtered_metadata$Color)
legend('topleft', legend = unique(tumorLabels), fill = unique(filtered_metadata$Color), border=T, title='Subtype')
dev.off()

# Validation #
dataVal = as.data.frame(datat[which(rownames(datat) %in% rownames(metaVal)),])
dataVal = dataVal[order(match(rownames(dataVal), rownames(metaVal))),]

filtered_val = dataVal[,colnames(dataVal)%in%colnames(filtered_RF)]
filtered_metadata = merge(tumorMetadata, filtered_val, by.x=1, by.y=0)
rownames(filtered_metadata) = filtered_metadata[,1]
names(filtered_metadata) = c(as.character(colnames(tumorMetadata)), as.character(colnames(filtered_t)))
tumorLabels = as.character(filtered_metadata$Response)
plottingColors = c("Red", "cadetblue")
names(plottingColors) = unique(tumorLabels)

filtered_metadata_t = as.data.frame(as.matrix(t(filtered_metadata)))
filtered_Clean = filtered_metadata[,-c(1,2)]
filtered_Clean_t = as.data.frame(as.matrix(t(filtered_Clean)))
filtered_RF = filtered_Clean
filtered_RF$Label = filtered_metadata$Response

filtered_metadata$Color[filtered_metadata$Response=="PD"]="red"
filtered_metadata$Color[filtered_metadata$Response=="ORR"] = "yellow"
colv = as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(filtered_Clean))))))
rowv = as.dendrogram(hclust(as.dist(1-cor(as.matrix(filtered_Clean)))))

library(gplots)
pdf("Heatmap_All_DEG_Validation.pdf")
heatmap.2(as.matrix(t(filtered_Clean)), cexCol=0.7,
          col = rev(redblue(256)), Rowv=rowv, Colv=colv, scale = "row", trace=NULL, tracecol=NULL, ColSideColors = filtered_metadata$Color)
legend('topleft', legend = unique(tumorLabels), fill = unique(filtered_metadata$Color), border=T, title='Subtype')
dev.off()

pdf("UMAP All DEG Validation Cohort.pdf")
umap(filtered_Clean_t,labels=as.factor(filtered_metadata$Response),controlscale=TRUE,scale=3)
dev.off()

# Whole Cohort #
filtered_val = datat[,colnames(datat)%in%colnames(filtered_t)]
filtered_val = data.frame(apply(filtered_val, 2, function(x) as.numeric(as.character(x))), row.names = rownames(filtered_val))
filtered_metadata = merge(tumorMetadata, filtered_val, by.x=1, by.y=0)
rownames(filtered_metadata) = filtered_metadata[,1]
names(filtered_metadata) = c(as.character(colnames(tumorMetadata)), as.character(colnames(filtered_t)))
tumorLabels = as.character(filtered_metadata$Response)
plottingColors = c("Red", "cadetblue")
names(plottingColors) = unique(tumorLabels)

filtered_metadata_t = as.data.frame(as.matrix(t(filtered_metadata)))
filtered_Clean = filtered_metadata[,-c(1,2)]
filtered_Clean_t = as.data.frame(as.matrix(t(filtered_Clean)))
filtered_RF = filtered_Clean
filtered_RF$Label = filtered_metadata$Response

filtered_metadata$Color[filtered_metadata$Response=="PD"]="red"
filtered_metadata$Color[filtered_metadata$Response=="ORR"] = "yellow"
colv = as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(filtered_Clean))))))
rowv = as.dendrogram(hclust(as.dist(1-cor(as.matrix(filtered_Clean)))))

library(gplots)
pdf("Heatmap_All_DEG_Whole_Cohort.pdf")
heatmap.2(as.matrix(t(filtered_Clean)), cexCol=0.7,
          col = rev(redblue(256)), Rowv=rowv, Colv=colv, scale = "row", trace=NULL, tracecol=NULL, ColSideColors = filtered_metadata$Color)
legend('topleft', legend = unique(tumorLabels), fill = unique(filtered_metadata$Color), border=T, title='Subtype')
dev.off()

#### TCGA Survival ####
setwd("D:/R/BRCA/DiNome/2- Ethnic Analysis")
clin = read.table("1- Clinical Data BC Summary.csv", header = T, sep = ",")
clin = clin[which(clin$Unified.Subtyping == "4-HRneg/HER2neg"),]
clin = clin[which(clin$CPE.Purity > 0.5),]

setwd("D:/R/BRCA")
data = read.table("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", sep = "\t", header = T)
colnames(data) = substr(colnames(data),1,12)
rownames(data) <- make.names(data[,1], unique = TRUE)
data = data[,-1]
data = data[,colnames(data) %in% clin$Patient.ID]
data = log2(data+1)
dataSD = apply(data, 1, sd)
data = data[which(!dataSD == 0),]
dataMean = apply(data, 1, mean)
data = data[which(dataMean > 0.1),]
datat = as.data.frame(scale(as.matrix(t(data))))
data = as.data.frame(as.matrix(t(datat)))

selected_genes2 = selected_genes[selected_genes %in% rownames(data)]
filt_rev = filtered_Clean_t[which(rownames(filtered_Clean_t) %in% selected_genes2),]
filt_revt = as.data.frame(as.matrix(t(filt_rev)))
filt_revt = filt_revt[order(match(rownames(filt_revt), rownames(metaTrain))),]
filt_revt$Label = metaTrain$Response

fit.def <- randomForest(factor(Label) ~ ., data = filt_revt, importance = TRUE, type = "classification", ntree = 10000)
ROCR.simple_Best1 = list(predictions=fit.def$votes[,2], labels=filt_revt$Label)
pred_Best1 = prediction(ROCR.simple_Best1$predictions, ROCR.simple_Best1$labels)
perf_Best1 = performance(pred_Best1,"tpr","fpr")
perfAUC_Best1 = performance(pred_Best1,"tpr","fpr", measure="auc")
perfAUC_Best1@y.values

dataCl = data[rownames(data) %in% selected_genes,]
dataCl = dataCl[order(match(rownames(dataCl), colnames(filt_revt))),]
dataClt = as.data.frame(as.matrix(t(dataCl)))

GSEresult = predict(fit.def, dataClt, type = "prob")
clinRes = merge(clin, GSEresult, by.x=1, by.y=0)
clinRes$Disease.Free.Status[clinRes$Disease.Free.Status == "Recurred/Progressed"] = 1
clinRes$Disease.Free.Status[clinRes$Disease.Free.Status == "DiseaseFree"] = 0

write.table(clinRes, file = "Clinical TCGA for kmplot.txt", sep = "\t")

#### ROC Curves using each gene ####
dataFilt = datat[,colnames(datat) %in% selected_genes]

for(i in 1:ncol(dataFilt)){
  pdf(paste("ROC per gene/ROC Curve ", colnames(dataFilt)[i], ".pdf", sep = ""))
  ROCgene = roc(meta$Response, dataFilt[,i], plot = F, direction = "<")
  if(auc(ROCgene) > 0.5){
    plot(ROCgene, main = paste(colnames(dataFilt)[i], "AUC =", round(auc(ROCgene), 2)), col = "blue", lwd = 10) # restart a new plot
    par(new=T)
    plot(ROCwhole, col = "black", lwd = 10)
  } else{
    ROCgene = roc(meta$Response, dataFilt[,i], plot = F, direction = ">")
    plot(ROCgene, main = paste(colnames(dataFilt)[i], "AUC =", round(auc(ROCgene), 2)), col = "blue", lwd = 10) # restart a new plot
    par(new=T)
    plot(ROCwhole, col = "black", lwd = 10)
  }
  dev.off()
}
