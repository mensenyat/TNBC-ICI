# Installation steps #
install.packages("pROC")
install.packages("splitTools")
install.packages("tidyverse")
install.packages("varSelRF")
install.packages("verification")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("M3C")
BiocManager::install("sva")

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

load("~/Initial data.RData")

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

# THIS PART IS RESOURCE AND TIME CONSUMING - DON'T DO IT IF YOU ONLY WANT TO CHECK THE CODE, OR CHANGE THE NUMBER OF ITERATIONS #
N = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
AUC = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
AUCval = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
AUCvalSum = as.data.frame(matrix(NA, ncol = 31, nrow = 1000))
SelGenes = list(NA)

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
# THIS PART IS RESOURCE AND TIME CONSUMING - DON'T DO IT IF YOU ONLY WANT TO CHECK THE CODE, OR CHANGE THE NUMBER OF ITERATIONS #
N2 = as.data.frame(matrix(NA, ncol = length(unique(Genes1$Freq)), nrow = 30))
AUC2 = as.data.frame(matrix(NA, ncol = length(unique(Genes1$Freq)), nrow = 30))
AUCval2 = as.data.frame(matrix(NA, ncol = length(unique(Genes1$Freq)), nrow = 30))

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
