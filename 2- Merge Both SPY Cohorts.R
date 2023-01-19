library(dplyr)
library(M3C)
library(tidyr)
library(ROCR)
library(ggplot2)
library(Rtsne)
library(sva)

Ola = read.table("GSE173839_data.txt", sep = "\t", header = T, row.names = 1)
genes = read.table("Probes_to_genes_updated.txt", sep = "\t", header = T)
Ola = merge(genes, Ola, by.x = 1, by.y = 0)
rownames(Ola) <- make.names(Ola$GeneName, unique = TRUE)
Ola = Ola[,-c(1:4)]
Ola = data.frame(apply(Ola, 2, function(x) as.numeric(as.character(x))), row.names = rownames(Ola))

metaOla = read.table("Metadata.txt", sep = "\t", header = T)
metaOla = metaOla[which(metaOla$Hormone == "TNBC"),]
metaOla = metaOla[which(metaOla$X.Sample_characteristics_ch1.4 == "arm: durvalumab/olaparib"),]

Ola = Ola[,which(colnames(Ola) %in% metaOla$Patient)]

setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature")
Pembro = read.table("GeneExp_SPY2_Neoadjuvant_GPL20078_Corrected.txt", sep = "\t", header = T, row.names = 1)
metaPembro = read.table("Clinical_Data_GPL20078.txt", sep = "\t", header = T)
metaPembro = metaPembro[which(metaPembro$HR == 0 & metaPembro$HER2 == 0),]
metaPembro = metaPembro[which(metaPembro$Treatment == "Paclitaxel+Pembrolizumab"),]
rownames(metaPembro) = paste("X",metaPembro$Patient.ID, sep = "")

Pembro = Pembro[,order(match(colnames(Pembro), rownames(metaPembro)))]
Pembro = Pembro[,which(colnames(Pembro) %in% rownames(metaPembro))]
Pembro = data.frame(apply(Pembro, 2, function(x) as.numeric(as.character(x))), row.names = rownames(Pembro))

# Select only genes in both #
Ola = Ola[which(rownames(Ola) %in% rownames(Pembro)),]
Pembro = Pembro[which(rownames(Pembro) %in% rownames(Ola)),]
Ola = Ola[order(match(rownames(Ola), rownames(Pembro))),]
Pembrot = as.data.frame(as.matrix(t(Pembro)))
Olat = as.data.frame(as.matrix(t(Ola)))

#### Correct batch effect ####
# Merge all datasets #
spy = merge(Pembro, Pembro, by=0)
rownames(spy) = spy[,1]
spy = spy[,-1]
spy = na.omit(spy)

cb.df.mdata <- cbind.data.frame("sample" = colnames(spy), # exclude uid column, c(1)
                                "batch" = c(rep("SPY-Ola", ncol(spy)), rep("SPY-Pembro", ncol(Pembro))))

cb.corr.model <- model.matrix(~1, data = cb.df.mdata)

data_corrected = ComBat(dat=spy,
                        batch=cb.df.mdata$batch,
                        mod=cb.corr.model,
                        par.prior=TRUE,
                        prior.plot=FALSE)
data_corrected = na.omit(data_corrected)

tumorLabels=as.character(cb.df.mdata$batch)
plottingColors=c("Red", "cadetblue")
names(plottingColors)=unique(tumorLabels)

Tsne = Rtsne(as.matrix(t(data_corrected)),perplexity=10)
pdf('tSNE Both SPY cohorts Corrected.pdf')
plot(Tsne$Y, main="Corrected Methylation", col=plottingColors[tumorLabels], pch=19)
legend('bottomleft', legend=unique(tumorLabels), fill=plottingColors[unique(tumorLabels)], border=T, title='Subtype')
dev.off()

Tsne = Rtsne(as.matrix(t(spy)),perplexity=10)
pdf('tSNE Both SPY cohorts uncorrected.pdf')
plot(Tsne$Y, main="Uncorrected Methylation", col=plottingColors[tumorLabels], pch=19)
legend('bottomleft', legend=unique(tumorLabels), fill=plottingColors[unique(tumorLabels)], border=T, title='Subtype')
dev.off()

write.table(data_corrected, file = "Both SPY corrected.txt", sep = "\t")

#### Merge clinical data ####
setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature")
metaOla = read.table("Metadata.txt", sep = "\t", header = T)
metaOla = metaOla[which(metaOla$Hormone == "TNBC"),]
metaOla = metaOla[which(metaOla$X.Sample_characteristics_ch1.4 == "arm: durvalumab/olaparib"),]

metaBoth = as.data.frame(cbind("Patient" = c(metaSPY$Patient, rownames(meta)),
                           "Response" = c(metaSPY$Response, meta$Response),
                           "Batch" = c(rep("Olaparib", nrow(metaSPY)), rep("Pembro", nrow(meta)))))
metaBoth = metaBoth[order(match(metaBoth$Patient, colnames(data_corrected))),]

write.table(metaBoth, file = "Meta_both_SPY_cohorts.txt", sep = "\t", row.names = F)
