library(dplyr)
library(M3C)
library(tidyr)
library(ROCR)
library(ggplot2)
library(Rtsne)
library(sva)

Durva = read.table("GSE173839_data.txt", sep = "\t", header = T, row.names = 1)
genes = read.table("Probes_to_genes_updated.txt", sep = "\t", header = T)
Durva = merge(genes, Durva, by.x = 1, by.y = 0)
rownames(Durva) <- make.names(Durva$GeneName, unique = TRUE)
Durva = Durva[,-c(1:4)]
Durva = data.frame(apply(Durva, 2, function(x) as.numeric(as.character(x))), row.names = rownames(Durva))

metaDurva = read.table("Metadata.txt", sep = "\t", header = T)
metaDurva = metaDurva[which(metaDurva$Hormone == "TNBC"),]
metaDurva = metaDurva[which(metaDurva$X.Sample_characteristics_ch1.4 == "arm: durvalumab/olaparib"),]

Durva = Durva[,which(colnames(Durva) %in% metaDurva$Patient)]

Pembro = read.table("GeneExp_SPY2_Neoadjuvant_GPL20078_Corrected.txt", sep = "\t", header = T, row.names = 1)
metaPembro = read.table("Clinical_Data_GPL20078.txt", sep = "\t", header = T)
metaPembro = metaPembro[which(metaPembro$HR == 0 & metaPembro$HER2 == 0),]
metaPembro = metaPembro[which(metaPembro$Treatment == "Paclitaxel+Pembrolizumab"),]
rownames(metaPembro) = paste("X",metaPembro$Patient.ID, sep = "")

Pembro = Pembro[,order(match(colnames(Pembro), rownames(metaPembro)))]
Pembro = Pembro[,which(colnames(Pembro) %in% rownames(metaPembro))]
Pembro = data.frame(apply(Pembro, 2, function(x) as.numeric(as.character(x))), row.names = rownames(Pembro))

# Select only genes in both #
Durva = Durva[which(rownames(Durva) %in% rownames(Pembro)),]
Pembro = Pembro[which(rownames(Pembro) %in% rownames(Durva)),]
Durva = Durva[order(match(rownames(Durva), rownames(Pembro))),]
Pembrot = as.data.frame(as.matrix(t(Pembro)))
Durvat = as.data.frame(as.matrix(t(Durva)))

#### Correct batch effect ####
# Merge all datasets #
spy = merge(Pembro, Pembro, by=0)
rownames(spy) = spy[,1]
spy = spy[,-1]
spy = na.omit(spy)

cb.df.mdata <- cbind.data.frame("sample" = colnames(spy), # exclude uid column, c(1)
                                "batch" = c(rep("SPY-Durva", ncol(spy)), rep("SPY-Pembro", ncol(Pembro))))

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
metaDurva = read.table("Metadata SPY-Durvalumab.txt", sep = "\t", header = T)
metaDurva = metaDurva[which(metaDurva$Hormone == "TNBC"),]
metaDurva = metaDurva[which(metaDurva$X.Sample_characteristics_ch1.4 == "arm: durvalumab/olaparib"),]

metaBoth = as.data.frame(cbind("Patient" = c(metaSPY$Patient, rownames(meta)),
                               "Response" = c(metaSPY$Response, meta$Response),
                               "Batch" = c(rep("Durva", nrow(metaSPY)), rep("Pembro", nrow(meta)))))
metaBoth = metaBoth[order(match(metaBoth$Patient, colnames(data_corrected))),]

write.table(metaBoth, file = "Meta_both_SPY_cohorts.txt", sep = "\t", row.names = F)
