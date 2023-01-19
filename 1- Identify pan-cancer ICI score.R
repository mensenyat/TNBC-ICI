setwd("D:/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer")

#### IMvigor210 - Melanoma ####
library("IMvigor210CoreBiologies")

data(cds)

meta = cds@phenoData@data
data = as.data.frame(cds@assayData$counts)
feat = cds@featureData@data

# Check that meta, feat and data are in the same order and convert Entrez ID to Gene Symbol #
table(feat$entrez_id == rownames(data))
table(colnames(data) == rownames(meta)) 

rownames(data) = make.names(feat$Symbol, unique = T)
write.table(data, "IMvigor210/IMvigor210 data.txt", sep = "\t")
write.table(meta, "IMvigor210/IMvigor210 metadata.txt", sep = "\t")

#### E-MTAB-3218 ####
library(affy)
setwd("D:/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer/E-MTAB-3218")

Data <- ReadAffy()
eset <- rma(Data)
data = as.data.frame(eset@assayData$exprs)

u219 = read.table("U219 Codes.txt", sep = "\t", header = T)
u219 = u219[order(match(u219$Array.ID, rownames(data))),]

table(rownames(data) == u219$Array.ID)

rownames(data) = make.names(u219$Gene.Symbol, unique = T)

meta = read.table("E-MTAB-3218.sdrf.txt", sep = "\t", header = T)

meta = meta[order(match(meta$Array.Data.File, colnames(data))),]
table(meta$Array.Data.File == colnames(data))

colnames(data) = meta$Source.Name

meta = meta[,-c(2,3,28,41:53)]

write.table(data, "E-MTAB-3218 data.txt", sep = "\t")
write.table(meta, "E-MTAB-3218 metadata.txt", sep = "\t")

#### E-MTAB-4030 ####
library(affy)
setwd("D:/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer/E-MTAB-4030")

Data <- ReadAffy()
eset <- rma(Data)
data = as.data.frame(eset@assayData$exprs)

u219 = read.table("U219 Codes.txt", sep = "\t", header = T)
u219 = u219[order(match(u219$Array.ID, rownames(data))),]

table(rownames(data) == u219$Array.ID)

rownames(data) = make.names(u219$Gene.Symbol, unique = T)

meta = read.table("E-MTAB-4030.sdrf.txt", sep = "\t", header = T)

meta = meta[meta$Array.Data.File %in% colnames(data),]
meta = meta[order(match(meta$Array.Data.File, colnames(data))),]
table(meta$Array.Data.File == colnames(data))

colnames(data) = meta$Source.Name

meta = meta[,-c(2,3,28,20:27)]

write.table(data, "E-MTAB-4030 data.txt", sep = "\t")
write.table(meta, "E-MTAB-4030 metadata.txt", sep = "\t")

#### GSE111636 ####
setwd("D:/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer/GSE111636")
data = read.table("GSE111636 data.txt", sep = "\t", header = T)

hta2 = read.table("GPL17586-45144.txt", sep = "\t", header = T)
hta2 = hta2[hta2$ID %in% data$ID_REF,]
hta2 = hta2[order(match(hta2$ID, data$ID_REF)),]
table(hta2$ID == data$ID_REF)
data = data[!is.na(hta2$GeneSymbol),]
hta2 = hta2[!is.na(hta2$GeneSymbol),]
table(hta2$ID == data$ID_REF)

data$ID_REF = make.names(hta2$GeneSymbol, unique = T)

write.table(data, file = "GSE111636 data.txt", sep = "\t", row.names = F)

#### GSE165252 ####
setwd("D:/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer/GSE165252")
data = read.table("GSE165252 data.txt", sep = "\t", header = T)

library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ens=getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"),  ### The type of data that you want to obtain
          filters = 'ensembl_gene_id',  ### The type of data you are submitting
          values = substr(data$X, 1, 15), ### The data you are submitting
          mart = ensembl)
ens = ens[which(ens$ensembl_gene_id %in% substr(data$X, 1, 15)),]
ens = ens[which(duplicated(ens$ensembl_gene_id) == "FALSE"),]
data = data[which(substr(data$X, 1, 15) %in% ens$ensembl_gene_id),]
data = data[which(duplicated(substr(data$X, 1, 15)) == "FALSE"),]

ens = ens[order(match(ens$ensembl_gene_id, substr(data$X, 1, 15))),]

data$X = make.names(ens$hgnc_symbol, unique = T)

write.table(data, file = "GSE165252 data.txt", sep = "\t", row.names = F)

#### GSE67501 ####
setwd("D:/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer/GSE67501")
data = read.table("GSE67501 data.csv", sep = ",", header = T)

ht12 = read.table("GPL14951-11332.txt", sep = "\t", header = T)
ht12 = ht12[ht12$ID %in% data$ID_REF,]
data = data[data$ID_REF %in% ht12$ID,]
ht12 = ht12[order(match(ht12$ID, data$ID_REF)),]
table(ht12$ID == data$ID_REF)
data = data[!is.na(ht12$ILMN_Gene),]
ht12 = ht12[!is.na(ht12$ILMN_Gene),]
table(ht12$ID == data$ID_REF)

data$ID_REF = make.names(ht12$ILMN_Gene, unique = T)

write.table(data, file = "GSE67501 data.txt", sep = "\t", row.names = F)


#### Create mega-list containing all normalized datasets ####
setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer")
dataGSE = list(NA)
metaGSE = list(NA)
namesGSE = c( "E-MTAB-4030", "GSE11636", "GSE165252", "GSE181815", "GSE183924", "GSE67501", "GSE78220", "IMvigor210",
              "E-MTAB-3218", "GSE100797")

#### E-MTAB-4030 ####
dataEMTAB4030 = read.table("E-MTAB-4030/E-MTAB-4030 data.txt", sep = "\t", header = T, row.names = 1)
metaEMTAB4030 = read.table("E-MTAB-4030/E-MTAB-4030 metadata.txt", sep = "\t", header = T)
metaEMTAB4030 = metaEMTAB4030[metaEMTAB4030$Characteristics.sampling.time.point. == "screen",]
metaEMTAB4030$Source.Name = paste("X", metaEMTAB4030$Source.Name, sep = "")
metaEMTAB4030 = metaEMTAB4030[!metaEMTAB4030$Characteristics.BestOverallResponse.RECIST. == "NonEvaluable",]

metaEMTAB4030$Response = NA
metaEMTAB4030$Response[metaEMTAB4030$Characteristics.BestOverallResponse.RECIST. %in% 
                         c("Complete/Partial Response", "Stable Disease")] = "ORR"
metaEMTAB4030$Response[metaEMTAB4030$Characteristics.BestOverallResponse.RECIST. == "Progressive Disease"] = "PD"

dataEMTAB4030 = dataEMTAB4030[,colnames(dataEMTAB4030) %in% metaEMTAB4030$Source.Name]
metaEMTAB4030 = metaEMTAB4030[order(match(metaEMTAB4030$Source.Name, colnames(dataEMTAB4030))),]

dataEMTAB4030t = as.data.frame(as.matrix(t(dataEMTAB4030)))
dataEMTAB4030t = as.data.frame(scale(dataEMTAB4030t))
dataEMTAB4030t[dataEMTAB4030t == "NaN"] = NA
dataEMTAB4030t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataEMTAB4030t))))))

rownames(metaEMTAB4030) = metaEMTAB4030$Source.Name

dataGSE[[1]] = dataEMTAB4030t
metaGSE[[1]] = metaEMTAB4030

#### GSE11636 ####
dataGSE11636 = read.table("GSE111636/GSE111636 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE11636 = read.table("GSE111636/GSE111636 metadata.csv", sep = ",", header = T)
rownames(metaGSE11636) = metaGSE11636$Sample_geo_accession

dataGSE11636t = as.data.frame(scale(as.matrix(t(dataGSE11636))))
dataGSE11636t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataGSE11636t))))))

dataGSE[[2]] = dataGSE11636t
metaGSE[[2]] = metaGSE11636

#### GSE165252 ####
dataGSE165252 = read.table("GSE165252/GSE165252 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE165252 = read.table("GSE165252/GSE165252 metadata.csv", sep = ",", header = T)
metaGSE165252 = metaGSE165252[metaGSE165252$Time == "baseline",]
metaGSE165252$Response[metaGSE165252$Response == "responder"] = "ORR"
metaGSE165252$Response[metaGSE165252$Response == "non-responder"] = "PD"

dataGSE165252 = dataGSE165252[,colnames(dataGSE165252) %in% metaGSE165252$Sample_description]
metaGSE165252 = metaGSE165252[order(match(metaGSE165252$Sample_description, colnames(dataGSE165252))),]

dataGSE165252t = as.data.frame(as.matrix(t(dataGSE165252)))
dataGSE165252t = as.data.frame(scale(dataGSE165252t))
dataGSE165252t[dataGSE165252t == "NaN"] = NA
dataGSE165252t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataGSE165252t))))))

rownames(metaGSE165252) = metaGSE165252$Sample_description

dataGSE[[3]] = dataGSE165252t
metaGSE[[3]] = metaGSE165252

#### GSE181815 ####
dataGSE181815 = read.table("GSE181815/GSE181815 data.csv", sep = ",", header = T)
rownames(dataGSE181815) = make.names(dataGSE181815$GeneSymbol, unique = T)
dataGSE181815 = dataGSE181815[,-1]

colnames(dataGSE181815) = substr(colnames(dataGSE181815), 1, 7)
metaGSE181815 = read.table("GSE181815/GSE181815 metadata.txt", sep = "\t", header = T)
metaGSE181815$characteristics..Response.to.pembrolizumab[metaGSE181815$characteristics..Response.to.pembrolizumab %in% 
                                                           c("Progressive Disease")] = "PD"
metaGSE181815$characteristics..Response.to.pembrolizumab[metaGSE181815$characteristics..Response.to.pembrolizumab %in% 
                                                           c("Partial Response", "Complete Response", "Stable Disease")] = "ORR"
metaGSE181815$Sample.name = gsub("-", ".", metaGSE181815$Sample.name)
metaGSE181815 = metaGSE181815[metaGSE181815$Sample.name %in% colnames(dataGSE181815),]
metaGSE181815 = metaGSE181815[order(match(metaGSE181815$Sample.name, colnames(dataGSE181815))),]

dataGSE181815 = dataGSE181815[,colnames(dataGSE181815) %in% metaGSE181815$Sample.name]
dataGSE181815t = as.data.frame(scale(as.matrix(t(dataGSE181815))))
dataGSE181815t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataGSE181815t))))))

rownames(metaGSE181815) = metaGSE181815$Sample.name
metaGSE181815$Response = metaGSE181815$characteristics..Response.to.pembrolizumab

dataGSE[[4]] = dataGSE181815t
metaGSE[[4]] = metaGSE181815

#### GSE183924 ####
dataGSE183924 = read.table("GSE183924/GSE183924 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE183924 = read.table("GSE183924/GSE183924 metadata.txt", sep = "\t", header = T)

metaGSE183924$Response = NA
metaGSE183924$Response[metaGSE183924$Relapse == "Yes"] = "PD"
metaGSE183924$Response[metaGSE183924$Relapse == "No"] = "ORR"
metaGSE183924 = metaGSE183924[!is.na(metaGSE183924$Response),]

dataGSE183924 = dataGSE183924[,colnames(dataGSE183924) %in% metaGSE183924$RNA.Seq.ID]
metaGSE183924 = metaGSE183924[order(match(metaGSE183924$RNA.Seq.ID, colnames(dataGSE183924))),]

dataGSE183924t = as.data.frame(as.matrix(t(dataGSE183924)))
dataGSE183924t = as.data.frame(scale(dataGSE183924t))
dataGSE183924t[dataGSE183924t == "NaN"] = NA
dataGSE183924t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataGSE183924t))))))

rownames(metaGSE183924) = metaGSE183924$RNA.Seq.ID

dataGSE[[5]] = dataGSE183924t
metaGSE[[5]] = metaGSE183924

#### GSE67501 ####
dataGSE67501 = read.table("GSE67501/GSE67501 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE67501 = read.table("GSE67501/GSE67501 metadata.csv", sep = ",", header = T)

dataGSE67501 = dataGSE67501[,colnames(dataGSE67501) %in% metaGSE67501$Sample_geo_accession]
metaGSE67501 = metaGSE67501[order(match(metaGSE67501$Sample_geo_accession, colnames(dataGSE67501))),]

dataGSE67501t = as.data.frame(as.matrix(t(dataGSE67501)))
dataGSE67501t = as.data.frame(scale(dataGSE67501t))
dataGSE67501t[dataGSE67501t == "NaN"] = NA
dataGSE67501t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataGSE67501t))))))

rownames(metaGSE67501) = metaGSE67501$Sample_geo_accession

dataGSE[[6]] = dataGSE67501t
metaGSE[[6]] = metaGSE67501

#### GSE78220 ####
dataGSE78220 = read.table("GSE78220/GSE78220 data.txt", sep = "\t", header = T, row.names = 1)
metaGSE78220 = read.table("GSE78220/GSE78220 metadata.csv", sep = ",", header = T)
metaGSE78220$Response[metaGSE78220$Response == "Progressive Disease"] = "PD"
metaGSE78220$Response[metaGSE78220$Response %in% c("Partial Response", "Complete Response", "Stable Disease")] = "ORR"

dataGSE78220 = dataGSE78220[,colnames(dataGSE78220) %in% metaGSE78220$Patient]
metaGSE78220 = metaGSE78220[order(match(metaGSE78220$Patient, colnames(dataGSE78220))),]

dataGSE78220t = as.data.frame(as.matrix(t(dataGSE78220)))
dataGSE78220t = as.data.frame(scale(dataGSE78220t))
dataGSE78220t[dataGSE78220t == "NaN"] = NA
dataGSE78220t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataGSE78220t))))))

rownames(metaGSE78220) = metaGSE78220$Patient

dataGSE[[7]] = dataGSE78220t
metaGSE[[7]] = metaGSE78220

#### IMvigor210 ####
dataIMvigor210 = read.table("IMvigor210/IMvigor210 data.txt", sep = "\t", header = T, row.names = 1)
metaIMvigor210 = read.table("IMvigor210/IMvigor210 metadata.txt", sep = "\t", header = T)
metaIMvigor210 = metaIMvigor210[!metaIMvigor210$Best.Confirmed.Overall.Response == "NE",]
metaIMvigor210$Best.Confirmed.Overall.Response[metaIMvigor210$Best.Confirmed.Overall.Response %in% c("PR", "CR", "SD")] = "ORR"
metaIMvigor210$Response = metaIMvigor210$Best.Confirmed.Overall.Response

dataIMvigor210 = dataIMvigor210[,colnames(dataIMvigor210) %in% rownames(metaIMvigor210)]
metaIMvigor210 = metaIMvigor210[order(match(rownames(metaIMvigor210), colnames(dataIMvigor210))),]

dataIMvigor210t = as.data.frame(as.matrix(t(dataIMvigor210)))
dataIMvigor210t = as.data.frame(scale(dataIMvigor210t))
dataIMvigor210t[dataIMvigor210t == "NaN"] = NA
dataIMvigor210t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataIMvigor210t))))))

dataGSE[[8]] = dataIMvigor210t
metaGSE[[8]] = metaIMvigor210

#### E-MTAB-3218 ####
dataEMTAB3218 = read.table("E-MTAB-3218/E-MTAB-3218 data.txt", sep = "\t", header = T, row.names = 1)
metaEMTAB3218 = read.table("E-MTAB-3218/E-MTAB-3218 metadata.txt", sep = "\t", header = T)
metaEMTAB3218 = metaEMTAB3218[metaEMTAB3218$Characteristics.biopsy.timepoint. == "Screen",]
metaEMTAB3218$Source.Name = paste("X", metaEMTAB3218$Source.Name, sep = "")
metaEMTAB3218 = metaEMTAB3218[!metaEMTAB3218$Characteristics.BOR. == "NE",]

metaEMTAB3218$Response = NA
metaEMTAB3218$Response[metaEMTAB3218$Characteristics.BOR3. %in% c("CRPR", "SD")] = "ORR"
metaEMTAB3218$Response[metaEMTAB3218$Characteristics.BOR3. == "PD"] = "PD"

dataEMTAB3218 = dataEMTAB3218[,colnames(dataEMTAB3218) %in% metaEMTAB3218$Source.Name]
metaEMTAB3218 = metaEMTAB3218[order(match(metaEMTAB3218$Source.Name, colnames(dataEMTAB3218))),]

dataEMTAB3218t = as.data.frame(as.matrix(t(dataEMTAB3218)))
dataEMTAB3218t = as.data.frame(scale(dataEMTAB3218t))
dataEMTAB3218t[dataEMTAB3218t == "NaN"] = NA
dataEMTAB3218t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataEMTAB3218t))))))

rownames(metaEMTAB3218) = metaEMTAB3218$Source.Name

dataGSE[[9]] = dataEMTAB3218t
metaGSE[[9]] = metaEMTAB3218


#### GSE100797 ####
dataGSE100797 = read.table("GSE100797/GSE100797_Data.txt", sep = "\t", header = T, row.names = 1)
metaGSE100797 = read.table("GSE100797/GSE100797_Metadata.txt", sep = "\t", header = T)

dataGSE100797 = dataGSE100797[,colnames(dataGSE100797) %in% metaGSE100797$Patient]
metaGSE100797 = metaGSE100797[order(match(metaGSE100797$Patient, colnames(dataGSE100797))),]
metaGSE100797$Response[metaGSE100797$Response %in% c("CR", "PR", "SD")] = "ORR"

dataGSE100797t = as.data.frame(as.matrix(t(dataGSE100797)))
dataGSE100797t = as.data.frame(scale(dataGSE100797t))
dataGSE100797t[dataGSE100797t == "NaN"] = NA
dataGSE100797t = as.data.frame(as.matrix(t(na.omit(as.matrix(t(dataGSE100797t))))))

rownames(metaGSE100797) = metaGSE100797$Patient

dataGSE[[10]] = dataGSE100797t
metaGSE[[10]] = metaGSE100797

names(dataGSE) = namesGSE
names(metaGSE) = namesGSE

#### Identify DEG and give score ####
setwd("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer")
load("C:/Users/miki7/OneDrive/R/BRCA/SMARCA4/TONIC/ORR Signature/Pan-cancer/Environment All Non-TNBC ICI Databases.RData")

genes = NA

for(i in 1:10){
  colnames(dataGSE[[i]]) = gsub("-", ".", colnames(dataGSE[[i]]))
  shortGene = substr(colnames(dataGSE[[i]]), nchar(colnames(dataGSE[[i]]))-1, nchar(colnames(dataGSE[[i]])))
  dataGSE[[i]] = dataGSE[[i]][,which(!shortGene %in% c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", ".10", ".11"))] # Remove genes with >1 version
  genes = c(genes, colnames(dataGSE[[i]]))
  print(i)
}

genes = unique(genes)
genes = na.omit(genes)

geneScore = as.data.frame(matrix(NA, nrow = length(genes), ncol = 10))
rownames(geneScore) = genes
colnames(geneScore) = namesGSE

for(l in 1:10){
  data = dataGSE[[l]]
  meta = metaGSE[[l]]
  
  metaORR = meta[meta$Response == "ORR",]
  metaPD = meta[meta$Response == "PD",]
  
  dataORR = data[rownames(data) %in% rownames(metaORR),]
  dataPD = data[rownames(data) %in% rownames(metaPD),]
  
  meanORR = apply(dataORR, 2, mean)
  meanPD = apply(dataPD, 2, mean)
  sdORR = apply(dataORR, 2, sd)
  sdORR[is.na(sdORR)] = 0
  sdPD = apply(dataPD, 2, sd)
  sdPD[is.na(sdPD)] = 0
  
  data = data[,sdORR>0 & sdPD>0]
  dataORR = data[rownames(data) %in% rownames(metaORR),]
  dataPD = data[rownames(data) %in% rownames(metaPD),]
  meanORR = apply(dataORR, 2, mean)
  meanPD = apply(dataPD, 2, mean)
  
  fold = (meanORR-meanPD)
  Zratio = (fold/sd(fold))
  
  # Compute statistical significance #
  pvalue = NULL
  for(i in 1 : ncol(dataORR)) {
    pvalue[i] = t.test(as.numeric(dataORR[,i]), as.numeric(dataPD[,i]))$p.value
  }
  
  score = NA
  score = -log10(pvalue) * log2(nrow(data)) 
  score[score > 2*log2(nrow(data))] = 2*log2(nrow(data)) # Equivalent to p-value = 0.01
  score2 = score
  
  for(i in 1:length(fold)){
    if(fold[i] < 0){
      score[i] = -score[i]
    }
  }
  names(score) = colnames(data)
  
  for(i in 1:length(score)){
    geneScore[which(rownames(geneScore) == names(score)[i]),l] = score[i]
  }
  print(l)
}

finalScore = rowSums(geneScore, na.rm = T)
View(finalScore[!abs(finalScore) < 40])

final2 = as.data.frame(cbind(names(finalScore), finalScore))
final2$finalScore = as.numeric(final2$finalScore)

write.table(final2, file = "Score immunotherapy pan-cancer.txt", sep = "\t")
