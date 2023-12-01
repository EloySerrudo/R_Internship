
# Installing DESeq2 -------------------------------------------------------

sessionInfo()
.libPaths()

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.16")
BiocManager::install("DESeq2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")


library(tidyverse)
library(DESeq2)
library(readxl)

# Refreshing how to use DESeq2 --------------------------------------------

rawCounts <- read.table("inData/SARSCoV2_rawCounts.txt")
sampleTable <- read.table("inData/SARSCoV2_sampleInfo.txt")

rawCounts2 <- read.table("inData/sarsCov2_rawCounts_smybols.txt")
colNamesList <- strsplit(colnames(rawCounts2), 'X')
colnames(rawCounts2) <- sapply(colNamesList, '[[', 2)

sampleTable2 <- read.table("inData/sarsCov2_sampleTable_smybols.txt")
str(sampleTable)
str(sampleTable2)
sampleTable2 <- sampleTable2 |> 
  mutate(timepoint = factor(timepoint, c('0h','3h','6h','12h','24h')))

ddsMat1 <- DESeqDataSetFromMatrix(countData = rawCounts, 
                                 colData = sampleTable, 
                                 design = ~virus)

ddsMat2 <- DESeqDataSetFromMatrix(countData = rawCounts2, 
                                  colData = sampleTable2, 
                                  design = ~timepoint)

dds1 <- DESeq(ddsMat1)
dds2 <- DESeq(ddsMat2)

res1 <- results(dds1)
res2 <- results(dds2)

summary(results(dds1, alpha = 0.05, lfcThreshold = 0.5))
summary(results(dds2, alpha = 0.05, lfcThreshold = 0.5))

normCounts <- estimateSizeFactors(ddsMat1)
normCounts <- counts(normCounts, normalized= TRUE)

ddsMat1

# DESeq2 ZIKA -------------------------------------------------------------

# Infected: JEG-3, U-251 MG, and HK-2 cell lines
# - HK-2: Epithelial, Human kidney, Cortex, Proximal tubule, Papilloma
# - U-251 MG: Human Brain, Gliobastoma
# - JEG-3: Epithelial, Human Placenta, Choriocarcinoma
# uninfected controls. 

# The RNA-seq was performed in infected cell lines under mock and ZIKV-infected
# conditions, including placental cell line (JEG-3), nerve cell line (U-251 MG),
# and kidney cell line (HK-2).
# The placental cell line was measured at 3 h,12 h, and 24 h post-ZIKV infection,
# while the other two cell lines were uniformly measured at 
# 24 h post infection (h.p.i.) merely. 

# Viral strain: ZIKV African strain MR766 (MOI = 1)

ddsMat1 <- DESeqDataSetFromMatrix(countData = rawCounts,
                                 colData = sampleTable,
                                 design = ~virus)
dds <- DESeq(ddsMat1)

# Creating the Data -----------------------------------------------------

U251MG.Table <- data.frame(sample = c("U1", "U2", "U3", "UZV_1", "UZV_2", "UZV_3"),
                           sampleNr = c(1, 2, 3, 4, 5, 6),
                           timepoint = factor(sampleTable1$timepoint),
                           virus = sampleTable1$virus,
                           row.names = c("U1", "U2", "U3", "UZV_1", "UZV_2", "UZV_3"))

HK2.Table <- data.frame(sample = c("HRP1", "HRP2", "HRP3", "HRZV1", "HRZV2", "HRZV3"),
                        sampleNr = c(1, 2, 3, 4, 5, 6),
                        timepoint = factor(sampleTable1$timepoint),
                        virus = sampleTable1$virus,
                        row.names = c("HRP1", "HRP2", "HRP3", "HRZV1", "HRZV2", "HRZV3"))

U251MG.FPKM <- read_excel("ZIKVData/12864_2022_8919_MOESM1_ESM.xlsx", 
                          sheet = "U-251 MG cell gene expression")
rownames(U251MG.FPKM) <- U251MG.FPKM$GeneName
length(unique(U251MG.FPKM$GeneName))
a <- duplicated(U251MG.FPKM$GeneName)
U251MG.FPKM[a,]

HK2.FPKM <- read_excel("ZIKVData/12864_2022_8919_MOESM1_ESM.xlsx", 
                       sheet = "HK-2 cells gene expression")
sum(duplicated(HK2.FPKM$GeneName))

ddsMat <- DESeqDataSetFromMatrix(countData = U251MG.FPKM[2:7],
                                 colData = U251MG.Table,
                                 design = ~virus)

aux.table <- round(U251MG.FPKM[2:7])
