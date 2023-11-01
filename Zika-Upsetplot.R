library(tidyverse)
library(readxl)
library(ComplexUpset)

# Dukhovny (2019) Los CSV no estan incluidos en la tabla ¿xq?
#my_data <- read.table(file = '1-s2-0.txt', sep = ',', header = FALSE)
#colnames(my_data) <- c("Gene ID","MAGeCK1",	"MAGeCK2",	"RIGER 1",	"RIGER 2")

# Loading Datasets --------------------------------------------------------

dataSavidis <- read_excel("1-s2.0-S2211124716307689-mmc8.xlsx", sheet = "Summary CME Compare")[-1,] # Savidis

dataLi <- read_excel("pnas.1900867116.sd01.xlsx", sheet = "zika_gw_5th_best") # Li
colnames(dataLi) <- as.vector(dataLi[1,])
dataLi <- dataLi[-1,]
cols.num <- c("sgRNAs", "zika.init.5th", "p.bh")
dataLi[cols.num] <- sapply(dataLi[cols.num],as.numeric)

dataWangGSC <- read_excel("NIHMS1553325-supplement-2.xlsx", sheet = "Ranking", col_names = FALSE, skip = 1)
dataWang293FT <- read_excel("NIHMS1553325-supplement-3.xlsx", sheet = "Sheet1", col_names = FALSE, skip = 2)
colnames(dataWangGSC) <- as.vector(dataWangGSC[1,])
colnames(dataWang293FT) <- as.vector(dataWang293FT[1,])
dataWangGSC <- dataWangGSC[-1,]
dataWang293FT <- dataWang293FT[-1,]

dataRother <- read_excel("1-s2.0-S0168170221000459-mmc1.xlsx", sheet = "(B) Gene ranking")

dataShue <- read_excel("jvi.00596-21-s0001.xls", sheet = "Genelist")
cols.num <- c("deseq2.FC", "deseq2.pval", "mageck.rank.pos", "mageck.fdr.pos", 
              "mageck.rank.neg", "mageck.fdr.neg")
dataShue[cols.num] <- sapply(dataShue[cols.num],as.numeric)
dataShue <- dataShue[order(dataShue$deseq2.pval),]


# Extracting Hits ---------------------------------------------------------


Hits.Savidis <- dataSavidis[,-2]
Hits.Li <- dataLi[,-(2:4)] # Está ordenado por defecto con zika.init.5th
Hits.WangGSC <- dataWang1[,-(2:3)] # Ya está ordenado
Hits.Wang293FT <- dataWang2[,-(2:5)] # Ya está ordenado. El primer Gen MMGT1 (EMC5)
Hits.Rother <- dataRother[,-(2:6)] # Utilizar la hoja Gene Ranking.Ya está ordenado por la columna Rank
Hits.Shue <- dataShue[,-(2:7)] # La columna deseq2.FC ordena los 500

rm(dataSavidis, dataLi, dataWang1, dataWang2, dataRother, dataShue)

Hits.Wang293FT[[1,1]] <- "MMGT1"
Hits.Li <- Hits.Li[1:500,]
Hits.Shue <- Hits.Shue[1:500,]

genes <- unlist(Hits.Li)
genes <- c(genes, unlist(Hits.Rother[!(unlist(Hits.Rother) %in% genes),]))
genes <- c(genes, unlist(Hits.Savidis[!(unlist(Hits.Savidis) %in% genes),]))
genes <- c(genes, unlist(Hits.Shue[!(unlist(Hits.Shue) %in% genes),]))
genes <- c(genes, unlist(Hits.Wang1[!(unlist(Hits.Wang1) %in% genes),]))
genes <- c(genes, unlist(Hits.Wang2[!(unlist(Hits.Wang2) %in% genes),]))
genes <- sort(genes)

Li <- genes %in% unlist(Hits.Li)
Rother <- genes %in% unlist(Hits.Rother)
Savidis <- genes %in% unlist(Hits.Savidis)
Shue <- genes %in% unlist(Hits.Shue)
Wang.GSC <- genes %in% unlist(Hits.Wang1)
Wang.293FT <- genes %in% unlist(Hits.Wang2)

genes.DF <- data.frame(genes, Savidis, Li, Wang.GSC, Wang.293FT, Rother, Shue,
                       row.names = NULL)
papers <- colnames(genes.DF)[-1]

upset(
  data = genes.DF,
  intersect = papers,
  name = 'Papers', 
  width_ratio=0.1
)

rm(Savidis, Li, Wang.GSC, Wang.293FT, Rother, Shue)

dataLi2 <- read_excel("pnas.1900867116.sd01.xlsx", sheet = "zika_gw_5th_best") # Li
colnames(dataLi2) <- as.vector(dataLi2[1,])
dataLi2 <- dataLi2[-1,]
dataLi <- dataLi[-1,]

dataLi2 <- dataLi2[order(dataLi2$p.bh),]
