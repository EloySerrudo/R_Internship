install.packages("readxl")

library("readxl")

# Dukhovny (2019) Los CSV no estan incluidos en la tabla ¿xq?
my_data <- read.table(file = '1-s2-0.txt', sep = ',', header = FALSE)
colnames(my_data) <- c("Gene ID","MAGeCK1",	"MAGeCK2",	"RIGER 1",	"RIGER 2")

dataSavidis <- read_excel("ZKVData/1-s2.0-S2211124716307689-mmc8.xlsx", sheet = "Summary CME Compare") # Savidis
dataLi <- read_excel("ZKVData/pnas.1900867116.sd01.xlsx", sheet = "zika_gw_5th_best") # Li
dataWang1 <- read_excel("ZKVData/NIHMS1553325-supplement-2.xlsx", sheet = "Ranking", col_names = FALSE, skip = 1) # Wang GSC
dataWang2 <- read_excel("ZKVData/NIHMS1553325-supplement-3.xlsx", sheet = "Sheet1", col_names = FALSE, skip = 2) # Wang 293FT
dataRother <- read_excel("ZKVData/1-s2.0-S0168170221000459-mmc1.xlsx", sheet = "(B) Gene ranking") # Rother
dataShue <- read_excel("ZKVData/jvi.00596-21-s0001.xls", sheet = "Genelist") # Shue

colnames(dataLi) <- as.vector(dataLi[1,])
colnames(dataWang1) <- as.vector(dataWang1[1,])
colnames(dataWang2) <- as.vector(dataWang2[1,])
dataShue <- dataShue[order(dataShue$deseq2.FC, decreasing = TRUE),]

Hits.Savidis <- dataSavidis[-1,-2] # Ya está ordenado
Hits.Li <- dataLi[-1,-(2:4)] # Está ordenado por defecto con zika.init.5th
Hits.Wang1 <- dataWang1[-1,-(2:3)] # Ya está ordenado
Hits.Wang2 <- dataWang2[-1,-(2:5)] # Ya está ordenado. El primer Gen MMGT1 (EMC5)
Hits.Rother <- dataRother[,-(2:6)] # Utilizar la hoja Gene Ranking.Ya está ordenado por la columna Rank
Hits.Shue <- dataShue[,-(2:7)] # La columna deseq2.FC ordena los 500

rm(dataSavidis)
rm(dataLi)
rm(dataWang1)
rm(dataWang2)
rm(dataRother)
rm(dataShue)

Hits.Wang2[[1,1]] <- "MMGT1"
Hits.Li <- Hits.Li[1:500,]
Hits.Shue <- Hits.Shue[1:500,]

sum(unlist(Hits_Wang) %in% unlist(Hits_Savidis))
sum(unlist(Hits_Wang) %in% unlist(Hits_Rother))
sum(unlist(Hits_Wang) %in% unlist(Hits_Dukhovny))

sum(unlist(Hits.Li) %in% unlist(Hits.Shue))

genes <- unlist(Hits.Li)

namesRother <- unlist(Hits.Rother[!(unlist(Hits.Rother) %in% genes),])


genes2 <- c(genes, namesRother)
