install.packages("readxl")

library("readxl")

path <- "1-s2.0-S2211124716307689-mmc8.xlsx"
my_data <- read_excel(path, sheet = "Summary CME Compare")
my_data <- read_excel("NIHMS1553325-supplement-2.xlsx", sheet = "Ranking")
my_data <- read_excel("1-s2.0-S0168170221000459-mmc1.xlsx", sheet = "(B) Gene ranking")
my_data <- read.table(file = '1-s2-0.txt', sep = ',', header = FALSE)
colnames(my_data) <- c("Gene ID","MAGeCK1",	"MAGeCK2",	"RIGER 1",	"RIGER 2")

Hits_Savidis <- my_data[-1,-2]
Hits_Wang <- my_data[-1,-c(2,3)]
Hits_Rother <- my_data[,-(2:6)]
Hits_Dukhovny <- my_data[1]

sum(unlist(Hits_Wang) %in% unlist(Hits_Savidis))
sum(unlist(Hits_Wang) %in% unlist(Hits_Rother))
sum(unlist(Hits_Wang) %in% unlist(Hits_Dukhovny))
