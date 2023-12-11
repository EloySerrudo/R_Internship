library(tidyverse)
library(data.table)


dataCombined <- fread("ZIKVData/combinedFinal.csv", sep=',', 
                      header = TRUE, fill = TRUE)

dataCombined[17572:17583, 1:10]

typeof(dataCombined)
class(dataCombined)

