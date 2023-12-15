# Install -----------------------------------------------------------------


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db", version = "3.16")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.16")



# Data combined -----------------------------------------------------------

# library(gprofiler2)


dataCombined <- fread("ZIKVData/combinedFinal.csv", sep=',', 
                      header = TRUE, fill = TRUE)

dataCombined[17572:17583, 1:10]

typeof(dataCombined)
class(dataCombined)

# genes.HDF <- genes.HDF.DF |> dplyr::select(genes, ensembl.ids)

# Machine Learning --------------------------------------------------------

library(data.table)

genes.input <- fread("ZIKVData/genesInput.csv", sep = ',', 
                     header = TRUE, fill = TRUE)
names <- genes.input$Ensembl_ID

genes.input <- genes.input |> dplyr::select(-1)

genes.input <- as.data.frame(genes.input, row.names = names) |> 
  dplyr::select(-where(~all(. == 0)))

z.score <- as.data.frame(scale(genes.input), row.names = names)

coef.var <- data.frame(
  sd = apply(z.score, MARGIN = 2, FUN = sd),
  mean = apply(z.score, MARGIN = 2, FUN = mean)) |> 
  mutate(cv = sd / mean)

hist(coef.var$cv, breaks = 1)

sum(coef.var$cv < 0.01)

variance.rank <- rank(1/coef.var$sd)
plot(variance.rank, coef.var$cv, xlab = "Genes", ylab = "Variance")

sum(coef.var$sd > 0.01)
