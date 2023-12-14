# Install -----------------------------------------------------------------


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db", version = "3.16")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.16")


# Machine Learning --------------------------------------------------------


library(data.table)
# library(gprofiler2)


dataCombined <- fread("ZIKVData/combinedFinal.csv", sep=',', 
                      header = TRUE, fill = TRUE)

dataCombined[17572:17583, 1:10]

typeof(dataCombined)
class(dataCombined)

# genes.HDF <- genes.HDF.DF |> dplyr::select(genes, ensembl.ids)

genes.input <- dataCombined |> filter(Ensembl_ID %in% genes.HDF.DF$ensembl.ids)
z.score <- as.data.frame(scale(genes.input[,-1])) |> 
  mutate(Ensembl_ID = genes.input$Ensembl_ID, .before = 1)

std.dev <- apply(z.score[,-1], MARGIN = 2, FUN = sd)
mean.dt <- apply(z.score[,-1], MARGIN = 2, FUN = mean)
coef.var <- std.dev / mean.dt
coef.var[1:20]

coef.var.rank <- base::rank(1/coef.var)
plot(coef.var.rank, coef.var)

coef.var.rank[1:10]

write_csv(genes.input, "ZIKVData/genesInput.csv")
