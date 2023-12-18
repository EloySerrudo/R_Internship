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

# Pre-processing ----------------------------------------------------------

library(data.table)

genes.input <- fread("ZIKVData/genesInput.csv", sep = ',', 
                     header = TRUE) |> 
  as.data.frame()

names <- genes.input$Ensembl_ID

genes.input <- genes.input |> dplyr::select(-1)

genes.input <- as.data.frame(genes.input, row.names = names) |> 
  dplyr::select(-where(~all(. == 0)))

z.score <- as.data.frame(scale(genes.input), row.names = names)

coef.var <- data.frame(
  sd = apply(z.score, MARGIN = 2, FUN = sd),
  mean = apply(z.score, MARGIN = 2, FUN = mean)) |> 
  mutate(cv = abs(sd / mean))

# hist(coef.var$cv, breaks = 1, main = "CV", xlab = "CV")

features.rank <- rank(1/coef.var$cv)
plot(features.rank, 
     coef.var$cv, 
     xlab = "Features", ylab = "Coefficient of variation ",
     xlim = c(0,800))

sum(coef.var$cv > 5e16)

features.sel <- which(coef.var$cv>quantile(coef.var$cv,0.99)) # cv > 0.01

genes.fil.1 <- z.score |> dplyr::select(all_of(features.sel))

corr.matrix <- cor(genes.fil.1)

correlation.sel <- which(abs(corr.matrix) > 0.7, arr.ind = TRUE)

features.pairs.sel <- cbind(rownames(corr.matrix)[correlation.sel[, 1]],
                            colnames(corr.matrix)[correlation.sel[, 2]])


### Set up data structure and filter highly correlated columns
fb <- corr.matrix
result <- vector(mode = "list", length = 605)
for ( i in 1:605) {
  result[[i]] <- {
    if ( i == nrow(fb)) {
      which(abs(unlist(fb[i, 1:i])) > 0.7)
    } else {
      which(abs(c(unlist(fb[i ,1:i]), unlist(fb[(i+1):nrow(fb), i]))) > 0.7)
    }
  }
}

# correlation.sel es lo mismo que result

### Iteratively remove most correlated feature
maximum = 1000
corr.list <- result
while(maximum > 1) {
  lengths <- as.numeric(map(result, function(y) length(y))) # Calcula la cantidad de otras columnas con un pearson mayor a 0.7 con esa columna
  maximum <- max(lengths) # Encuentra la columna con cantidad de correlaciones >0.7
  longest <- which(lengths == maximum) # Toma su índice
  #print(maximum)
  if (maximum == 1) { break }
  ### in case there is more than one maximum correlated feature choose random
  j <- sample(1:length(longest), 1)
  index <- longest[[j]]
  ### remove feature from all other rows
  for (k in result[[index]]) {
    result[[k]] <- result[[k]][result[[k]] != index]
  }
  ### remove row of most correlated feature
  result[[index]] <- c(0) # Elimina la columna con más correlaciones de la tabla
}
### Calculate columns to be removed
result <- which(unlist(result) == 0) + 1 # index shift to fit original data

rm(fb, i, index, j, k, lengths, longest, maximum)

genes.fil.2 <- genes.fil.1 |> dplyr::select(all_of(result))

rm(coef.var, corr.matrix, correlation.sel, features.pairs.sel, genes.fil.1,
   genes.input, z.score, variance.rank, variance.sel)

result <- result - 1

# Machine Learning --------------------------------------------------------

gene.data <- fread("ZIKVData/combinedFinal.csv") %>% as.data.frame
gene.data <- gene.data |> dplyr::select(all_of(result))

GS1.of.7 <- GS1.of.7[GS1.of.7 %in% names]
GS2.of.7 <- GS2.of.7[GS2.of.7 %in% names]
GS3.of.7 <- GS3.of.7[GS3.of.7 %in% names]
