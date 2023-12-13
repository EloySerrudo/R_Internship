if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db", version = "3.16")

library(tidyverse)
library(data.table)
library(gprofiler2)

library(org.Hs.eg.db)

dataCombined <- fread("ZIKVData/dataReducedReduced.csv", sep=',', 
                      header = TRUE, fill = TRUE)

dataCombined[17572:17583, 1:10]

typeof(dataCombined)
class(dataCombined)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.16")


library(biomaRt)

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_genes <- getBM(attributes=c('ensembl_gene_id', 
                                    'hgnc_symbol', 
                                    'external_gene_name', 
                                    'external_synonym'), mart = ensembl)

head(chr1_genes)

listAttributes(ensembl)[1:38,1:2]

# Not important -----------------------------------------------------------



simbolo_hgnc <- mapIds(org.Hs.eg.db, 
                       keys = descripcion_del_gen, 
                       column = "SYMBOL", 
                       keytype = "GENENAME")

keytypes(org.Hs.eg.db)

ensembl_ids <- mapIds(org.Hs.eg.db, 
                      keys = genes.HDF, 
                      column = "ENSEMBL", 
                      keytype = "SYMBOL")

ens.Savidis <- mapIds(org.Hs.eg.db, 
                      keys = Hits.Savidis, 
                      column = "ENSEMBL", 
                      keytype = "ALIAS")


unique(ensembl_ids)

gconvert(
  query = genes.HDF,
  organism = "hsapiens",
  target = "ENSG",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE
)
