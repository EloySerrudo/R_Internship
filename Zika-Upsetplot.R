library(tidyverse)
library(readxl)
#library(writexl)
#library(ComplexUpset)
#library(biomaRt)
#library(org.Hs.eg.db)
#library(gprofiler2)



# Loading Datasets --------------------------------------------------------

dataSavidis <- read_excel("ZIKVData/1-s2.0-S2211124716307689-mmc8.xlsx", 
                          sheet = "Summary CME Compare")[-1,]

dataLi <- read_excel("ZIKVData/pnas.1900867116.sd01.xlsx", 
                     sheet = "zika_gw_5th_best")
colnames(dataLi) <- as.vector(dataLi[1,])
dataLi <- dataLi[-1,]
cols.num <- c("sgRNAs", "zika.init.5th", "p.bh")
dataLi[cols.num] <- sapply(dataLi[cols.num],as.numeric)
dataLi <- dataLi[order(dataLi$p.bh),]

dataDukhovny <- read.table("ZIKVData/BIOGRID-ORCS-SCREEN_1209-1.1.14.screen.tab.txt",
                           header = TRUE, sep = "\t", quote = "", strip.white = TRUE,
                           fill = TRUE, comment.char = "") |> 
  filter(!duplicated(OFFICIAL_SYMBOL)) |> 
  filter(IDENTIFIER_TYPE != "UNKNOWN") |> 
  arrange(desc(row_number()))
#indexDukhovny <- read.table("ZIKVData/BIOGRID-ORCS-SCREEN_INDEX-1.1.14.index.tab.txt",
#                            header = TRUE, sep = "\t", quote = "", strip.white = TRUE,
#                            fill = TRUE, comment.char = "")

# filter(!grepl("^[0-9]+$", OFFICIAL_SYMBOL)) |> 
# filter(grepl("^[A-Z0-9]+$", OFFICIAL_SYMBOL)) |> 

# mutate(Gene = ifelse(grepl("^[A-Z0-9]+$", id), id, str_extract(id, "\\(([^)]+)\\)$")), 
#        .after = id) |> 
# mutate(Gene = gsub("\\(|\\)", "", Gene)) |> 
# arrange(neg.rank) |> 
# filter(!is.na(Gene))

dataWangGSC <- read_excel("ZIKVData/NIHMS1553325-supplement-2.xlsx", 
                          sheet = "Ranking", col_names = FALSE, skip = 1)
dataWang293FT <- read_excel("ZIKVData/NIHMS1553325-supplement-3.xlsx", 
                            sheet = "Sheet1", col_names = FALSE, skip = 2)
colnames(dataWangGSC) <- as.vector(dataWangGSC[1,])
colnames(dataWang293FT) <- as.vector(dataWang293FT[1,])
dataWangGSC <- dataWangGSC[-1,]
dataWang293FT <- dataWang293FT[-1,]
dataWang293FT[[1,1]] <- "MMGT1"

dataRother <- read_excel("ZIKVData/1-s2.0-S0168170221000459-mmc1.xlsx", 
                         sheet = "(B) Gene ranking")

dataShue <- read_excel("ZIKVData/jvi.00596-21-s0001.xls", sheet = "Genelist")
cols.num <- c("deseq2.FC", "deseq2.pval", "mageck.rank.pos", "mageck.fdr.pos", 
              "mageck.rank.neg", "mageck.fdr.neg")
dataShue[cols.num] <- sapply(dataShue[cols.num], as.numeric)
dataShue <- dataShue |> arrange(deseq2.pval, desc(deseq2.FC))
rm(cols.num)


# Extracting Hits ---------------------------------------------------------

Hits.Savidis <- sort(dataSavidis$`Zika CRISPR Screen Hits - Top 100`)
Hits.Li <- sort(dataLi$Gene[1:500])
Hits.Dukhovny <- sort(dataDukhovny$OFFICIAL_SYMBOL[1:500])
Hits.WangGSC <- sort(dataWangGSC$`Gene Symbol`)
Hits.Wang293FT <- sort(dataWang293FT$`Gene ID`)
Hits.Rother <- sort(dataRother$`Gene symbol`)
Hits.Shue <- sort(dataShue$genes[1:500])

papers <- c("Savidis", "Li", "Dukhovny", "WangGSC", 
            "Wang293FT", "Rother", "Shue")

genes.HDF <- c(Hits.Li, Hits.Rother, Hits.Dukhovny, Hits.Savidis, Hits.Shue, 
               Hits.WangGSC, Hits.Wang293FT)
genes.HDF <- genes.HDF |> unique() |> sort()

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_genes <- getBM(attributes=c('ensembl_gene_id', 
                                    'hgnc_symbol', 
                                    'external_gene_name', 
                                    'external_synonym'), mart = ensembl)

# write_csv(ensembl_genes, "ZIKVData/ensembl_genes.csv")

ensembl_ids <- mapIds(org.Hs.eg.db, 
                      keys = genes.HDF, 
                      column = "ENSEMBL", 
                      keytype = "ALIAS")

genes.HDF.DF <- tibble(
  genes = genes.HDF, 
  ensembl.ids = ensembl_ids
)

IDX <- is.na(genes.HDF.DF$ensembl.ids)
genes.HDF.DF$ensembl.ids[IDX] <- ensembl_genes$ensembl_gene_id[match(toupper(genes.HDF[IDX]), ensembl_genes$external_synonym)]
IDX <- is.na(genes.HDF.DF$ensembl.ids)
genes.HDF.DF$ensembl.ids[IDX] <- ensembl_genes$ensembl_gene_id[match(genes.HDF[IDX], ensembl_genes$external_gene_name)]
IDX <- is.na(genes.HDF.DF$ensembl.ids)
genes.HDF.DF$ensembl.ids[IDX] <- c("ENSG00000264066", NA, "ENSG00000271043", "ENSG00000213029")

rm(ensembl, ensembl_genes, ensembl_ids, IDX, genes.HDF)

genes.HDF.DF <- genes.HDF.DF |> 
  filter(!(is.na(ensembl.ids))) |> 
  mutate(Savidis = genes %in% Hits.Savidis, 
         Li = genes %in% Hits.Li, 
         Dukhovny = genes %in% Hits.Dukhovny, 
         WangGSC = genes %in% Hits.WangGSC, 
         Wang293FT = genes %in% Hits.Wang293FT, 
         Rother = genes %in% Hits.Rother, 
         Shue = genes %in% Hits.Shue) |> 
  rowwise() |> 
  mutate(Intersection = sum(c_across(3:9)))

# genes.HDF.DF$ensembl.ids[duplicated(genes.HDF.DF$ensembl.ids)]

genes.upset <- genes.HDF.DF |> 
  dplyr::select(ensembl.ids:Intersection) |> 
  group_by(ensembl.ids) |> 
  summarise_all(sum)

genes.upset[2:8] <- genes.upset  |> 
  dplyr::select(Savidis:Shue) |> 
  mutate_all(~ . == 1)


# Plotting the Upset plot -------------------------------------------------

upset(
  data = genes.upset,
  intersect = papers,
  name = 'Papers', 
  width_ratio=0.2,
  sort_intersections_by = c('degree', 'cardinality'),
  sort_intersections = 'ascending',
) +
  ggtitle("Genes from 6 studies")

# All Genes ---------------------------------------------------------------

gene.names <- c(dataSavidis$`Zika CRISPR Screen Hits - Top 100`, 
                dataLi$Gene, 
                dataDukhovny$OFFICIAL_SYMBOL, 
                dataWangGSC$`Gene Symbol`, 
                dataWang293FT$`Gene ID`, 
                dataRother$`Gene symbol`, 
                dataShue$genes) |> unique() |> sort()

genes <- data.frame(Symbol = gene.names, HDF = FALSE) |> 
  mutate(HDF = ifelse(Symbol %in% genes.HDF, TRUE, FALSE))

# Selecting Groups --------------------------------------------------------

GS3.of.7 <- (genes.upset |> 
  filter(Intersection >= 3))$ensembl.ids
GS2.of.7 <- (genes.upset |> 
  filter(Intersection == 2))$ensembl.ids
GS1.of.7 <- (genes.upset |> 
  filter(Intersection == 1))$ensembl.ids

# Gen set enrichment analysis ---------------------------------------------

# ALL GENES
gsea <- gost(query = genes, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaUp <- gsea$result |> filter(term_size < 500, term_size > 10)
# gostplot(gsea, interactive = FALSE)

gseaUp.mod <- gseaUp[, c("query", "source", "term_id", "term_name", "p_value", 
                         "query_size", "intersection_size", "term_size", 
                         "effective_domain_size", "intersection")]

gseaUp.mod$GeneRatio = paste0(gseaUp.mod$intersection_size,  "/", gseaUp.mod$query_size)
gseaUp.mod$BgRatio = paste0(gseaUp.mod$term_size, "/", gseaUp.mod$effective_domain_size)

names(gseaUp.mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                      "query_size", "Count", "term_size", "effective_domain_size", 
                      "geneID", "GeneRatio", "BgRatio")
gseaUp.mod$geneID = gsub(",", "/", gseaUp.mod$geneID)

# row.names(gp_mod) = gp_mod$ID

# write_xlsx(gseaUp.mod, "ZIKVData/GSEAup_Zika.xlsx")


dataCombined <- read.csv("ZIKVData/combinedFinal.csv", header=TRUE, 
                         stringsAsFactors=FALSE)



# Info --------------------------------------------------------------------

# Para calcular el 'inclusive_intersection' a partir del exclusive sólo sumar la
# cantidad de genes que coinciden en cada intersección
# Gene Ontology (GO): Biological Process (BP), Molecular Function (MF) and 
# Cellular Component (CC)
# Term Size: Genes que son del BP, MF o CC (Gene set)
# Query Size: Genes que reaccionaron independientemente del BP, MF o CC
# Count: Genes que son del BP, MF o CC y Que también reacionaron
# Effective domain size: Total de genes
#dataLi <- dataLi[order(dataLi$p.bh),]
