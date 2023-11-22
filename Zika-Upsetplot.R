library(tidyverse)
library(readxl)
library(writexl)
library(ComplexUpset)
library(gprofiler2)


# Loading Datasets --------------------------------------------------------

dataSavidis <- read_excel("1-s2.0-S2211124716307689-mmc8.xlsx", 
                          sheet = "Summary CME Compare")[-1,]

dataLi <- read_excel("pnas.1900867116.sd01.xlsx", sheet = "zika_gw_5th_best")
colnames(dataLi) <- as.vector(dataLi[1,])
dataLi <- dataLi[-1,]
cols.num <- c("sgRNAs", "zika.init.5th", "p.bh")
dataLi[cols.num] <- sapply(dataLi[cols.num],as.numeric)
dataLi <- dataLi[order(dataLi$p.bh),]

dataDukhovny <- read.csv("inData/JVI.00211-19-sd003.csv", header=TRUE, 
                         sep = "\t", stringsAsFactors=FALSE)
dataDukhovny <- dataDukhovny |> 
  mutate(Gene = ifelse(grepl("^[A-Z0-9]+$", id), id, str_extract(id, "\\(([^)]+)\\)$")), 
         .after = id) |> 
  mutate(Gene = gsub("\\(|\\)", "", Gene)) |> 
  arrange(neg.rank) |> 
  filter(!is.na(Gene))

dataWangGSC <- read_excel("NIHMS1553325-supplement-2.xlsx", 
                          sheet = "Ranking", col_names = FALSE, skip = 1)
dataWang293FT <- read_excel("NIHMS1553325-supplement-3.xlsx", 
                            sheet = "Sheet1", col_names = FALSE, skip = 2)
colnames(dataWangGSC) <- as.vector(dataWangGSC[1,])
colnames(dataWang293FT) <- as.vector(dataWang293FT[1,])
dataWangGSC <- dataWangGSC[-1,]
dataWang293FT <- dataWang293FT[-1,]
dataWang293FT[[1,1]] <- "MMGT1"

dataRother <- read_excel("1-s2.0-S0168170221000459-mmc1.xlsx", 
                         sheet = "(B) Gene ranking")

dataShue <- read_excel("jvi.00596-21-s0001.xls", sheet = "Genelist")
cols.num <- c("deseq2.FC", "deseq2.pval", "mageck.rank.pos", "mageck.fdr.pos", 
              "mageck.rank.neg", "mageck.fdr.neg")
dataShue[cols.num] <- sapply(dataShue[cols.num], as.numeric)
dataShue <- dataShue |> arrange(deseq2.pval, desc(deseq2.FC))
rm(cols.num)


# Extracting Hits ---------------------------------------------------------


Hits.Savidis <- sort(dataSavidis$`Zika CRISPR Screen Hits - Top 100`)
Hits.Li <- sort(dataLi$Gene[1:500])
Hits.Dukhovny <- sort(dataDukhovny[1:500, 2])
Hits.WangGSC <- sort(dataWangGSC$`Gene Symbol`)
Hits.Wang293FT <- sort(dataWang293FT$`Gene ID`)
Hits.Rother <- sort(dataRother$`Gene symbol`)
Hits.Shue <- sort(dataShue$genes[1:500])

papers <- c("Wang.GSC", "Wang.293FT", "Shue", "Savidis", 
            "Rother", "Li", "Dukhovny")

genes <- c(Hits.Li, Hits.Rother, Hits.Dukhovny, Hits.Savidis, Hits.Shue, 
           Hits.WangGSC, Hits.Wang293FT)
genes <- genes |> unique() |> sort()

genes.DF <- data.frame(genes, 
                       genes %in% Hits.WangGSC, 
                       genes %in% Hits.Wang293FT, 
                       genes %in% Hits.Shue, 
                       genes %in% Hits.Savidis, 
                       genes %in% Hits.Rother,
                       genes %in% Hits.Li, 
                       genes %in% Hits.Dukhovny,
                       row.names = NULL)
colnames(genes.DF)[-1] <- papers


# Plotting the Upset plot -------------------------------------------------

upset(
  data = genes.DF,
  intersect = papers,
  name = 'Papers', 
  width_ratio=0.2,
  sort_intersections_by = 'degree',
  sort_intersections = 'ascending',
  sort_sets = FALSE
) +
  ggtitle("Genes from 6 studies")


# Selecting Groups --------------------------------------------------------

genes.DF.int <- genes.DF |> 
  rowwise() |>
  mutate(Intersection = sum(c_across(2:8))) |> 
  filter(Intersection > 1)

genes.int <- genes.DF.int$genes

genes.DF.int |> 
  filter(Li & Wang.293FT & Rother & Savidis & Shue & Wang.GSC) |> 
  select(genes) |> 
  unlist(use.names = FALSE)

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
