library(tidyverse)
library(readxl)
library(writexl)
library(ComplexUpset)
library(gprofiler2)


# Loading Datasets --------------------------------------------------------

dataSavidis <- read_excel("ZKVData/1-s2.0-S2211124716307689-mmc8.xlsx", 
                          sheet = "Summary CME Compare")[-1,]

dataLi <- read_excel("ZKVData/pnas.1900867116.sd01.xlsx", sheet = "zika_gw_5th_best")
colnames(dataLi) <- as.vector(dataLi[1,])
dataLi <- dataLi[-1,]
cols.num <- c("sgRNAs", "zika.init.5th", "p.bh")
dataLi[cols.num] <- sapply(dataLi[cols.num],as.numeric)
dataLi <- dataLi[order(dataLi$p.bh),]

dataDukhovny <- read.csv("ZKVData/JVI.00211-19-sd003.csv", header=TRUE, 
                         sep = "\t", stringsAsFactors=FALSE)
dataDukhovny <- dataDukhovny |> 
  mutate(Gene = ifelse(grepl("^[A-Z0-9]+$", id), id, str_extract(id, "\\(([^)]+)\\)$")), 
         .after = id) |> 
  mutate(Gene = gsub("\\(|\\)", "", Gene)) |> 
  arrange(pos.rank) |> 
  filter(!is.na(Gene))

dataWangGSC <- read_excel("ZKVData/NIHMS1553325-supplement-2.xlsx", 
                          sheet = "Ranking", col_names = FALSE, skip = 1)
dataWang293FT <- read_excel("ZKVData/NIHMS1553325-supplement-3.xlsx", 
                            sheet = "Sheet1", col_names = FALSE, skip = 2)
colnames(dataWangGSC) <- as.vector(dataWangGSC[1,])
colnames(dataWang293FT) <- as.vector(dataWang293FT[1,])
dataWangGSC <- dataWangGSC[-1,]
dataWang293FT <- dataWang293FT[-1,]
dataWang293FT[[1,1]] <- "MMGT1"

dataRother <- read_excel("ZKVData/1-s2.0-S0168170221000459-mmc1.xlsx", 
                         sheet = "(B) Gene ranking")

dataShue <- read_excel("ZKVData/jvi.00596-21-s0001.xls", sheet = "Genelist")
cols.num <- c("deseq2.FC", "deseq2.pval", "mageck.rank.pos", "mageck.fdr.pos", 
              "mageck.rank.neg", "mageck.fdr.neg")
dataShue[cols.num] <- sapply(dataShue[cols.num], as.numeric)
dataShue <- dataShue |> arrange(deseq2.pval, desc(deseq2.FC))
rm(cols.num)


# Extracting Hits ---------------------------------------------------------


Hits.Savidis <- dataSavidis$`Zika CRISPR Screen Hits - Top 100`
Hits.Li <- dataLi$Gene[1:500]
Hits.Dukhovny <- dataDukhovny[1:500, 2]
Hits.WangGSC <- dataWangGSC$`Gene Symbol`
Hits.Wang293FT <- dataWang293FT$`Gene ID`
Hits.Rother <- dataRother$`Gene symbol`
Hits.Shue <- dataShue$genes[1:500]

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
colnames(genes.DF)[-1] <- c("Wang.GSC", "Wang.293FT", "Shue", "Savidis", 
                            "Rother", "Li", "Dukhovny")
papers <- colnames(genes.DF)[-1]


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
  filter(Li & Wang.293FT & Rother & Savidis & Shue) |> 
  select(genes) |> 
  unlist(use.names = FALSE)

# Gen set enrichment analysis ---------------------------------------------

# ALL GENES
gsea <- gost(query = genes.int, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaUp <- gsea$result |> filter(term_size < 500, term_size > 10)
gostplot(gsea, interactive = FALSE)

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

write_xlsx(gseaUp.mod, "ZIKVData/GSEAup_Zika.xlsx")


# SAVIDIS DATA -----------------------------------------------------------

gseaSavidis <- gost(query = Hits.Savidis, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaUpSavidis <- gseaSavidis$result |> filter(term_size < 500, term_size > 10)
gseaUpSavidis.mod <- gseaUpSavidis[,c("query", "source", "term_id", "term_name", "p_value", 
                                "query_size", "intersection_size", "term_size", 
                                "effective_domain_size", "intersection")]

gseaUpSavidis.mod$GeneRatio = paste0(gseaUpSavidis.mod$intersection_size,  "/", gseaUpSavidis.mod$query_size)
gseaUpSavidis.mod$BgRatio = paste0(gseaUpSavidis.mod$term_size, "/", gseaUpSavidis.mod$effective_domain_size)

names(gseaUpSavidis.mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                          "query_size", "Count", "term_size", "effective_domain_size", 
                          "geneID", "GeneRatio", "BgRatio")
gseaUpSavidis.mod$geneID = gsub(",", "/", gseaUpSavidis.mod$geneID)


SavidisGeneSets <- data.frame(genes, 
                           genes %in% str_split(gseaUpSavidis.mod$geneID[1], "/")[[1]], 
                           genes %in% str_split(gseaUpSavidis.mod$geneID[2], "/")[[1]], 
                           genes %in% str_split(gseaUpSavidis.mod$geneID[3], "/")[[1]],
                           genes %in% str_split(gseaUpSavidis.mod$geneID[4], "/")[[1]], 
                           genes %in% str_split(gseaUpSavidis.mod$geneID[5], "/")[[1]], 
                           row.names = NULL)
colnames(SavidisGeneSets)[-1] <- paste0("Savidis_", gseaUpSavidis.mod$ID[1:11])

gene.sets <- colnames(SavidisGeneSets)[-1]

upset(
  data = SavidisGeneSets,
  intersect = gene.sets,
  name = 'Savidis Gene Sets', 
  width_ratio=0.3, 
  min_degree=1
) +
  ggtitle("Pathways Savidis (2016). Biological Processes")


write_xlsx(gseaUpShue.mod, "ZIKVData/GSEAupShue_Zika.xlsx")





# SHUE DATA -----------------------------------------------------------

gseaShue <- gost(query = Hits.Shue, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaUpShue <- gseaShue$result |> filter(term_size < 500, term_size > 10)
gseaUpShue.mod <- gseaUpShue[,c("query", "source", "term_id", "term_name", "p_value", 
                                "query_size", "intersection_size", "term_size", 
                                "effective_domain_size", "intersection")]

gseaUpShue.mod$GeneRatio = paste0(gseaUpShue.mod$intersection_size,  "/", gseaUpShue.mod$query_size)
gseaUpShue.mod$BgRatio = paste0(gseaUpShue.mod$term_size, "/", gseaUpShue.mod$effective_domain_size)

names(gseaUpShue.mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                          "query_size", "Count", "term_size", "effective_domain_size", 
                          "geneID", "GeneRatio", "BgRatio")
gseaUpShue.mod$geneID = gsub(",", "/", gseaUpShue.mod$geneID)


ShueGeneSets <- data.frame(genes, 
                       genes %in% str_split(gseaUpShue.mod$geneID[1], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[2], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[3], "/")[[1]],
                       genes %in% str_split(gseaUpShue.mod$geneID[4], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[5], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[6], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[7], "/")[[1]],
                       genes %in% str_split(gseaUpShue.mod$geneID[8], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[9], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[10], "/")[[1]], 
                       genes %in% str_split(gseaUpShue.mod$geneID[11], "/")[[1]],
                       row.names = NULL)
colnames(ShueGeneSets)[-1] <- paste0("Shue_", gseaUpShue.mod$ID[1:11])

gene.sets <- colnames(ShueGeneSets)[-1]

upset(
  data = ShueGeneSets,
  intersect = gene.sets,
  name = 'Shue Gene Sets', 
  width_ratio=0.3, 
  mode='inclusive_intersection',
  min_degree=1
) +
  ggtitle("Pathways Shue (2021). Biological Processes")


write_xlsx(gseaUpShue.mod, "ZIKVData/GSEAupShue_Zika.xlsx")



# Later -------------------------------------------------------------------


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE"), force = TRUE)

library(clusterProfiler)
library(enrichplot)
library(DOSE)

# modify the g:Profiler data frame


# define as compareClusterResult object
gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
enrichplot::dotplot(gp_mod_cluster)

# define as enrichResult object
gp_mod_enrich  = new("enrichResult", result = gp_mod)

barplot(gp_mod_enrich, showCategory = 40, font.size = 16) + 
  ggplot2::facet_grid(~Cluster) +
  ggplot2::ylab("Intersection size")

# Gene Ontology (GO): Biological Process (BP), Molecular Function (MF) and Cellular Component (CC)
# Term Size: Genes que son del BP, MF o CC (Gene set)
# Query Size: Genes que reaccionaron independientemente del BP, MF o CC
# Count: Genes que son del BP, MF o CC y Que también reacionaron
# Effective domain size: Total de genes


# LI DATA ---------------------------------------------------------------

gseaLi <- gost(query = Hits.Li, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaUpLi <- gseaLi$result |> filter(term_size < 500, term_size > 10)
gseaUpLi.mod <- gseaUpLi[,c("query", "source", "term_id", "term_name", "p_value", 
                            "query_size", "intersection_size", "term_size", 
                            "effective_domain_size", "intersection")]

gseaUpLi.mod$GeneRatio = paste0(gseaUpLi.mod$intersection_size,  "/", gseaUpLi.mod$query_size)
gseaUpLi.mod$BgRatio = paste0(gseaUpLi.mod$term_size, "/", gseaUpLi.mod$effective_domain_size)

names(gseaUpLi.mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                        "query_size", "Count", "term_size", "effective_domain_size", 
                        "geneID", "GeneRatio", "BgRatio")
gseaUpLi.mod$geneID = gsub(",", "/", gseaUpLi.mod$geneID)

LiGeneSets <- gseaUpLi.mod[gseaUpLi.mod$Category == "GO:BP", 10]

a <- str_split(LiGeneSets, "/")
b <- genes %in% a$value

c <- data.frame(genes)
c <- c |> mutate()

LiGeneSets <- data.frame(genes, 
                         genes %in% str_split(gseaUpSavidis.mod$geneID[1], "/")[[1]], 
                         genes %in% str_split(gseaUpSavidis.mod$geneID[2], "/")[[1]], 
                         genes %in% str_split(gseaUpSavidis.mod$geneID[3], "/")[[1]],
                         genes %in% str_split(gseaUpSavidis.mod$geneID[4], "/")[[1]], 
                         genes %in% str_split(gseaUpSavidis.mod$geneID[5], "/")[[1]], 
                         row.names = NULL)
colnames(SavidisGeneSets)[-1] <- paste0("Savidis_", gseaUpSavidis.mod$ID[1:11])

gene.sets <- colnames(SavidisGeneSets)[-1]

upset(
  data = SavidisGeneSets,
  intersect = gene.sets,
  name = 'Savidis Gene Sets', 
  width_ratio=0.3, 
  min_degree=1
) +
  ggtitle("Pathways Savidis (2016). Biological Processes")


write_xlsx(gseaUpShue.mod, "ZIKVData/GSEAupShue_Zika.xlsx")
