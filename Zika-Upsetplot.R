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


Hits.Savidis <- dataSavidis[,-2]
Hits.Li <- dataLi[1:500, -(2:4)]
Hits.Dukhovny <- dataDukhovny[1:500, -c(1, 3:15)]
Hits.WangGSC <- dataWangGSC[,-(2:3)]
Hits.Wang293FT <- dataWang293FT[,-(2:5)]
Hits.Rother <- dataRother[,-(2:6)]
Hits.Shue <- dataShue[1:500,-(2:7)]

genes <- c(unlist(Hits.Li), unlist(Hits.Rother), Hits.Dukhovny, unlist(Hits.Savidis), 
           unlist(Hits.Shue), unlist(Hits.WangGSC), unlist(Hits.Wang293FT))
genes <- genes |> unique() |> sort()

genes.DF <- data.frame(genes, 
                       genes %in% unlist(Hits.Savidis), 
                       genes %in% unlist(Hits.Li), 
                       genes %in% Hits.Dukhovny,
                       genes %in% unlist(Hits.WangGSC), 
                       genes %in% unlist(Hits.Wang293FT), 
                       genes %in% unlist(Hits.Rother), 
                       genes %in% unlist(Hits.Shue),
                       row.names = NULL)
colnames(genes.DF)[-1] <- c("Savidis", "Li", "Dukhovny", "Wang.GSC", 
                            "Wang.293FT", "Rother", "Shue")
papers <- colnames(genes.DF)[-1]


# Plotting the Upset plot -------------------------------------------------

upset(
  data = genes.DF,
  intersect = papers,
  name = 'Papers', 
  width_ratio=0.3, 
  mode='inclusive_intersection',
  min_degree=2
)


# Selecting Groups --------------------------------------------------------

genes.DF.int <- genes.DF |> 
  rowwise() |>
  mutate(Intersection = sum(c_across(2:8))) |> 
  filter(Intersection > 1)

genes.int <- unlist(genes.DF.int$genes, use.names = FALSE)

genes.DF.int |> 
  filter(Li & Wang.293FT & Rother & Savidis & Shue) |> 
  select(genes) |> 
  unlist(use.names = FALSE)

# Gen set enrichment analysis ---------------------------------------------

# ALL GENES
gsea <- gost(query = genes.int, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaUp <- gsea$result |> filter(term_size < 500, term_size > 10)
write_xlsx(gseaUp, "ZIKVData/GSEAup_Zika.xlsx")

gostplot(gsea, interactive = FALSE)

GO.MF <- gseaUp$intersection[20:22]

# SHUE GENES
gseaShue <- gost(query = Hits.Shue$genes, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaUpShue <- gseaShue$result |> filter(term_size < 500, term_size > 10)





GO.BP <- unique(str_split(paste(gseaUp$intersection[1:12], collapse = ","), ",")[[1]])
GO.BP

GO.CC <- unique(str_split(paste(gseaUp$intersection[13:19], collapse = ","), ",")[[1]])
GO.CC

GO.MF <- unique(str_split(paste(gseaUp$intersection[20:22], collapse = ","), ",")[[1]])
GO.MF

genes2 <- sort(unique(c(GO.BP, GO.CC, GO.MF)))

genes2.DF <- data.frame(genes2, 
                       genes2 %in% GO.BP, 
                       genes2 %in% GO.CC, 
                       genes2 %in% GO.MF,
                       row.names = NULL)
colnames(genes2.DF)[-1] <- c("GO.BP", "GO.CC", "GO.MF")
gen.set <- colnames(genes2.DF)[-1]

upset(
  data = genes2.DF,
  intersect = gen.set,
  name = 'gen.set', 
  width_ratio=0.3, 
  mode='inclusive_intersection',
)

BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE"))

library(clusterProfiler)
library(enrichplot)
library(DOSE)

# modify the g:Profiler data frame
gp_mod <- gseaUp[,c("query", "source", "term_id",
                            "term_name", "p_value", "query_size", 
                            "intersection_size", "term_size", 
                            "effective_domain_size", "intersection")]
gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)

names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                  "query_size", "Count", "term_size", "effective_domain_size", 
                  "geneID", "GeneRatio", "BgRatio")
gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
#row.names(gp_mod) = gp_mod$ID

# define as compareClusterResult object
gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)

# define as enrichResult object
gp_mod_enrich  = new("enrichResult", result = gp_mod)
