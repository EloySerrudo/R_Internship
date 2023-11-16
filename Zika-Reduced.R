# install.packages('/opt/gurobi1003/linux64/R/gurobi_10.0-3_R_4.2.0.tar.gz', repos=NULL)
# install.packages('slam')
# install.packages("~/Downloads/ontologyLIP_1.3.tar.gz", repos=NULL, type = "source")

library(gurobi)
library(ontologyLIP)


# FUNCTIONS ---------------------------------------------------------------

adaptGSEA <- function(DF, nameDF) {
  DF <- DF$result |> filter(term_size < 500, term_size > 10)
  DF.mod <- DF[,c("query", "source", "term_id", "term_name", "p_value", "query_size", 
                  "intersection_size", "term_size", "effective_domain_size", "intersection")]
  DF.mod$GeneRatio <- paste0(DF.mod$intersection_size,  "/", DF.mod$query_size)
  DF.mod$BgRatio <- paste0(DF.mod$term_size, "/", DF.mod$effective_domain_size)
  names(DF.mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust", 
                     "query_size", "Count", "term_size", "effective_domain_size", 
                     "geneID", "GeneRatio", "BgRatio")
  #DF.mod$geneID <- gsub(",", "/", DF.mod$geneID)
  DF.mod$ID <- paste0(DF.mod$ID, nameDF)
  DF.mod
}


removeRedundancy <- function(input, threshold) {
  input.short <- input
  input.short$geneID = gsub("/", ",", input.short$geneID)
  input.short <- input.short[,c("ID", "p.adjust", "geneID")]
  output.short <- removeRedundant(inputData = input.short, 
                                  jacCutoff = threshold, 
                                  model = "weighted",
                                  gprofilerOutput = F, 
                                  parallel = F)
  output <- input[input$ID %in% output.short$ID,]
  #output$geneID = gsub(",", "/", output$geneID)
  output
}


gmt.2.DataFrame <- function(gmt_file, go) {
  DF <- data.frame(ID = gmt_file$geneset.names, 
                   Category = go,
                   Descriptions = gmt_file$geneset.descriptions,
                   geneID = NA)
  for (i in 1:length(gmt_file$genesets)) {
    aux <- gmt_file$genesets[[i]]
    DF$geneID[i] <- paste(aux, collapse = ",")
  }
  rm(i, aux)
  DF
}

# Gen set enrichment analysis ---------------------------------------------

gseaSavidis <- gost(query = sort(Hits.Savidis), organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaLi <- gost(query = sort(Hits.Li), organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaDukhovny <- gost(query = sort(Hits.Dukhovny), organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaWangGSC <- gost(query = sort(Hits.WangGSC), organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaWang293FT <- gost(query = sort(Hits.Wang293FT), organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaRother <- gost(query = sort(Hits.Rother), organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaShue <- gost(query = sort(Hits.Shue), organism = "hsapiens", sources = c("GO"), evcodes = T)

gseaSavidis.mod <- adaptGSEA(gseaSavidis, "_Sa")
gseaLi.mod <- adaptGSEA(gseaLi, "_Li")
gseaDukhovny.mod <- adaptGSEA(gseaDukhovny, "_Du")
gseaWangGSC.mod <- adaptGSEA(gseaWangGSC, "_WG")
gseaWang293FT.mod <- adaptGSEA(gseaWang293FT, "_W2")
gseaRother.mod <- adaptGSEA(gseaRother, "_Ro")
gseaShue.mod <- adaptGSEA(gseaShue, "_Sh")


GSL <- rbind(
  gseaSavidis.mod[,c("ID", "p.adjust", "geneID")],
  gseaLi.mod[,c("ID", "p.adjust", "geneID")],
  gseaDukhovny.mod[,c("ID", "p.adjust", "geneID")],
  gseaWangGSC.mod[,c("ID", "p.adjust", "geneID")],
  gseaWang293FT.mod[,c("ID", "p.adjust", "geneID")],
  gseaRother.mod[,c("ID", "p.adjust", "geneID")],
  gseaShue.mod[,c("ID", "p.adjust", "geneID")]
)

# REMOVE REDUNDANCY -------------------------------------------------------

reduced.GSL <- removeRedundant(inputData = GSL, 
                                jacCutoff = 0.3, 
                                model = "weighted",
                                gprofilerOutput = F, 
                                parallel = F)
row.names(reduced.GSL) <- NULL

reduced.GSL <- reduced.GSL |> 
  mutate(origin = substr(ID, 12, 13), .after = ID) |> 
  mutate(ID = substr(ID, 1, 10))


# GMT FILE ----------------------------------------------------------------

install.packages("GSA")
library(GSA)

#https://baderlab.github.io/CBW_Pathways_2020/gprofiler-lab.html

GO.BP_path <- "ZIKVData/gprofiler_hsapiens.name/hsapiens.GO:BP.name.gmt"
GO.MF_path <- "ZIKVData/gprofiler_hsapiens.name/hsapiens.GO:MF.name.gmt"
GO.CC_path <- "ZIKVData/gprofiler_hsapiens.name/hsapiens.GO:CC.name.gmt"

gmt_GO.CC <- GSA.read.gmt(GO.CC_path); # 1976
gmt_GO.BP <- GSA.read.gmt(GO.BP_path); # 15704
gmt_GO.MF <- GSA.read.gmt(GO.MF_path); # 5035
rm(GO.BP_path, GO.MF_path, GO.CC_path)

gmt_GO.CC <- gmt.2.DataFrame(gmt_GO.CC, "GO:CC")
gmt_GO.BP <- gmt.2.DataFrame(gmt_GO.BP, "GO:BP")
gmt_GO.MF <- gmt.2.DataFrame(gmt_GO.MF, "GO:MF")

GO.Terms <- bind_rows(gmt_GO.CC, gmt_GO.BP, gmt_GO.MF) |> arrange(ID)
rm(gmt_GO.CC, gmt_GO.BP, gmt_GO.MF)

# write.csv(GO.Terms, file = "GO_Terms.csv", row.names = FALSE)

# NEW DATA ----------------------------------------------------------------

new.dataLi <- GO.Terms[(GO.Terms$ID %in% reduced.GSL$ID[reduced.GSL$origin == 'Li']), ]
rownames(new.dataLi) <- NULL
new.dataWangGSC <- GO.Terms[(GO.Terms$ID %in% reduced.GSL$ID[reduced.GSL$origin == 'WG']), ]
rownames(new.dataWangGSC) <- NULL
new.dataWang293FT <- GO.Terms[(GO.Terms$ID %in% reduced.GSL$ID[reduced.GSL$origin == 'W2']), ]
rownames(new.dataWang293FT) <- NULL
new.dataDukhovny <- GO.Terms[(GO.Terms$ID %in% reduced.GSL$ID[reduced.GSL$origin == 'Du']), ]
rownames(new.dataDukhovny) <- NULL
new.dataRother <- GO.Terms[(GO.Terms$ID %in% reduced.GSL$ID[reduced.GSL$origin == 'Ro']), ]
rownames(new.dataRother) <- NULL
new.dataShue <- GO.Terms[(GO.Terms$ID %in% reduced.GSL$ID[reduced.GSL$origin == 'Sh']), ]
rownames(new.dataShue) <- NULL

genes.Li <- paste(new.dataLi$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()

# Later -------------------------------------------------------------------

gseaSavidis.mod.red <- removeRedundancy(gseaSavidis.mod, 0.3)
gseaLi.mod.red <- removeRedundancy(gseaLi.mod, 0.3)
gseaDukhovny.mod.red <- removeRedundancy(gseaDukhovny.mod, 0.3)
gseaWangGSC.mod.red <- removeRedundancy(gseaWangGSC.mod, 0.3)
gseaWang293FT.mod.red <- removeRedundancy(gseaWang293FT.mod, 0.3)
gseaRother.mod.red <- removeRedundancy(gseaRother.mod, 0.3)
gseaShue.mod.red <- removeRedundancy(gseaShue.mod, 0.3)

reduced.GSL <- sort(c(
  gseaSavidis.mod.red$ID, 
  gseaLi.mod.red$ID, 
  gseaDukhovny.mod.red$ID, 
  gseaWangGSC.mod.red$ID, 
  gseaWang293FT.mod.red$ID, 
  gseaRother.mod.red$ID, 
  gseaUpShue.mod.red$ID
))
unique(reduced.GSL)



str_split(reduced.GSL2$geneID[1], ",")[[1]]

reducedGeneSet <- data.frame(genes, 
                             genes %in% str_split(reduced.GSL2$geneID[1], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[2], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[3], ",")[[1]],
                             genes %in% str_split(reduced.GSL2$geneID[4], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[5], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[6], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[7], ",")[[1]],
                             genes %in% str_split(reduced.GSL2$geneID[8], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[9], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[10], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[11], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[12], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[13], ",")[[1]],
                             genes %in% str_split(reduced.GSL2$geneID[14], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[15], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[16], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[17], ",")[[1]],
                             genes %in% str_split(reduced.GSL2$geneID[18], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[19], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[20], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[21], ",")[[1]], 
                             genes %in% str_split(reduced.GSL2$geneID[22], ",")[[1]],
                             genes %in% str_split(reduced.GSL2$geneID[23], ",")[[1]], 
                             row.names = NULL)

colnames(reducedGeneSet)[-1] <- substr(reduced.GSL2$ID, 1, 10)
gene.sets.2 <- colnames(reducedGeneSet)[-1]

upset(
  data = reducedGeneSet,
  intersect = gene.sets.2,
  name = 'Reduced Gene Sets', 
  width_ratio=0.3, 
  min_degree=1,
  sort_intersections_by = 'degree',
  sort_intersections = 'ascending'
)


gseaUpShue.mod.red <- gseaUpShue.mod[gseaUpShue.mod$ID %in% output$ID,]

gseaUpShue.mod$geneID = gsub(",", "/", gseaUpShue.mod$geneID)
gseaUpShue.mod.red$geneID = gsub(",", "/", gseaUpShue.mod.red$geneID)

write_xlsx(gseaUpShue.mod, "ZIKVData/GSEAupShue_Zika.xlsx")
write_xlsx(gseaUpShue.mod.red, "ZIKVData/GSEAupShueRed_Zika.xlsx")
