library(gurobi)
library(ontologyLIP)


# FUNCTIONS ---------------------------------------------------------------

adaptGSEA <- function(DF) {
  DF <- DF$result |> filter(term_size < 500, term_size > 10)
  
  DF.mod <- DF[,c("query", "source", "term_id", "term_name", "p_value", "query_size", 
                  "intersection_size", "term_size", "effective_domain_size", "intersection")]
  DF.mod$GeneRatio <- paste0(DF.mod$intersection_size,  "/", DF.mod$query_size)
  DF.mod$BgRatio <- paste0(DF.mod$term_size, "/", DF.mod$effective_domain_size)
  names(DF.mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust", 
                     "query_size", "Count", "term_size", "effective_domain_size", 
                     "geneID", "GeneRatio", "BgRatio")
  #DF.mod$geneID <- gsub(",", "/", DF.mod$geneID)
  DF.mod$ID <- paste0(DF.mod$ID)
  DF.mod
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

bg.Li <- sort(dataLi$Gene)
bg.Dukhovny <- sort(dataDukhovny$Gene)
bg.Shue <- sort(dataShue$genes)

sum(is.na(bg.Dukhovny))


gseaSavidis <- gost(query = Hits.Savidis, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaLi <- gost(query = Hits.Li, organism = "hsapiens", sources = c("GO"), evcodes = T, custom_bg = bg.Li)

gseaDukhovny <- gost(query = Hits.Dukhovny, organism = "hsapiens", sources = c("GO"), evcodes = T)#, custom_bg = bg.Dukhovny)

gseaWangGSC <- gost(query = Hits.WangGSC, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaWang293FT <- gost(query = Hits.Wang293FT, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaRother <- gost(query = Hits.Rother, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaShue <- gost(query = Hits.Shue, organism = "hsapiens", sources = c("GO"), evcodes = T, custom_bg = bg.Shue)

gseaSavidis.mod <- adaptGSEA(gseaSavidis)
gseaLi.mod <- adaptGSEA(gseaLi)
gseaDukhovny.mod <- adaptGSEA(gseaDukhovny)
gseaWangGSC.mod <- adaptGSEA(gseaWangGSC)
gseaWang293FT.mod <- adaptGSEA(gseaWang293FT)
gseaRother.mod <- adaptGSEA(gseaRother)
gseaShue.mod <- adaptGSEA(gseaShue)

rm(gseaSavidis, gseaLi, gseaDukhovny, gseaWangGSC, 
   gseaWang293FT, gseaRother, gseaShue)

GSL.1 <- rbind(
  gseaSavidis.mod[,c("ID", "p.adjust", "geneID")],
  gseaLi.mod[,c("ID", "p.adjust", "geneID")],
  gseaDukhovny.mod[,c("ID", "p.adjust", "geneID")],
  gseaWangGSC.mod[,c("ID", "p.adjust", "geneID")],
  gseaWang293FT.mod[,c("ID", "p.adjust", "geneID")],
  gseaRother.mod[,c("ID", "p.adjust", "geneID")],
  gseaShue.mod[,c("ID", "p.adjust", "geneID")]
)

GSL.2 <- GSL.1[match(unique(GSL.1$ID), GSL.1$ID), ]

pre.gene.sets <- sort(GSL.2$ID) |> sort()

# PRE ANALYSIS ------------------------------------------------------------

pre.gene.sets.DF <- data.frame(pre.gene.sets, 
                               pre.gene.sets %in% gseaWangGSC.mod$ID, 
                               pre.gene.sets %in% gseaWang293FT.mod$ID, 
                               pre.gene.sets %in% gseaShue.mod$ID, 
                               pre.gene.sets %in% gseaSavidis.mod$ID, 
                               pre.gene.sets %in% gseaRother.mod$ID,
                               pre.gene.sets %in% gseaLi.mod$ID, 
                               pre.gene.sets %in% gseaDukhovny.mod$ID,
                               row.names = NULL)
colnames(pre.gene.sets.DF)[-1] <- c("Wang.GSC", "Wang.293FT", "Shue", "Savidis",
                                    "Rother", "Li", "Dukhovny")

# REMOVE REDUNDANCY -------------------------------------------------------

reduced.GSL.1 <- removeRedundant(inputData = GSL.1, 
                                   jacCutoff = 0.3, 
                                   model = "weighted",
                                   gprofilerOutput = F, 
                                   parallel = F)
row.names(reduced.GSL.1) <- NULL

reduced.GSL.2 <- removeRedundant(inputData = GSL.2, 
                                   jacCutoff = 0.3, 
                                   model = "weighted",
                                   gprofilerOutput = F, 
                                   parallel = F)
row.names(reduced.GSL.2) <- NULL

# Al final si tomamos los únicos de reduced.GSL.1 se parece a reduced.GSL.2
# porque toma el primero más largo

# GMT FILE ----------------------------------------------------------------

#install.packages("GSA")
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
GO.Terms <- read.csv("ZIKVData/GO_Terms.csv")


# Build a new GMT File ----------------------------------------------------

GO.Terms.Reduced <- GO.Terms[GO.Terms$ID %in% reduced.GSL.2$ID, ]
rownames(GO.Terms.Reduced) <- NULL
write.csv(GO.Terms.Reduced, file = "ZIKVData/GO_Terms_Reduced.csv", row.names = FALSE)

new_gmt <- apply(GO.Terms.Reduced[, -2], 1, function(row) {
  paste(c(row[1], row[2], unlist(strsplit(as.character(row[3]), ","))), collapse = "\t")
})
writeLines(new_gmt, con = "ZIKVData/new_gmt.name.gmt")

# POST ANALYSIS -----------------------------------------------------------
# gp__EhA0_PFup_VAE
new.gseaSavidis <- gost(query = Hits.Savidis, organism = "gp__avQf_XJ5q_AEk", sources = c("GO"), evcodes = T)
new.gseaLi <- gost(query = Hits.Li, organism = "gp__avQf_XJ5q_AEk", sources = c("GO"), evcodes = T, custom_bg = bg.Li)
new.gseaDukhovny <- gost(query = Hits.Dukhovny, organism = "gp__avQf_XJ5q_AEk", sources = c("GO"), evcodes = T)
new.gseaWangGSC <- gost(query = Hits.WangGSC, organism = "gp__avQf_XJ5q_AEk", sources = c("GO"), evcodes = T)
new.gseaWang293FT <- gost(query = Hits.Wang293FT, organism = "gp__avQf_XJ5q_AEk", sources = c("GO"), evcodes = T)
new.gseaRother <- gost(query = Hits.Rother, organism = "gp__avQf_XJ5q_AEk", sources = c("GO"), evcodes = T)
new.gseaShue <- gost(query = Hits.Shue, organism = "gp__avQf_XJ5q_AEk", sources = c("GO"), evcodes = T, custom_bg = bg.Shue)
rm(Hits.Dukhovny, Hits.Li, Hits.Rother, Hits.Wang293FT, Hits.WangGSC, Hits.Savidis, Hits.Shue)

new.gseaSavidis.mod <- adaptGSEA(new.gseaSavidis) |> mutate(Origin = "Savidis")
new.gseaLi.mod <- adaptGSEA(new.gseaLi) |> mutate(Origin = "Li")
new.gseaDukhovny.mod <- adaptGSEA(new.gseaDukhovny) |> mutate(Origin = "Dukhovny")
new.gseaWangGSC.mod <- adaptGSEA(new.gseaWangGSC) |> mutate(Origin = "WangGSC")
new.gseaWang293FT.mod <- adaptGSEA(new.gseaWang293FT) |> mutate(Origin = "Wang293FT")
new.gseaRother.mod <- adaptGSEA(new.gseaRother) |> mutate(Origin = "Rother")
new.gseaShue.mod <- adaptGSEA(new.gseaShue) |> mutate(Origin = "Shue")
rm(new.gseaSavidis, new.gseaLi, new.gseaDukhovny, new.gseaWangGSC, 
   new.gseaWang293FT, new.gseaRother, new.gseaShue)

GSL.Origin <- rbind(new.gseaSavidis.mod[c(3,13)],
           new.gseaLi.mod[c(3,13)],
           new.gseaDukhovny.mod[c(3,13)],
           new.gseaWangGSC.mod[c(3,13)],
           new.gseaWang293FT.mod[c(3,13)],
           new.gseaRother.mod[c(3,13)],
           new.gseaShue.mod[c(3,13)])

pos.gene.sets <- c(
  new.gseaSavidis.mod$ID,
  new.gseaLi.mod$ID,
  new.gseaDukhovny.mod$ID,
  new.gseaWangGSC.mod$ID,
  new.gseaWang293FT.mod$ID,
  new.gseaRother.mod$ID,
  new.gseaShue.mod$ID
) |> unique() |> sort()

pos.gene.sets.DF <- data.frame(pos.gene.sets, 
                           pos.gene.sets %in% new.gseaWangGSC.mod$ID, 
                           pos.gene.sets %in% new.gseaWang293FT.mod$ID, 
                           pos.gene.sets %in% new.gseaShue.mod$ID, 
                           pos.gene.sets %in% new.gseaSavidis.mod$ID, 
                           pos.gene.sets %in% new.gseaRother.mod$ID,
                           pos.gene.sets %in% new.gseaLi.mod$ID, 
                           pos.gene.sets %in% new.gseaDukhovny.mod$ID,
                           row.names = NULL)
colnames(pos.gene.sets.DF)[-1] <- c("Wang.GSC", "Wang.293FT", "Shue", "Savidis", 
                                "Rother", "Li", "Dukhovny")




# Plotting the Upset plot -------------------------------------------------

upset(
  data = pre.gene.sets.DF,
  intersect = papers,
  name = 'Papers', 
  width_ratio=0.2,
  sort_intersections_by = c('degree', 'cardinality'),
  sort_intersections = 'ascending',
  sort_sets = 'descending'
) +
  ggtitle("Pre Reduction Gene Sets")

upset(
  data = pos.gene.sets.DF,
  intersect = papers,
  name = 'Papers', 
  width_ratio=0.2,
  sort_intersections_by = c('degree', 'cardinality'),
  sort_intersections = 'ascending',
  sort_sets = 'descending'
) +
  ggtitle("Reduced Gene Sets")



Reduced.Gene.Sets <- GO.Terms |> filter(ID %in% pos.gene.sets) |> mutate(Origin = NA)

for (i in 1:nrow(Reduced.Gene.Sets) ) {
  Reduced.Gene.Sets[i, 5] <- GSL.Origin |> filter(ID == Reduced.Gene.Sets[i, 1]) |> select(Origin) |> unlist() |> paste(collapse = "/")
}

# Later -------------------------------------------------------------------

GO.set <- z.gseaSavidis.mod$ID[z.gseaSavidis.mod$ID %in% z.reduced.GSL.2$ID]
z.new.gseaSavidis.mod <- GO.Terms[GO.Terms$ID %in% GO.set, ]
rownames(z.new.gseaSavidis.mod) <- NULL
GO.set <- z.gseaLi.mod$ID[z.gseaLi.mod$ID %in% z.reduced.GSL.2$ID]
z.new.gseaLi.mod <- GO.Terms[GO.Terms$ID %in% GO.set, ]
rownames(z.new.gseaLi.mod) <- NULL
z.new.gseaDukhovny.mod <- GO.Terms[GO.Terms$ID %in% z.gseaDukhovny.mod$ID[z.gseaDukhovny.mod$ID %in% z.reduced.GSL.2$ID], ]
rownames(z.new.gseaDukhovny.mod) <- NULL
z.new.gseaWangGSC.mod <- GO.Terms[GO.Terms$ID %in% z.gseaWangGSC.mod$ID[z.gseaWangGSC.mod$ID %in% z.reduced.GSL.2$ID], ]
rownames(z.new.gseaWangGSC.mod) <- NULL
z.new.gseaWang293FT.mod <- GO.Terms[GO.Terms$ID %in% z.gseaWang293FT.mod$ID[z.gseaWang293FT.mod$ID %in% z.reduced.GSL.2$ID], ]
rownames(z.new.gseaWang293FT.mod) <- NULL
z.new.gseaRother.mod <- GO.Terms[GO.Terms$ID %in% z.gseaRother.mod$ID[z.gseaRother.mod$ID %in% z.reduced.GSL.2$ID], ]
rownames(z.new.gseaRother.mod) <- NULL
z.new.gseaShue.mod <- GO.Terms[GO.Terms$ID %in% z.gseaShue.mod$ID[z.gseaShue.mod$ID %in% z.reduced.GSL.2$ID], ]
rownames(z.new.gseaShue.mod) <- NULL

GO.Terms.Reduced <- GO.Terms[GO.Terms$ID %in% reduced.GSL.2$ID, ]
rownames(GO.Terms.Reduced) <- NULL
new.gseaSavidis.mod$Category <- GO.Terms.Reduced$Category[match(new.gseaSavidis.mod$ID, GO.Terms.Reduced$ID)]
new.gseaLi.mod$Category <- GO.Terms.Reduced$Category[match(new.gseaLi.mod$ID, GO.Terms.Reduced$ID)]
new.gseaDukhovny.mod$Category <- GO.Terms.Reduced$Category[match(new.gseaDukhovny.mod$ID, GO.Terms.Reduced$ID)]
new.gseaWangGSC.mod$Category <- GO.Terms.Reduced$Category[match(new.gseaWangGSC.mod$ID, GO.Terms.Reduced$ID)]
new.gseaWang293FT.mod$Category <- GO.Terms.Reduced$Category[match(new.gseaWang293FT.mod$ID, GO.Terms.Reduced$ID)]
new.gseaRother.mod$Category <- GO.Terms.Reduced$Category[match(new.gseaRother.mod$ID, GO.Terms.Reduced$ID)]
new.gseaShue.mod$Category <- GO.Terms.Reduced$Category[match(new.gseaShue.mod$ID, GO.Terms.Reduced$ID)]
rm(GO.Terms, GO.Terms.Reduced)

pos.genes.DF |> 
  rowwise() |>
  mutate(Intersection = sum(c_across(2:8))) |> 
  filter(Li & Wang.293FT & Rother & Savidis & Shue)

GO.Terms.Reduced[GO.Terms.Reduced$ID == pos.genes.DF |> 
                   rowwise() |>
                   mutate(Intersection = sum(c_across(2:8))) |> 
                   filter(Li & Wang.293FT & Rother & Savidis & Shue),]
