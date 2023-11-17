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

z.gseaSavidis.mod <- adaptGSEA(gseaSavidis, "")
z.gseaLi.mod <- adaptGSEA(gseaLi, "")
z.gseaDukhovny.mod <- adaptGSEA(gseaDukhovny, "")
z.gseaWangGSC.mod <- adaptGSEA(gseaWangGSC, "")
z.gseaWang293FT.mod <- adaptGSEA(gseaWang293FT, "")
z.gseaRother.mod <- adaptGSEA(gseaRother, "")
z.gseaShue.mod <- adaptGSEA(gseaShue, "")

z.GSL.1 <- rbind(
  z.gseaSavidis.mod[,c("ID", "p.adjust", "geneID")],
  z.gseaLi.mod[,c("ID", "p.adjust", "geneID")],
  z.gseaDukhovny.mod[,c("ID", "p.adjust", "geneID")],
  z.gseaWangGSC.mod[,c("ID", "p.adjust", "geneID")],
  z.gseaWang293FT.mod[,c("ID", "p.adjust", "geneID")],
  z.gseaRother.mod[,c("ID", "p.adjust", "geneID")],
  z.gseaShue.mod[,c("ID", "p.adjust", "geneID")]
)

z.GSL.2 <- z.GSL.1[match(unique(z.GSL.1$ID), z.GSL.1$ID), ]

# REMOVE REDUNDANCY -------------------------------------------------------

z.reduced.GSL.1 <- removeRedundant(inputData = z.GSL.1, 
                                   jacCutoff = 0.3, 
                                   model = "weighted",
                                   gprofilerOutput = F, 
                                   parallel = F)
row.names(z.reduced.GSL.1) <- NULL

z.reduced.GSL.2 <- removeRedundant(inputData = z.GSL.2, 
                                   jacCutoff = 0.3, 
                                   model = "weighted",
                                   gprofilerOutput = F, 
                                   parallel = F)
row.names(z.reduced.GSL.2) <- NULL

# Al final si tomamos los únicos de z.reduced.GSL.1 se parece a
# z.reduced.GSL.2 porque toma el primero más largo

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
# GO.Terms <- read.csv("ZIKVData/GO_Terms.csv")


# Build a new GMT File ----------------------------------------------------

GO.Terms.Reduced <- GO.Terms[GO.Terms$ID %in% z.reduced.GSL.2$ID, ]
rownames(GO.Terms.Reduced) <- NULL

new_gmt <- apply(GO.Terms.Reduced[, -2], 1, function(row) {
  paste(c(row[1], row[2], unlist(strsplit(as.character(row[3]), ","))), collapse = "\t")
})
writeLines(new_gmt, con = "ZIKVData/new_gmt.name.gmt")

gostres <- gost(query = sort(Hits.Shue), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)
View(gostres$result)

# NEW DATA ----------------------------------------------------------------

new.gseaSavidis <- gost(query = sort(Hits.Savidis), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)
new.gseaLi <- gost(query = sort(Hits.Li), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)
new.gseaDukhovny <- gost(query = sort(Hits.Dukhovny), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)
new.gseaWangGSC <- gost(query = sort(Hits.WangGSC), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)
new.gseaWang293FT <- gost(query = sort(Hits.Wang293FT), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)
new.gseaRother <- gost(query = sort(Hits.Rother), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)
new.gseaShue <- gost(query = sort(Hits.Shue), organism = "gp__EhA0_PFup_VAE", sources = c("GO"), evcodes = T)

z.new.gseaSavidis.mod <- adaptGSEA(new.gseaSavidis, "")
z.new.gseaLi.mod <- adaptGSEA(new.gseaLi, "")
z.new.gseaDukhovny.mod <- adaptGSEA(new.gseaDukhovny, "")
z.new.gseaWangGSC.mod <- adaptGSEA(new.gseaWangGSC, "")
z.new.gseaWang293FT.mod <- adaptGSEA(new.gseaWang293FT, "")
z.new.gseaRother.mod <- adaptGSEA(new.gseaRother, "")
z.new.gseaShue.mod <- adaptGSEA(new.gseaShue, "")

z.new.gseaSavidis.mod$Category <- GO.Terms.Reduced$Category[match(z.new.gseaSavidis.mod$ID, GO.Terms.Reduced$ID)]
z.new.gseaLi.mod$Category <- GO.Terms.Reduced$Category[match(z.new.gseaLi.mod$ID, GO.Terms.Reduced$ID)]
z.new.gseaDukhovny.mod$Category <- GO.Terms.Reduced$Category[match(z.new.gseaDukhovny.mod$ID, GO.Terms.Reduced$ID)]
z.new.gseaWangGSC.mod$Category <- GO.Terms.Reduced$Category[match(z.new.gseaWangGSC.mod$ID, GO.Terms.Reduced$ID)]
z.new.gseaWang293FT.mod$Category <- GO.Terms.Reduced$Category[match(z.new.gseaWang293FT.mod$ID, GO.Terms.Reduced$ID)]
z.new.gseaRother.mod$Category <- GO.Terms.Reduced$Category[match(z.new.gseaRother.mod$ID, GO.Terms.Reduced$ID)]
z.new.gseaShue.mod$Category <- GO.Terms.Reduced$Category[match(z.new.gseaShue.mod$ID, GO.Terms.Reduced$ID)]

z.genes.Savidis <- paste(z.new.gseaSavidis.mod$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()
z.genes.Li <- paste(z.new.gseaLi.mod$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()
z.genes.Dukhovny <- paste(z.new.gseaDukhovny.mod$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()
z.genes.WangGSC <- paste(z.new.gseaWangGSC.mod$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()
z.genes.Wang293FT <- paste(z.new.gseaWang293FT.mod$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()
z.genes.Rother <- paste(z.new.gseaRother.mod$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()
z.genes.Shue <- paste(z.new.gseaShue.mod$geneID, collapse = ",") |> str_split(',') |> 
  unlist() |> unique() |> sort()

z.genes.2 <- c(z.genes.Li, z.genes.Rother, z.genes.Dukhovny, z.genes.Shue, z.genes.WangGSC, 
               z.genes.Wang293FT)
z.genes.2 <- z.genes.2 |> unique() |> sort()

z.genes.2.DF <- data.frame(z.genes.2, 
                           z.genes.2 %in% z.genes.WangGSC, 
                           z.genes.2 %in% z.genes.Wang293FT, 
                           z.genes.2 %in% z.genes.Shue, 
                           z.genes.2 %in% z.genes.Savidis, 
                           z.genes.2 %in% z.genes.Rother,
                           z.genes.2 %in% z.genes.Li, 
                           z.genes.2 %in% z.genes.Dukhovny,
                           row.names = NULL)
colnames(z.genes.2.DF)[-1] <- c("Wang.GSC", "Wang.293FT", "Shue", "Savidis", 
                              "Rother", "Li", "Dukhovny")

# Plotting the Upset plot -------------------------------------------------

upset(
  data = z.genes.2.DF,
  intersect = papers,
  name = 'Papers', 
  width_ratio=0.2,
  sort_intersections_by = 'degree',
  sort_intersections = 'ascending',
  sort_sets = FALSE
) +
  ggtitle("New Genes from 6 studies")


z.genes.DF.int <- z.genes.2.DF |> 
  rowwise() |>
  mutate(Intersection = sum(c_across(2:8))) |> 
  filter(Intersection > 1)

z.genes.DF.int |> 
  filter(Li & Wang.293FT & Rother & Savidis & Shue & Wang.GSC) |> 
  select(z.genes.2) |> 
  unlist(use.names = FALSE)



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
