# DAY 5
# Related with Day 2 Part 1
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer) # Here there are palettes of colors
library(VennDiagram)
library(ComplexUpset)


rawCounts <- read.table("inData/sarsCov2_rawCounts_smybols.txt")
colnames(rawCounts)
str(rawCounts)
head(rawCounts)

# splits the string
colNamesList <- strsplit(colnames(rawCounts), 'X')
# keeps only the second part and use it as column names
colnames(rawCounts) <- sapply(colNamesList, '[[', 2)

# The data originally are integers, we change them to floats
rawCounts[] <- lapply(rawCounts, as.numeric) # Creo que esto no es necesario

sampleTable <- read.table("inData/sarsCov2_sampleTable_smybols.txt")
sampleTable <- mutate(sampleTable, # Here we change the type of the column
                      timepoint = factor(sampleTable$timepoint,
                                         c('0h','3h','6h','12h','24h')))

ddsMat <- DESeqDataSetFromMatrix(countData = rawCounts, 
                                 colData = sampleTable, 
                                 design = ~timepoint)

# PCA Plot ----------------------------------------------------------------

rld <- rlog(ddsMat, blind = F) # Data is transformed to logarithmic scale

data <- plotPCA(rld, intgroup=c("virus", "timepoint"), returnData=T)
percentVar <- round(100*attr(data,"percentVar"))

ggplot(
  data = data, 
  mapping = aes(PC1, PC2, color=timepoint, shape=virus)
) + 
  geom_point(size = 3) +
  labs(
    title = "SARS-Cov2 Data",
    subtitle = "PCA Data",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
  ) + 
  scale_color_manual(
    breaks = c('0h','3h','6h','12h','24h'),
    values=c("#B3E5FC", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
  )

# The Heatmap -------------------------------------------------------------

ddsDEG <- DESeq(ddsMat)
summary(results(ddsDEG, alpha = 0.05, lfcThreshold = 0.5))

select <- order(rowMeans(counts(ddsDEG, normalized=TRUE)), 
                decreasing=TRUE)[1:50]
# We select the 50 genes with the highest values

DF <- as.data.frame(colData(ddsDEG)[,c("virus","timepoint")])
annColors <- list(timepoint = c('0h' = "#B3E5FC",'3h' = "#66C2A5", 
                                '6h' = "#FC8D62", '12h' = "#8DA0CB",
                                '24h'="#E78AC3"),
                  virus = c(mock= "#ded4d4", infected="#b9b0b0"))
pheatmap(assay(rld)[select,], cluster_rows=T,
         show_rownames=FALSE, cluster_cols=T, 
         annotation_col=DF, annotation_colors=annColors)

# Bar plots ---------------------------------------------------------------

res <- read.table("inData/sarsCov2_DEG_smybols.txt", header = T)
res <- res |> 
  filter(log2FoldChange < (-1) | log2FoldChange > 1)

res$direction <- NA
res$direction[res$log2FoldChange < 0] <- "down"
res$direction[res$log2FoldChange > 0] <- "up"

res <- res |> 
  mutate(condition = factor(condition, levels = c("timepoint_3h_vs_0h",
                                                  "timepoint_6h_vs_0h",
                                                  "timepoint_12h_vs_0h",
                                                  "timepoint_24h_vs_0h")))

res$TotalNumDEGs <- NA
for (i in 1:length(unique(res$condition))) {
  res$TotalNumDEGs[res$condition == unique(res$condition)[i]] <- 
    length(res$Gene[res$condition == unique(res$condition)[i]])
}

res$NumDEGs <- NA
for (i in 1:length(unique(res$condition))) {
  res$NumDEGs[
    res$condition == unique(res$condition)[i] & res$direction == "up"] <-
    length(res$Gene[
      res$condition == unique(res$condition)[i] & res$direction == "up"])
  
  res$NumDEGs[
    res$condition == unique(res$condition)[i] & res$direction == "down"] <-
    length(res$Gene[
      res$condition == unique(res$condition)[i] & res$direction == "down"])
}

ggplot(
  data = res, 
  mapping = aes(x=condition, fill=condition)
) +
  geom_bar() +
  labs(
    title = "Total Number of DEGs",
    x ='Condition', y = 'Counts'
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="none")

ggplot(
  data = res, 
  mapping = aes(fill=direction, y=NumDEGs, x=condition)
) +
  geom_bar(position='dodge', stat='identity') +
  theme_minimal() +
  labs(
    title = "Number of DEGs (Up & Down)",
    x ='Condition', y = 'Counts'
  ) +
  scale_fill_manual('Direction', values=c('coral2', 'steelblue'))

# Pie charts --------------------------------------------------------------

res024 <- res |> 
  filter(condition == "timepoint_24h_vs_0h") |> 
  select(condition, direction, NumDEGs) |> 
  unique()

ggplot(
  data = res024, 
  mapping = aes(x="", y=NumDEGs, fill=direction)
  ) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_brewer(palette = "Dark2")

res012 <- res |> 
  filter(condition == "timepoint_12h_vs_0h") |> 
  select(condition, direction, NumDEGs) |> 
  unique()

ggplot(
  data = res012, 
  mapping = aes(x="", y=NumDEGs, fill=direction)
) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_brewer(palette = "Dark2")

# Box-(Whisker)-Plot ------------------------------------------------------

resUp <- res |> filter(direction == "up")

ggplot(
  data = resUp, 
  mapping = aes(x=condition, y=log2FoldChange, fill=condition)
) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set2") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("A boxplot with jitter") +
  xlab("") # sin nombre en el eje X


ggplot(
  data = resUp,
  mapping = aes(x=condition, y=log2FoldChange, fill=condition)
) +
  geom_violin() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("A violin plot") +
  xlab("")

# Venn Diagrams -----------------------------------------------------------

brewer.pal(n=4,"Set2")

grid.newpage()
venn_object <- venn.diagram(
  x = list(
    unname(unlist(select(filter(res, condition == "timepoint_3h_vs_0h"), Gene))),
    unname(unlist(select(filter(res, condition == "timepoint_6h_vs_0h"), Gene))),
    unname(unlist(select(filter(res, condition == "timepoint_12h_vs_0h"), Gene))),
    unname(unlist(select(filter(res, condition == "timepoint_24h_vs_0h"), Gene)))),
  category.names = c("3h" , "6h" , "12h", "24h"),
  filename = 'outData/timepoint_venn.png', # output file name
  #output = TRUE ,
  #imagetype="png" , # creates output as PNG file
  height = 480 , # height of the picture
  width = 480 , # width of the picture
  resolution = 300, # resolution
  compression = "lzw",
  lwd = 1,
  col=brewer.pal(n=4,"Set2"),
  fill = c(
    alpha("#66C2A5",0.3), 
    alpha('#FC8D62',0.3), 
    alpha('#8DA0CB',0.3), 
    alpha('#E78AC3',0.3)),
  cex = 0.5, # label size
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = brewer.pal(n=4,"Set2"),
)
grid.draw(venn_object)

typeof(venn_object)
class(venn_object)

# Upset plot --------------------------------------------------------------

upsetList <- list(
  timepoint_3h = unname(unlist(select(filter(res, condition =="timepoint_3h_vs_0h"), Gene))),
  timepoint_6h = unname(unlist(select(filter(res, condition == "timepoint_6h_vs_0h"), Gene))),
  timepoint_12h = unname(unlist(select(filter(res, condition == "timepoint_12h_vs_0h"), Gene))),
  timepoint_24h = unname(unlist(select(filter(res, condition == "timepoint_24h_vs_0h"), Gene))))

sets <- c()
for (i in 1:length(upsetList)) {
  sets <- unique(c(sets, upsetList[[i]]))
}

upsetData <- c()
for (i in 1:length(upsetList)) {
  upsetData <- rbind(upsetData, sets %in% upsetList[[i]])
}

upsetData <- as.data.frame(t(as.matrix(upsetData)))
colnames(upsetData) <- names(upsetList) # sets names to the columns

upset(data = upsetData, intersect = colnames(upsetData), name = "", min_size=0)

upset(
  data = upsetData, 
  intersect = colnames(upsetData), 
  name = "", 
  min_size=0,
  base_annotations = list('Intersection size'=intersection_size(counts = T)),
  matrix = intersection_matrix(
    geom = geom_point(shape='circle filled', size=3),
    segment = geom_segment(size=.8)
  ) +
    scale_color_manual(values = c('3h'='#66C2A5','6h'='#FC8D62','12h'='#8DA0CB','24h'='#E78AC3'), 
                       guide=guide_legend(override.aes = list(shape='circle'))),
    queries = list(
      upset_query(set="timepoint_3h",fill="#66C2A5"),
      upset_query(set="timepoint_6h",fill="#FC8D62"),
      upset_query(set="timepoint_12h",fill="#8DA0CB"),
      upset_query(set="timepoint_24h",fill="#E78AC3")
    )
)

ggsave("outData/upsetPlot_upregulatedGenes.png", 
       width = 15, height = 10, units = "cm", dpi = 500)
