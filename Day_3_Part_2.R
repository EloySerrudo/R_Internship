# DAY 3


# Clustering Iris Set ---------------------------------------------------------

library(caret)

data(iris)

set.seed(9)

TrainingIndex <- createDataPartition(iris$Species, p=0.8, list = F)

TrainingSet <- iris[TrainingIndex,] # 80% of the data with class labels
TestingSet <- iris[-TrainingIndex,] # mask the species labels for the remaining

TrueSpecies <- TestingSet$Species

TestingSet$Species <- "unknown"

data_unsup <- as.data.frame(rbind(TrainingSet, TestingSet))


install.packages("installr")

library(dendextend)
library(installr)
library(colorspace)

#Dissimilarity matrix
d <- dist(data_unsup[,1:4], method = "euclidean")

#Hierarchical clustering using Complete Linkage
hc_iris <- hclust(d, method = "complete" )

#Plot the obtained dendrogram
plot(hc_iris, cex = 0.6, hang = -1)

# Why?
iris_species <- levels(iris[,5])

# Convert hclust into a dendrogram and plot
dend <- as.dendrogram(hc_iris)

typeof(dend)
class(dend)
is.dendrogram(dend)

# Rotate a set of XY coordinates by an angle:
dend <- rotate(dend, 1:150)

# Color the branches based on the clusters
dend <- color_branches(dend, k=3)

# Assign for each sample the corresponding iris species as color
labels_colors(dend) <- 
  rainbow_hcl(3)[
    sort_levels_values(as.numeric(iris[,5])[order.dendrogram(dend)])]

labels(dend) <- paste(as.character(iris[,5])[order.dendrogram(dend)], "(", 
                      labels(dend), ")", sep = "")

# Hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.5)

# And plot
par(mar = c(3,3,3,7))
png(filename = "day3_data/Dendo.png", width = 12, height = 12, 
    units = "cm", res = 300)
plot(dend,
     main = "Clustered Iris data set (the labels give the true flower species)",
     horiz = TRUE, nodePar = list(cex = .007))
legend("topleft", legend = iris_species, fill = rainbow_hcl(3))
dev.off()

install.packages('circlize')
library(circlize)

par(mar = rep(0,4))
circlize_dendrogram(dend)

# Gene Clustering ---------------------------------------------------------

install.packages('pvclust')

library(pvclust)
library(gplots)

# We'll perform a gene selection (“feature selection“)
# w/ genes w/ the most relevant information

matrix.org <- read.table("inData/covert.withsymbols.log.Dec05.tab")

# We see which genes have slight variance
variance <- apply(matrix.org, MARGIN = 1, FUN = sd)
variance.rank <- rank(1/variance)
png(filename = "outData/variance2.png", 
    width = 35, height = 20, units = "cm", res = 300)
plot(variance.rank, variance, xlab = "Genes", ylab = "Standard Deviation")
dev.off()

# We see which genes have higher expression
mean.expression <- apply(matrix.org, 1, mean)
int.rank <- rank(1/mean.expression)
png(filename = "day3_data/mean expression.png", 
    width = 35, height = 20, units = "cm", res = 300)
plot(int.rank,mean.expression,xlab = "Genes",ylab = "Expression")
dev.off()

# We select a small group of genes
variance.sel <- which(variance>quantile(variance,0.50))
expression.sel <- which(mean.expression>quantile(mean.expression,0.50))

gene.sel <- intersect(variance.sel, expression.sel)
matrix.sel <- matrix.org[gene.sel, ]

png(filename = "day3_data/gene selection.png", 
    width = 35, height = 20, units = "cm", res = 300)
plot(mean.expression, variance, xlab = "Expression", 
     ylab = "Standard deviation", pch=".")
points(mean.expression[gene.sel], variance[gene.sel], 
       col = "red", pch="+", cex = .7)
dev.off()

# Hierarchical clustering -------------------------------------------------

clustering <- hclust(dist(t(matrix.sel), method = "euclidean"), 
                     method = 'average')

png(filename = "day3_data/Simple cluster 1.png", 
    width = 35, height = 20, units = "cm", res = 300)
plot(clustering)
dev.off()

# Cut by height
cutree(clustering, h = 8)
png(filename = "day3_data/Simple cluster 2.png", 
    width = 35, height = 20, units = "cm", res = 300)
plot(clustering)
abline(h = 8, col = "red")
dev.off()

#Cut tree into 8 groups
png(filename = "day3_data/Simple cluster 3.png", 
    width = 35, height = 20, units = "cm", res = 300)
plot(clustering)
rect.hclust(clustering , k = 8, border = 1:6)
dev.off()

# we can color the branches
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(clustering)
avg_col_dend <- color_branches(avg_dend_obj, h = 8)

png(filename = "day3_data/Simple cluster 4.png", 
    width = 35, height = 20, units = "cm", res = 300)
plot(avg_col_dend)
dev.off()

# Simple heatmap
png(filename = "day3_data/Simple heatmap 1.png", 
    width = 35, height = 20, units = "cm", res = 300)
heatmap(as.matrix(matrix.sel))
dev.off()

# We normalize the rows since we need to absorb the variation between row
png(filename = "day3_data/Simple heatmap 2.png", 
    width = 35, height = 20, units = "cm", res = 300)
heatmap(as.matrix(matrix.sel), scale="row")
dev.off()

# Customizing the heatmap
png(filename = "day3_data/Simple heatmap 3.png", 
    width = 35, height = 20, units = "cm", res = 300)
heatmap(as.matrix(matrix.sel), scale="row" , col = cm.colors(256), 
        xlab="Sample", ylab="Genes", main="heatmap", labCol = F)
dev.off()

library(gplots) #for redgreen(75)

png(filename = "day3_data/variance heatmap.png", width = 35, height = 20, 
    units = "cm", res = 300)
heatmap.2(as.matrix(matrix.sel), Rowv=T, Colv=T, key=F, 
          scale="row", col = redgreen(75), symkey=T,
          density.info = "histogram", trace="none")
dev.off()

# k-means clustering ------------------------------------------------------

kmclustering <- hclust(dist(t(matrix.sel), method = "euclidean"), method = 'average')

number_of_clusters <- 2
set.seed(5)
km <- kmeans(t(matrix.sel), centers=number_of_clusters)

# The number of samples for each cluster can be seen
table(km$cluster)

my.col <- palette()[2:3]
View(my.col)
col <- my.col[as.vector(km$cluster)]
samp.col <- col[order(km$cluster)]

heatmap.2(as.matrix(matrix.sel), Rowv=FALSE, Colv=as.dendrogram(clustering), scale="row",
          dendrogram="col", ColSideColors=samp.col, col = redgreen(75), density.info="histogram",
          trace="none", labRow=FALSE)

col <- my.col[as.vector(km$cluster)]
samp.col <- col[order(km$cluster)]
