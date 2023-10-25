# DAY 4
library(pvclust)
library(gplots)

matrix.org <- read.table("inData/covert.withsymbols.log.Dec05.tab")

typeof(matrix.org)
class(matrix.org )
dim(matrix.org)

# We'll rank the genes according to their standard deviation and mean

# See low variance of expression genes
variance <- apply(matrix.org, MARGIN = 1, FUN = sd)
variance.rank <- rank(1/variance)
plot(variance.rank,variance,xlab = "Genes",ylab = "Standard Deviation")

# See high expression genes
mean.expression <- apply(matrix.org, 1, mean)
int.rank <- rank(1/mean.expression)
plot(int.rank,mean.expression,xlab = "Genes",ylab = "Expression")

# Selected genes ----------------------------------------------------------

# Remove low variance and expression genes
variance.sel <- which(variance>quantile(variance,0.50))
expression.sel <- which(mean.expression>quantile(mean.expression,0.50))

length(variance.sel)
length(expression.sel)

# Here we select the genes that accomplish both criteria
gene.sel <- intersect(variance.sel, expression.sel)
matrix.sel <- matrix.org[gene.sel, ]

plot(mean.expression, variance, xlab = "Expression", 
     ylab = "Standard deviation", pch=".")
points(mean.expression[gene.sel], variance[gene.sel], 
       col = "red", pch="+", cex = .7)

# Exercise ----------------------------------------------------------------

percentage <- 0.80
variance.ex <- which(variance>quantile(variance, percentage))
expression.ex <- which(mean.expression>quantile(mean.expression, percentage))
length(variance.ex)
length(expression.ex)

gen.ex <- intersect(variance.ex, expression.ex)
length(gene.ex)

# Hierarchical clustering -------------------------------------------------

# We choose this distance measures and methods
method <- "average"
distance <- "euclidean"

clustering <- hclust(dist(t(matrix.sel), method = distance), method = method)
plot(clustering)

# Splitting the number of clusters
groups <- cutree(clustering,6)
samples <- colnames(matrix.sel)
print(samples[groups==1])
print(samples[groups==2])
print(samples[groups==3])
print(samples[groups==4])
print(samples[groups==5])
print(samples[groups==6])

heatmap.2(as.matrix(matrix.sel), Rowv=TRUE, Colv=as.dendrogram(clustering), 
          key=TRUE, scale="row", col = redgreen(75), symkey=FALSE, 
          density.info = "histogram", trace="none", labRow=FALSE)
# labRow shows the name of the genes

# Cluster stability -------------------------------------------------------

set.seed(5)
pv <- pvclust(matrix.sel, method.hclust=method, method.dist=distance, nboot=10)
# nboot means number of bootstrapings

p.value <- 0.05

plot(pv)
pvrect(pv, alpha = 1-p.value)

# What did we just do here?
# A p.value of 5% means that we want to be sure that the clusters formed at the
# end to be in at least 95% (1-p) of the bootstraped clusters (that means that
# the clusters are Significant). That's the reason for the alpha value. The AU
# values give us this degree of certainly.

outliers <- c(grep("ec_aer_fnr_nO_a", colnames(matrix.sel)),
              grep("ec_aer_appY_O_c", colnames(matrix.sel)),
              grep("ec_aer_oxyR_nO_c", colnames(matrix.sel)),
              grep("ec_aer_appY_nO_c", colnames(matrix.sel)))

# grep() gives us the indexes of the outliers

new.matrix <- matrix.sel[,-outliers]

dim(matrix.sel)
dim(new.matrix)

pv0 <- pvclust(new.matrix, method.hclust=method, method.dist=distance, nboot=10)
plot(pv0)
pvrect(pv0, alpha = 1-p.value)

# Now we're gonna check the stability of the gene clusters

# t() This function transpose de Matrix or Dataframe
pv1 <- pvclust(t(matrix.sel), method.hclust=method, method.dist=distance, nboot=10)
# We' won't plot pv1 because there're too many things on the x-axis. So, we'll
# just display a numerical value to see the clusters formed.

pv1.pp <- pvpick(pv1, alpha = 1-p.value)

k <- length(pv1.pp$clusters)
print(k) # Number of stable clusters. Ther're just two clusters, because the low
# value of boostraps that we did. If we increase this number to 100 for instance
# this value will be better, but it will take a long.

# Interpretation (functional analysis) ------------------------------------

clustering.genes = hclust(dist(matrix.sel, method = distance), method = method)
plot(clustering.genes)
groups = cutree(clustering.genes, k=2) # k=2 porque sólo tenemos 2 clústers
# Note: group numbers are different in 'pvclust()'
gene = rownames(matrix.sel)

cluster.sel <- 2
genes_cluster <- gene[!groups==cluster.sel] # vector con el nombre de los genes
print(length(genes_cluster))

# K-Means Clustering ------------------------------------------------------

number_of_clusters <- 2
km <- kmeans(t(matrix.sel), centers=number_of_clusters)
table(km$cluster)

my.col = palette()[2:3]
col = my.col[as.vector(km$cluster)]
samp.col = col[order(km$cluster)]

heatmap.2(as.matrix(matrix.sel), Rowv=FALSE, Colv=as.dendrogram(clustering), 
          scale="row", dendrogram="col", ColSideColors=samp.col, 
          col=redgreen(75), density.info="histogram", trace="none", labRow=FALSE)

number_of_clusters <- 8
km <- kmeans(t(new.matrix), centers=number_of_clusters)
table(km$cluster)

my.col_25 = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", 
             "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", 
             "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", 
             "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", 
             "darkorange4", "brown")

ncol(as.matrix(new.matrix))

library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

my.col_39 <- col_vector[1:39]

col = my.col_39[as.vector(km$cluster)]
samp.col <- col[order(km$cluster)]

clustering.new <- hclust(dist(t(new.matrix), method = distance), method = method)
heatmap.2(as.matrix(new.matrix), Rowv=FALSE, Colv=as.dendrogram(clustering.new), 
          scale="row", dendrogram="col", ColSideColors=samp.col, 
          col=redgreen(75), density.info="histogram", trace="none", labRow=FALSE)

