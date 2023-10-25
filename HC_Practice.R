library(tidyverse)
library(caret)
library(cluster)


# Importing the dataset ---------------------------------------------------

data(iris)

# Splitting the dataset into the Training set and Test set ----------------

set.seed(9)
TrainingIndex <- createDataPartition(iris$Species, p=0.8, list = F)
TrainingSet <- iris[TrainingIndex,] # 80% of the data with class labels
TestingSet <- iris[-TrainingIndex,] # mask the species labels for the remaining
TrueSpecies <- TestingSet$Species
TestingSet$Species <- "unknown"
data_unsup <- as.data.frame(rbind(TrainingSet, TestingSet))

X <- data_unsup[1:4]

# Using the dendrogram to find the optimal number of clusters -------------

#Dissimilarity matrix
d <- dist(X, method = "euclidean")

dendrogram <- hclust(d = d, method = "complete")

plot(dendrogram,
     main = paste('Dendrogram'),
     xlab = 'Genes',
     ylab = 'Euclidean distances', 
     cex = 0.6, 
     hang = -1)

# Rotate a set of XY coordinates by an angle:
dend <- rotate(dendrogram, 1:150)

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

iris_species <- levels(iris[,5])
plot(dend,
     main = "Clustered Iris data set (the labels give the true flower species)",
     horiz = TRUE, nodePar = list(cex = .007))
legend("topleft", legend = iris_species, fill = rainbow_hcl(3))

# Fitting Hierarchical Clustering to the dataset --------------------------

hc <- hclust(d = d, method = 'complete')
y_hc <- cutree(hc, 3)

Species <- factor(y_hc,
                  levels = c(1, 2, 3),
                  labels = c('setosa', 'virginica', 'versicolor'))

X2 <- X
X2 <- X2 |> mutate(Species)

sum(as.character(data_unsup$Species) == as.character(X2$Species))
sum(as.character(data_unsup$Species) != as.character(X2$Species))


# Gene selection ----------------------------------------------------------

