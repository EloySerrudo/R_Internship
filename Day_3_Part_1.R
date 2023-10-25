# DAY 3
install.packages("caret", dependencies = TRUE)
install.packages("grid", dependencies = TRUE)
install.packages("gridExtra", dependencies = TRUE)
install.packages("factoextra", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("rpart.plot", dependencies = TRUE)
install.packages("ggdendro", dependencies = TRUE)


library(tidyverse)
library(caret)
library(factoextra) # because of ‘fviz_eig()’ function
library(grid)
library(gridExtra)


# Importing the dataset ---------------------------------------------------

data(iris)
head(iris)
tail(iris)
dim(iris)
summary(iris)
str(iris)

sum(is.na(iris)) # How many NAs there are

#newdata <- na.omit(mydata) # If we want to create new data w/o NAs


# Exploratory analysis ----------------------------------------------------

set.seed(7)
iris.shuffled <- iris[sample(1:nrow(iris), size = nrow(iris), replace = F),]

head(iris.shuffled)

iris.shuffled$ID <- as.character(rownames(iris.shuffled)) # Added a new col. w/ char IDs

scales = list(x=list(relation="free"), y=list(relation="free"))
png(filename = "day3_data/featurePlot_Density.png", 
    width = 15, height = 15, units = "cm",res = 300)
featurePlot(x = iris.shuffled[, 1:4], y = iris.shuffled$Species, plot = "pairs",
            scales = scales, auto.key = list(columns = 3))
dev.off()

png(filename = "day3_data/featurePlot_Density.png", 
    width = 15, height = 15, units = "cm", res = 300)
featurePlot(x = iris.shuffled[, 1:4], y = iris.shuffled$Species, plot = "box", 
            scales = scales, par.strip.text=list(cex=0.6)) #!!!!!!!!!!!!!!
dev.off()

# Unsupervised classification ---------------------------------------------

set.seed(9)

# Creating the training set (70%)
TrainingIndex1 <- createDataPartition(iris.shuffled$Species, p=0.7, list = F)
TrainingSet1 <- iris.shuffled[TrainingIndex1,]
dim(TrainingSet1)

# Creating the training set (30%)
TestingSet1 <- iris.shuffled[-TrainingIndex1,]
TestingSet1$Species <- "unknown"
dim(TestingSet1)

# Mine
library(caTools)
split <- sample.split(iris.shuffled$Species, SplitRatio = 0.7)
TrainingSet2 <- subset(iris.shuffled, split == TRUE)
dim(TrainingSet2)

TestingSet2 <- subset(iris.shuffled, split == FALSE)
TestingSet2$Species <- "unknown"
dim(TestingSet2)

#


# Unsupervised Data
data_unsup <- as.data.frame(rbind(TrainingSet, TestingSet))

# PCA ---------------------------------------------------------------------

pca <- prcomp(data_unsup[,1:4], center = T, scale. = T)

typeof(pca)
class(pca)

png(filename = "day3_data/Scree.Plot.png", width = 12, height = 12, 
    units = "cm", res = 300)
fviz_eig(pca, addlabels = T)
 dev.off()

df_out <- as.data.frame(pca$x)
head(df_out)

df_out$Species <- as.character(data_unsup[,5])

p1 <- 
  ggplot(df_out, aes(x=PC1, y=PC2, color=Species, label=rownames(df_out))) +
  geom_point() + 
  geom_text(aes(label=rownames(df_out)), hjust=0, vjust=0) +
  theme_bw()

p2 <- 
  ggplot(df_out, aes(x=PC1, y=PC3, color=Species, label=rownames(df_out))) +
  geom_point() +
  geom_text(aes(label=rownames(df_out)), hjust=0, vjust=0) +
  theme_bw()

p3 <- 
  ggplot(df_out, aes(x=PC2, y=PC3, color=Species, label=rownames(df_out))) +
  geom_point() +
  geom_text(aes(label=rownames(df_out)), hjust=0, vjust=0) +
  theme_bw()

pFin <- grid.arrange(p1,p2,p3, ncol=2)

ggsave(pFin, filename = "day3_data/PCA.png", device = "png", 
       dpi = 600, width = 30, height = 30, units = "cm")

# Supervised Classification -----------------------------------------------

# Setosa flowers are easy identifiable, so we gonna remove it
iris <- iris[-c(which(iris$Species == "setosa")),]
iris$Species <- factor(iris$Species)

set.seed(9)
# Then we separate the data again (70%):
training_indices <- createDataPartition(iris$Species, p=0.7, list = F)

irisTrain <- iris[training_indices,] # Train set
irisTest <- iris[-training_indices,] # Test set

# Training the model w/ 10 fold cross validation and 3 repeats
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                        summaryFunction = multiClassSummary, classProbs = T, 
                        savePredictions = T)
metric <- "Accuracy"

# K Nearest Neighbors ----------------------------------------------------

set.seed(7)

fit.knn <- train(x=irisTrain[,1:4], y=irisTrain[,5], method = "knn", 
                 metric = metric, trControl = control, 
                 preProcess = c("center", "scale"))

# Decision Tree -----------------------------------------------------------

set.seed(7)

fit.cart <- train(x=irisTrain[,1:4], y=irisTrain[,5], method = "rpart", 
                  metric = metric, trControl = control, 
                  preProcess = c("center","scale"))

library(rpart.plot)

rpart.plot(fit.cart$finalModel, type = 5)

# Random Forest -----------------------------------------------------------

set.seed(7)

fit.rf <- train(x=irisTrain[,1:4], y=irisTrain[,5], method = "rf", 
                metric = metric, trControl = control,
                preProcess = c("center", "scale"))

# Support Vector Machines -------------------------------------------------

set.seed(7)

fit.svm_Lin <- train(x=irisTrain[,1:4], y=irisTrain[,5], method = "svmLinear", 
                     metric = metric, trControl = control, 
                     preProcess = c("center", "scale"))

# Artificial neuronal network ---------------------------------------------

set.seed(7)

fit.nnet <- train(x=irisTrain[,1:4], y=irisTrain[,5], method = "nnet", 
                  metric = metric, trControl = control, 
                  preProcess = c("center", "scale"))

# LDA: Linear Discriminant Analysis ---------------------------------------

set.seed(7)

fit.lda <- train(x=irisTrain[,1:4], y=irisTrain[,5], method = "lda", 
                 metric = metric, trControl = control, 
                 preProcess = c("center", "scale"))

# Summarize all ML model performance --------------------------------------

results <- resamples(list(lda=fit.lda, knn=fit.knn, svm_Lin=fit.svm_Lin, 
                          rf=fit.rf, cart=fit.cart, nnet=fit.nnet))
summary(results) # svm_Rad=fit.svm_Rad is missing!!!!!!!!!!!!!!!!!!!!!

scales <- list(x=list(relation="free"), y=list(relation="free"))
png(filename = "day3_data/ML_performance.png", width = 20, height = 20, 
    units = "cm", res = 300)
dotplot(results, scales=scales, par.strip.text=list(cex=0.76), 
        par.settings = list(par.xlab.text = list(cex = 0)))
dev.off()

# Testing the Models ------------------------------------------------------

predictions <- predict(fit.cart, irisTest[,1:4]) # We predict the outcomes
confusionMatrix(predictions, irisTest$Species) # And we compare them

gbmImp <- varImp(fit.svm_Lin, scale = T)
png(filename = "day3_data/Iris_ML_FeatImp.png", width = 10, height = 10, 
    units = "cm", res = 300)
plot(gbmImp, main="Iris Variable Importance (SVM)", top=4)
dev.off()