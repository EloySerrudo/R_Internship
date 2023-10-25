# German Credit Excercise

data(GermanCredit, package = "caret")

dim(GermanCredit)
str(GermanCredit)
head(GermanCredit)
summary(GermanCredit)

sum(is.na(GermanCredit)) # How many NAs there are

colnames(GermanCredit)

DF <- GermanCredit
DF <- DF |> relocate(Class, .after = Job.Management.SelfEmp.HighlyQualified)

NZV_indexes <- nearZeroVar(DF)
DF <- DF[,-NZV_indexes]

set.seed(9)

split <- sample.split(DF$Class, SplitRatio = 0.7)
train_set <- subset(DF, split == T)
test_set <- subset(DF, split == F)
rm(split)

control <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 10,
                        summaryFunction = multiClassSummary, 
                        classProbs = T, 
                        savePredictions = T)
metric <- "Accuracy"

set.seed(7)

# K Nearest Neighbors ----------------------------------------------------

fit.knn <- train(x = train_set[,1:48], y = train_set[,49], method = "knn", 
                 metric = metric, trControl = control, 
                 preProcess = c("center", "scale"))

# Support Vector Machines -------------------------------------------------

fit.svm_Lin <- train(x = train_set[,1:48], y = train_set[,49], method = "svmLinear", 
                     metric = metric, trControl = control, 
                     preProcess = c("center", "scale"))

# Random Forest -----------------------------------------------------------

fit.rf <- train(x = train_set[,1:48], y = train_set[,49], method = "rf", 
                metric = metric, trControl = control,
                preProcess = c("center", "scale"))

# Decision Tree -----------------------------------------------------------

fit.cart <- train(x = train_set[,1:48], y = train_set[,49], method = "rpart", 
                  metric = metric, trControl = control, 
                  preProcess = c("center","scale"))

library(rpart.plot)

rpart.plot(fit.cart$finalModel, type = 5)

# LDA: Linear Discriminant Analysis ---------------------------------------

fit.lda <- train(x = train_set[,1:48], y = train_set[,49], method = "lda", 
                 metric = metric, trControl = control, 
                 preProcess = c("center", "scale"))

# Artificial neuronal network ---------------------------------------------

fit.nnet <- train(x = train_set[,1:48], y = train_set[,49], method = "nnet", 
                  metric = metric, trControl = control, 
                  preProcess = c("center", "scale"))

# Testing the Models ------------------------------------------------------

predictions <- predict(fit.nnet, test_set[,1:48]) # We predict the outcomes
confusionMatrix(predictions, test_set$Class) # And we compare them

gbmImp <- varImp(fit.svm_Lin, scale = T)
plot(gbmImp, main="Iris Variable Importance (SVM)", top=10)

# Summarize all ML model performance --------------------------------------

results <- resamples(list(lda=fit.lda, knn=fit.knn, svm_Lin=fit.svm_Lin, 
                          rf=fit.rf, cart=fit.cart, nnet=fit.nnet))
summary(results)

scales <- list(x=list(relation="free"), y=list(relation="free"))

dotplot(results, scales=scales, par.strip.text=list(cex=0.76), 
        par.settings = list(par.xlab.text = list(cex = 0)))
