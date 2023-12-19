
# Extract the non HDF -----------------------------------------------------

nHDF.Dukhovny <- dataDukhovny |> 
  arrange(desc(row_number())) |> 
  dplyr::select(genes = OFFICIAL_SYMBOL, p.value = SCORE.2) |> 
  mutate(p.value = 1-p.value)

nHDF.Li <- dataLi |> 
  arrange(desc(row_number())) |> 
  dplyr::select(genes = Gene, p.value = p.bh)

nHDF.Shue <- dataShue |> 
  arrange(desc(row_number())) |> 
  dplyr::select(genes, p.value = deseq2.pval)

nHDF.Dukhovny <- nHDF.Dukhovny[1:1000,]
nHDF.Li <- nHDF.Li[1:1000,]
nHDF.Shue <- nHDF.Shue[1:1000,]


genes.nHDF <- rbind(nHDF.Dukhovny, nHDF.Li, nHDF.Shue)
rm(nHDF.Dukhovny, nHDF.Li, nHDF.Shue)

genes.nHDF <- genes.nHDF |> 
  mutate(ensembl.ids = mapIds(org.Hs.eg.db, keys = genes.nHDF$genes, 
                             column = "ENSEMBL", keytype = "ALIAS"), .after = 1)

IDX <- is.na(genes.nHDF$ensembl.ids)
genes.nHDF$ensembl.ids[IDX] <- ensembl_genes$ensembl_gene_id[match(toupper(genes.nHDF$genes[IDX]), ensembl_genes$external_synonym)]
IDX <- is.na(genes.nHDF$ensembl.ids)
genes.nHDF$ensembl.ids[IDX] <- ensembl_genes$ensembl_gene_id[match(genes.nHDF$genes[IDX], ensembl_genes$external_gene_name)]
IDX <- is.na(genes.nHDF$ensembl.ids)
genes.nHDF$ensembl.ids[IDX] <- c(NA, NA, "ENSG00000256892", NA, NA, "ENSG00000256045", 
                                   NA, "ENSG00000221961", "ENSG00000262180")
genes.nHDF.DF <- genes.nHDF |> filter(!is.na(genes.nHDF$ensembl.ids))

genes.nHDF <- genes.nHDF.DF[,2:3] |> 
  group_by(ensembl.ids) |> 
  summarise_all(prod) |> 
  arrange(desc(p.value))

genes.nHDF <- genes.nHDF$ensembl.ids[1:1000]
genes.nHDF <- genes.nHDF[genes.nHDF %in% gene.data$Ensembl_ID]

# Extract the HDF ---------------------------------------------------------

gene.data <- fread("ZIKVData/combinedFinal.csv") |> 
  as.data.frame() |> 
  dplyr::select(all_of(c("Ensembl_ID", finalFeatures)))

genes.Hits <- genes.upset |> 
  dplyr::select(c(1,9)) |> 
  filter(ensembl.ids %in% names)

GS1.of.7 <- genes.Hits$ensembl.ids
GS2.of.7 <- genes.Hits$ensembl.ids[genes.Hits$Intersection >= 2]
GS3.of.7 <- genes.Hits$ensembl.ids[genes.Hits$Intersection >= 3]

# Labeling ----------------------------------------------------------------

gene.data <- gene.data |> mutate(Class = NA)

GS3.data <- gene.data |> 
  mutate(Class = ifelse(Ensembl_ID %in% GS3.of.7, "HDF", 
                        ifelse(Ensembl_ID %in% genes.nHDF, "nHDF", NA))) |> 
  filter(!is.na(Class)) |> 
  mutate(Class = as.factor(Class))

GS2.data <- gene.data |> 
  mutate(Class = ifelse(Ensembl_ID %in% GS2.of.7, "HDF", 
                        ifelse(Ensembl_ID %in% genes.nHDF, "nHDF", NA))) |> 
  filter(!is.na(Class)) |> 
  mutate(Class = as.factor(Class))

GS1.data <- gene.data |> 
  mutate(Class = ifelse(Ensembl_ID %in% GS1.of.7, "HDF", 
                        ifelse(Ensembl_ID %in% genes.nHDF, "nHDF", NA))) |> 
  filter(!is.na(Class)) |> 
  mutate(Class = as.factor(Class))

# Machine Learning --------------------------------------------------------

install.packages('randomForest')
install.packages('caTools')
install.packages("caret", dependencies = TRUE)

library(caTools)
library(caret)


# GS3 ---------------------------------------------------------------------

X3 <- GS3.data[,2:263]
split <- sample.split(X3$Class, SplitRatio = 0.8)
X3_train <- subset(X3, split == TRUE)
X3_test <- subset(X3, split == FALSE)

y3_train <- X3_train[,262]
y3_test <- X3_test[,262]
X3_train <- X3_train[,1:261]
X3_test <- X3_test[,1:261]

ctrl <- trainControl(method = "cv", number = 10, 
                     summaryFunction = twoClassSummary, 
                     classProbs = T, 
                     savePredictions = T)

classifier <- train(x=X3_train, y=y3_train, method = "rf", metric = "ROC", 
                    trControl = ctrl)

y3_pred <- predict(classifier, newdata = X3_test)

cm3 <- confusionMatrix(y3_pred, y3_test)
cm3

# GS2 ---------------------------------------------------------------------

X2 <- GS2.data[,2:263]
split <- sample.split(X2$Class, SplitRatio = 0.8)
X2_train <- subset(X2, split == TRUE)
X2_test <- subset(X2, split == FALSE)

y2_train <- X2_train[,262]
y2_test <- X2_test[,262]
X2_train <- X2_train[,1:261]
X2_test <- X2_test[,1:261]

classifier <- train(x=X2_train, y=y2_train, method = "rf", metric = "ROC", 
                    trControl = ctrl)

y2_pred <- predict(classifier, newdata = X2_test)

cm2 <- confusionMatrix(y2_pred, y2_test)
cm2

# GS1 ---------------------------------------------------------------------

X1 <- GS1.data[,2:263]
split <- sample.split(X1$Class, SplitRatio = 0.8)
X1_train <- subset(X1, split == TRUE)
X1_test <- subset(X1, split == FALSE)

y1_train <- X1_train[,262]
y1_test <- X1_test[,262]
X1_train <- X1_train[,1:261]
X1_test <- X1_test[,1:261]

classifier <- train(x=X1_train, y=y1_train, method = "rf", metric = "ROC", 
                    trControl = ctrl)

y1_pred <- predict(classifier, newdata = X1_test)

cm1 <- confusionMatrix(y1_pred, y1_test)
cm1


