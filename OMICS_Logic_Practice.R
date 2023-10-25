library(tidyverse)

path <- 'inData/Final_cell_lines_RNA-expression_FPKM_values_1000genes_with_NA.txt'
DF <- read.table(path, header=TRUE) # row.names=1) if we want name de rows with the first column

rm(path)

dim(DF)
colnames(DF)
str(DF)

colSums(is.na(DF))
colSums(DF[,-1]==0)

DF_new <- DF |> na.exclude()
rownames(DF_new) <- DF_new[,1]
DF_new <- DF_new[-1]
dim(DF_new)

# Descriptive statistics, i.e., compute the summary statistics of data
summ <- summary(DF_new)
# print summary statistics
print(summ)
# Write summary statistics of data  into a file and export it
write.table(summ, file="stat_sum-FPKM-data.txt", col.names=TRUE, sep="\t")

# Draw boxplot for all samples (for FPKM values)
boxplot(DF_new)
# Add the main title, axis titles, and color 
boxplot(DF_new, 
        main="Boxplot for FPKM data", 
        xlab="", ylab="Gene expression (FPKM)", 
        col="red", las=1, cex.axis = 0.65)

log_df <- log(DF_new+1)
#Check the dimension of data
dim(log_df)
# let’s check the summary statistics of data after log transformation
summ_log <- summary(log_df)
# Print summary statistics of data 
summ_log
# Write summary statistics to a file and export it
write.table(summ_log, file="stat_sum-log-data.txt", col.names=TRUE, sep="\t")
# Draw boxplot for all samples (for FPKM values)
boxplot(log_df, 
        main="Boxplot for log-transformed Data", 
        xlab="", ylab="Gene expression (log[FPKM])", 
        col="red", las=2, cex.axis = 0.7)

# Extract specific samples from data with their name
Non_ML1 <- as.data.frame(log_df[grep('^Non.malignant', names(log_df))])
Claudin1 <- as.data.frame(log_df[grep('^Claudin', names(log_df))])
#Compute the mean for each gene within each group
Non_ML<- rowMeans(Non_ML1) # Media de cada fila
Claudin <- rowMeans(Claudin1)

#Bind groups together into 2 columns
group <- cbind(Non_ML, Claudin) # Retorna una matriz
#Provide the column names to both groups
colnames(group) <- c("Non-malignant", "Claudin-low")

#Draw boxplot for groups of samples
boxplot(group, 
        main="Boxplot for groups", 
        xlab="Groups", ylab="Gene expression(log[FPKM])", 
        col=c('cyan', 'pink'))
