install.packages("gprofiler2")

library(boot)
library(BiocManager)
library(dplyr)
library(gprofiler2)

BiocManager::install()
BiocManager::install("Biobase")

library(Biobase)


# 4. Analyzing the Netea data set -----------------------------------------

fr <- read.table("part2_data/Neteadata.txt")
typeof(fr)
class(fr)
str(fr)

M <- as.matrix(fr[-1,])
h <- dim(M)

M[1:5, 1:5]
typeof(M)
class(M)

genename <- rownames(M)
sample <- colnames(M)

M <- as.numeric(M) # NOTE: This turns M to a 1D vector

dim(M) <- h
num.samples <- h[2]
num.genes <- h[1]
rownames(M) <- genename
colnames(M) <- sample

sampleclass <- as.vector(fr[1,]) #Here we put the names that we left in line 21
sampleclass[1:5]

# Check the how the data looks alike
plot(M[,1],M[,2], pch=".")
plot(M[,10],M[,20], pch=".")

hist(M,1000)

# Reducing dimensionality
variance <- rep(0, num.genes)
for(i in 1:num.genes){
  variance[i] <- sd(M[i,])
}

M <- M[variance>quantile(variance, 0.3),]
num.genes <- length(M[,1])

meanexpr <- rep(0, num.genes)
for(i in 1:num.genes){
  meanexpr[i] <- mean(M[i,])
}
M <- M[meanexpr>quantile(meanexpr, 0.2),]
num.genes <- length(M[,1])

genename <- rownames(M)

hist(M,1000)

# Checking consistency of the data
CC <- rep(0, num.samples**2) # Initialize the correlation matrix
dim(CC) <- c(num.samples,num.samples)

for(i in 1:(num.samples-1)){
  CC[i,i] <- 1
  CC[i+1,i+1] <- 1
  for(j in (i+1):num.samples){
    CC[i,j] <- corr(M[,c(i,j)])
    CC[j,i] <- CC[i,j]
  }
}

samplecol <- rep("black",num.samples)
samplecol[sampleclass=="RPMI"] <- "blue"
samplecol[sampleclass=="Candida"] <- "yellow"
samplecol[sampleclass=="MTB"] <- "red"
heatmap(CC,ColSideColors=samplecol)

# Identification of differentially expressed Genes
p <- rep(0,num.genes)
t <- p
# X is the matrix with the Candida samples.
# Y is the matrix with the control samples.
X <- M[,sampleclass=="Candida"]
Y <- M[,sampleclass=="RPMI"]

for(i in 1:num.genes){
  h <- t.test(X[i,],Y[i,])
  p[i] <- h[[3]]
  t[i] <- h[[1]]
}

p <- p.adjust(p,method="fdr")

r <- rank(-t, ties.method="random")
n <- 1:num.genes
for(i in 1:150){
  g <- n[r==i]
  s <- sprintf("%d %s %3.2e %3.2f",i,genename[g],p[g],t[g])
  print(s)
}

# 5. Gene Set Enrichment Analysis -----------------------------------------

h <-n[r<=500]
genelist <- genename[h]

hGsea <- gost(genelist, organism = "hsapiens",sources = "GO:BP", 
              correction_method = "fdr")

hGseaRes <- hGsea$result
hGseaRes <- filter(hGseaRes, term_size > 20, term_size < 300)

dim(hGseaRes)
hGseaRes[1:5, 11:16]

hGsea <- gost(genelist, organism = "hsapiens",sources = "GO:BP", 
              correction_method = "fdr", evcodes = T)
hGseaRes <- hGsea$result
hGseaRes <- filter(hGseaRes, term_size > 20, term_size < 300)

str(hGseaRes[1,]) # Here we take a look the details of one gene
