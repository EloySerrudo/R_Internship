# DAY 2

x <- c(-1,0,0,10,1000,33000) # creates a vector called x
y <- c(-2.3,18,220,100,0,1000)

x1 <- x
y1 <- y

x1[x1 < 1] <- 1 # If value of x1 < 1 then set it to 1
y1[y1 < 1] <- 1

plot(x1, y1, log="xy")

# 4 The Xenopus larvae experiment -----------------------------------------

Msignal <- read.table("part2_data/primR2.tab") # read in file
Mbgd <- read.table("part2_data/primbgdR2.tab")

head(Msignal)
str(Msignal)
str(Mbgd)

geneids <- rownames(Msignal) #creates a vector with gene IDs
color <- colnames(Msignal) #creates a vector of colors

geneids[1:5]
color[1:5]

# Scatterplot
Msignal_0 <- Msignal
Msignal[Msignal < 1] <- 1
plot(Msignal[,1], Msignal[,2], log="xy", pch=".", 
     xlab="cy3 (Stage 19)", ylab="cy5 (Stage 10)")
title(main="Signal raw data")
pdf("signalRawData.pdf")

# Normalizing Xenopus Data

#Each entry of Msignal is subtracted with the corresponding entry in Mbgd.
M <- Msignal - Mbgd
M[M < 1] <- 1
plot(M[,1], M[,2], log="xy", pch=".", 
     xlab="cy3 (Stage 19)", ylab="cy5 (Stage 10)")
title(main="Signal raw data without background")

q1 <- quantile(M[,1],0.05)
q2 <- quantile(M[,2],0.05)
q1
q2

# Global background subtraction
M[,1] <- M[,1] - q1
M[,2] <- M[,2] - q2
M[M < 1] <- 1
plot(M[,1], M[,2], log="xy", pch=".", 
     xlab="cy3 (Stage 19)", ylab="cy5 (Stage 10)")
title(main="Signal raw data without background & quantile")

plot(M[,1], M[,2], log="xy", pch=".")
m <- min(max(M[,1]),max(M[,2])) # determine minimum of maximums
m
lines(c(1,m),c(1,m), col="red") # draws the line in my favorite color
title(main="Signal raw data without background & quantile")

# Normalization of multiplicative influences
m1 <- median(M[,1])    #median of cy3 values
m2 <- median(M[,2])    #median of cy5 values
m12 <- mean(m1,m2)     #mean of medians
M[,1] <- M[,1]/m1*m12  #equation to equal out multiplicative factors
M[,2] <- M[,2]/m2*m12

plot(M[,1], M[,2], log="xy", pch=".", 
     xlab="cy3 (Stage 19)", ylab="cy5 (Stage 10)")
m <- min(max(M[,1]),max(M[,2])) # determine minimum of maximums
m
lines(c(1,m),c(1,m), col="red") # draws the line in my favorite color
title(main="Signal raw data without additive & multiplicative influences")

# Adding a constant to smooth the ratios for small signal intensities.
M <- M + 10 #adds 10 to all values in M
plot(M[,1], M[,2], log="xy", pch=".", 
     xlab="cy3 (Stage 19)", ylab="cy5 (Stage 10)")
m <- min(max(M[,1]),max(M[,2]))
m
lines(c(1,m),c(1,m), col="red")
lines(c(1,m),2 * c(1,m), col="yellow")
lines(c(1,m),1/2 * c(1,m), col="yellow")
lines(c(1,m),4 * c(1,m), col="green")
lines(c(1,m),1/4 * c(1,m), col="green")
title(main="Normalized Signal raw data")
#which genes might be differentially expressed?


# 5. VSN Normalization -------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("vsn")

library(vsn)

M0 <- Msignal - Mbgd
M <- as.matrix(M0)
typeof(M)
class(M)

data <- vsn2(M)
typeof(data)
class(data)

help('vsn2')
?vsn2

Mvsn <- exprs(data)
typeof(Mvsn)
class(Mvsn)

str(Mvsn)
head(Mvsn)

plot(Mvsn[,1], Mvsn[,2], pch=".")

# 6. The DESeq2 package ---------------------------------------------------

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")

library(DESeq2)

rawCounts <- read.table("part2_data/SARSCoV2_rawCounts.txt")
sampleTable <- read.table("part2_data/SARSCoV2_sampleInfo.txt")

# the data is normalized and analyzed all at once:
ddsMat <- DESeqDataSetFromMatrix(countData = rawCounts, colData = sampleTable, 
                                 design = ~virus)
dds <- DESeq(ddsMat)

summary(results(dds, alpha = 0.05, lfcThreshold = 0.5))

normCounts <- estimateSizeFactors(ddsMat)
normCounts <- counts(normCounts, normalized= TRUE)
