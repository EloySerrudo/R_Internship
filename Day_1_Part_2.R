# Day 1

x <- 1:20
length(x)
y <- x^2
plot(x, y, pch='x', col='red')
title(main="Exponential curve")

v <- -20:20
w <- exp(v)


# 13. Arrays --------------------------------------------------------------

x = c(2,5,8,3,8)
y = c(1,1,1,1,1)
z = c(5,4,3,2,1)
M = c(x,y,z)
dim(M) <- c(5,3)

smallM = M[1:4,2:3]
shortz2 = M[2,1:2]
small3 = M[2:5,3]

smallM
shortz2
small3

plot(M[,1], M[,3])
plot(M[1,], M[2,])

# 14. Tables --------------------------------------------------------------

DF <- read.table("table1.txt")
M <- as.matrix(DF)
str(DF) # Structure of the object
M
plot(M[,1],M[,2])

rownames <- row.names(DF)
columnnames <- names(DF)
rownames
columnnames

M <- matrix(M, dimnames=NULL, nrow=nrow(M), ncol=ncol(M))
M

# 15. Logarithmical plotting ----------------------------------------------

x <- seq(1,5, by=0.5) # generate a sequence from 1 to 5 in 0.5 steps
x <- 10^x
y <- x
z <- 2*x
plot(x,y)
lines(x,y, col="red")
lines(x,z, col="green")

plot(x,y, log="xy")
lines(x,y, col="red")
lines(x,z, col="green")

#### Task 2
y1 <- 1/4*x
y2 <- 1/2*x
y3 <- x
y4 <- 2*x
y5 <- 4*x

x0 <- c(3e+04, 4e+04)
y0 <- c(2e+04, 3e+04)
points(x0,y0)
lines(x,y1, col="yellow")
lines(x,y2, col="green")
lines(x,y3, col="red")
lines(x,y4, col="green")
lines(x,y5, col="yellow")

#### Data

DF <- read.table("primRshort.tab")
str(DF)
M <- as.matrix(DF)
str(M)

geneid <- row.names(M)
geneid

plot(M[,1],M[,2])
plot(M[,1]+1,M[,2]+1,log="xy") # Why did we add +1?

lines (c(1,100000), c(1,100000), col="red")
lines (c(1,100000), 2*c(1,100000), col="green")
lines (2*c(1,100000), c(1,100000), col="green")

x <- M[,1]+1
y <- M[,2]+1
while(1){
  res <- identify(x, y, n=1, pos=FALSE, plot=FALSE);
  print(geneid[[res]])
}

# 16.Statistical functions ------------------------------------------------

x <- c(1,3,5,4,2,3,3,3,3,3,3,3,1,7,4)
sort(x)
x
hist(x)
hist(M)
hist(log(M+1))

x <- rnorm(1000)
hist(x)

hist(rnorm(100000))

x = 1:10
mean(x)
sd(x)
median(x)
quantile(x,0.05) # the 5% quantile
# the value at which only 5 % of all data is above

# median seem to be more suitable against out-liers

print(x)
length(x) # should be 10, since you have 10 values
y=sort(x) # sort x by size
y[5:6]    # calculate the middle value of both values
mean(y[5:6]) # return the middle value
median(x) # is it the same as the median?

mean(M)
median(M)

mean(log(M+1))
median(log(M+1))

# Write data from R into a table ------------------------------------------

M = 1:20
dim(M) = c(4,5)
columns = c("c1","c2","c3","c4","c5")
rows = c("r1","r2","r3","r4")
fr = as.data.frame(M)
row.names(fr) = rows
names(fr) = columns
write.table(fr, "myTable.tab", sep="\t", quote=FALSE, col.names=NA)