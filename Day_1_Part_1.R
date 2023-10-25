# DAY 1

x <- 5

ls() # Lists all objects in the buffer.
str(x) # Shows the structure of an object
class(x) # Shows the class of an object
x # Shows the content of x
rm(x) # Deletes object
rm(x, y) # Deletes several objects
rm(list=ls()) # Deletes all objects (BE CAREFUL!)

getwd()# get the working directory)
setwd("/documents/")# set working directory)
source("test.R") # Calls and executes the script

help.start() # Calls the html-page via the console.


# 8. Scalars --------------------------------------------------------------

seq(3, 51, by=3)
x <- 20:1
x[4] <- 100

# 9. Plotting -------------------------------------------------------------

x <- seq(-pi,pi, by=0.2)
y <- sin(x)
z <- cos(x)

plot(x, y) # The plot function plot() is very versatile!
plot(x, y, pch=".") #The parameter pch defines the used sign.
plot(x, y, pch="+")

plot(x, y, pch='x', col="red")
plot(x, y, pch='x', col=rainbow(10))

lines(x,y,col="green")
points(x, z, col='blue')

# ******** What is the difference between plot and points? ********

plot(x,y,xlim=c(0,2)) # x- axis from 0 to 2
plot(x,y,xlim=c(0,2),ylim=c(-0.5,2))

# 9.1 Labelling plots -----------------------------------------------------

plot (x, y, xlab = "x-axis", ylab = "sinus", main = "The Sinus Curve")

# You CAN also choose a header AFTER PLOTTING the plot.
title(main="The Sinus Curve")
title(sub="Title below the headline")

# 9.2 Printing and generating output files from a plot --------------------

plot(x,y,xlab="x axis",ylab="sinus",main="The sinus curve")
png("sinuscurve.png")
pdf("sinuscurve.pdf")
dev.off()

# ******** Try different formats like png, jpeg and pdf ********
# Test what happens when one puts several plots between
# the commands pdf(filename.pdf) and dev.off().

# 10. Loops ---------------------------------------------------------------

x <- rnorm(1) # generates 1 random number around 0 from a normal distribution

#### TASTK 1:

x <- sample(1:10, 3)

sort_vector <- function(x) {
  if (x[1] > x[2]) {
    tmp <- x[2]
    x[2] <- x[1]
    x[1] <- tmp
  } 
  if (x[1] > x[3]) {
    tmp <- x[3]
    x[3] <- x[1]
    x[1] <- tmp
  }
  if (x[2] > x[3]) {
    tmp <- x[3]
    x[3] <- x[2]
    x[2] <- tmp
  }
  
  x
}
sort_vector(x)

#### TASTK 2:

X_vector <- sample(1:100, 10)
is_even_or_odd <- function(X_vector, x) {
  for (x in X_vector) {
    cat(x)
    if (x %% 2 == 0) {
      cat(" is even\n")
    } else {
      cat(" is odd\n")
    }
  }
}

#### TASK 3:

i = 1
while(i <= 10){
  print(i^2)
  i <- i + 1
}

# 11. Arrays --------------------------------------------------------------

x = c(2,5,8,3,8) # We create 3 1D vectors
y = c(1,1,1,1,1)
z = c(5,4,3,2,1)
M = c(x,y,z)     # Matrix is builted from vectors x, y and z
print(M)
# If the matrix ”doesn't know“ how many rows and columns it has, it
# will automatically be one-dimensional and the vectors x, y and z will simple
# be aligned on another.

dim(M) <- c(5,3)# You can change that by defining the dimension of M
print(M)# 5x3 matrix (5 rows, 3 columns)
# The elements of a matrix will be filled from top to bottom, from left to
# right, so x, y and z become column-vectors.

mean(z)
max(z)
min(z)

# 12. Bonus exercise ------------------------------------------------------

table <- read.table("table1.txt") # Is this a data frame?

M <- as.matrix(table)
plot(M[,1], M[,2])