##################################################################
##                             Flea                             ##
##################################################################

library(tourr)
library(colorspace)
data(flea)
plot(flea[,1:6],col=rainbow_hcl(3)[flea[,7]],pch=as.numeric(flea[,7]))

##################################################################
##                             ELKI                             ##
##################################################################

library(colorspace)
library('MixSim')
library('mclust')

## Extract various required variables
# Get dimensions
d <- 2
# Get the number of subpopulations
g <- 3

## Estimate mixture model parameters
# Proportions
Pi <- c(5/10,3/10,2/10)
# Mean vectors
Mu <- matrix(NA,3,2)
Mu[1,] <- c(0.3,0.3)
Mu[2,] <- c(0.85,0.35)
Mu[3,] <- c(0.45,0.85)
# Covariance matrices
Sigma <- array(NA,c(2,2,3))
Sigma[,,1] <- diag(c(0.09,0.09)^2)
Sigma[,,2] <- diag(c(0.05,0.1)^2)
Sigma[,,3] <- diag(c(0.035,0.035)^2)

NN <- 1000

set.seed(20190129)
# Simulate data
Pre_data <- simdataset(NN,Pi,Mu,Sigma)
IDs <- Pre_data$id
Data <- Pre_data$X

plot(Data,col=rainbow_hcl(3)[IDs],xlab=expression(y[1]),ylab=expression(y[2]),lwd=2)
grid()

##################################################################
##                           Quincunx                           ##
##################################################################

# Load libraries
library('MixSim')
library('mclust')
library(colorspace)

## Extract various required variables
# Get dimensions
d <- 3
# Get the number of subpopulations
g <- 9
## Estimate mixture model parameters
# Proportions
Pi <- c(2,1,1,1,1,1,1,1,1)/10
# Mean vectors
Mu <- matrix(NA,9,3)
Mu[1,] <- c(0,0,0)
Mu[2,] <- c(-1,-1,-1)
Mu[3,] <- c(-1,-1,1)
Mu[4,] <- c(-1,1,-1)
Mu[5,] <- c(-1,1,1)
Mu[6,] <- c(1,-1,-1)
Mu[7,] <- c(1,-1,1)
Mu[8,] <- c(1,1,-1)
Mu[9,] <- c(1,1,1)
# Covariance matrices
Sigma <- array(NA,c(3,3,9))
for (ii in 1:9) {
  Sigma[,,ii] <- diag(rep(1,3))/64
}

NN <- 1000

set.seed(20190129)
Pre_data <- simdataset(NN,Pi,Mu,Sigma)
IDs <- Pre_data$id
Data <- Pre_data$X
Data <- as.data.frame(Data)
names(Data) <- c('y1','y2','y3')
SAMPLE <- sample(1:nrow(Data))
plot(Data[SAMPLE,],col=rainbow_hcl(9)[IDs[SAMPLE]])
