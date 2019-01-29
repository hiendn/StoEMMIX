#################################################################
##                            Iris1                            ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Load in Iris data
data("iris", package = "datasets")

## Extract various required variables
# Get dimensions
d <- ncol(iris) - 1
# Get species label
id <- as.integer(iris[, 5])
# Get the number of subpopulations
g <- max(id)

## Estimate mixture model parameters
# Proportions
Pi <- prop.table(tabulate(id))
# Mean vectors
Mu <- t(sapply(1:g, function(k){ colMeans(iris[id == k, -5]) }))
# Covariance matrices
Sigma <- sapply(1:g, function(k){ var(iris[id == k, -5]) })
# Set the dimension of the covariance array
dim(Sigma) <- c(d, d, g)

## Setup parameters
# Number of observations to simulation
NN <- 10^6
# Set the number repetitions
Rep <- 100
# Number of components to fit
Groups <- 3
# Set a random seed
set.seed(20190129)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)

# Conduct simulation study
for (ii in 1:Rep) {
  # Simulate data
  Data <- simdataset(NN,Pi,Mu,Sigma)$X
  # Randomly generate labels for initialization
  Samp <- sample(1:Groups,NN,replace = T)
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  # Run batch EM algorithm 
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=10))
  # Get likelihood value for batch EM algorithm
  Results[ii,1] <- MC$loglik
  # Run minibatch algorithm with batch size 10000
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/10000,Groups,0.6,1-10^-10,10000)
  Results[ii,2] <- Sto$`reg_log-likelihood`
  Results[ii,3] <- Sto$`pol_log-likelihood`
  # Run minibatch algorithm with batch size 20000
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/20000,Groups,0.6,1-10^-10,20000)
  Results[ii,4] <- Sto$`reg_log-likelihood`
  Results[ii,5] <- Sto$`pol_log-likelihood`
  # Run truncated minibatch algorithm with batch size 10000
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/10000,Groups,0.6,1-10^-10,10000,
                           1000,1000,1000)
  Results[ii,6] <- Sto$`reg_log-likelihood`
  Results[ii,7] <- Sto$`pol_log-likelihood`
  # Run truncated minibatch algorithm with batch size 20000
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/20000,Groups,0.6,1-10^-10,20000,
                           1000,1000,1000)
  Results[ii,8] <- Sto$`reg_log-likelihood`
  Results[ii,9] <- Sto$`pol_log-likelihood`
  
  # Save and print outputs
  save(Results,file='./Iris1.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./Iris1.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}

#################################################################
##                            Iris2                            ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Load in Iris data
data("iris", package = "datasets")

## Extract various required variables
# Get dimensions
d <- ncol(iris) - 1
# Get species label
id <- as.integer(iris[, 5])
# Get the number of subpopulations
g <- max(id)

## Estimate mixture model parameters
# Proportions
Pi <- prop.table(tabulate(id))
# Mean vectors
Mu <- t(sapply(1:g, function(k){ colMeans(iris[id == k, -5]) }))
# Covariance matrices
Sigma <- sapply(1:g, function(k){ var(iris[id == k, -5]) })
# Set the dimension of the covariance array
dim(Sigma) <- c(d, d, g)

## Setup parameters
# Number of observations to simulation
NN <- 10^7
# Set the number repetitions
Rep <- 100
# Number of components to fit
Groups <- 3
# Set a random seed
set.seed(20190129)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)

# Conduct simulation study
for (ii in 1:Rep) {
  # Simulate data
  Data <- simdataset(NN,Pi,Mu,Sigma)$X
  # Randomly generate labels for initialization
  Samp <- sample(1:Groups,NN,replace = T)
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  # Run batch EM algorithm 
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=10))
  # Get likelihood value for batch EM algorithm
  Results[ii,1] <- MC$loglik
  # Run minibatch algorithm with batch size 10000
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/100000,Groups,0.6,1-10^-10,100000)
  Results[ii,2] <- Sto$`reg_log-likelihood`
  Results[ii,3] <- Sto$`pol_log-likelihood`
  # Run minibatch algorithm with batch size 20000
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/200000,Groups,0.6,1-10^-10,200000)
  Results[ii,4] <- Sto$`reg_log-likelihood`
  Results[ii,5] <- Sto$`pol_log-likelihood`
  # Run truncated minibatch algorithm with batch size 10000
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/100000,Groups,0.6,1-10^-10,100000,
                           1000,1000,1000)
  Results[ii,6] <- Sto$`reg_log-likelihood`
  Results[ii,7] <- Sto$`pol_log-likelihood`
  # Run truncated minibatch algorithm with batch size 20000
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/200000,Groups,0.6,1-10^-10,200000,
                           1000,1000,1000)
  Results[ii,8] <- Sto$`reg_log-likelihood`
  Results[ii,9] <- Sto$`pol_log-likelihood`
  
  # Save and print outputs
  save(Results,file='./Iris2.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./Iris2.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}