#################################################################
##                            Wreath1                          ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Load in Wreath data
data(wreath)

## Estimate mixture model parameters
MC_wreath <- Mclust(wreath, G = 14, modeNames = 'VVV')
# Set parameters 
g <- 14
d <-2
# True matrix
True_matrix <- matrix(NA,g,1+d+d+choose(d,2))
for (ii in 1:g) {
  True_matrix[ii,] <- c(Pi[ii],Mu[ii,],Sigma[,,ii][upper.tri(Sigma[,,ii],diag = T)])
}

## Setup parameters
# Number of observations to simulation
NN <- 10^6
# Set the number repetitions
Rep <- 100
# Number of components to fit
Groups <- 14
# Set a random seed
set.seed(20190130)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)

# Conduct simulation study
for (ii in 1:Rep) {
  # Simulate data
  Data <- simVVV(MC_wreath$parameters,NN)
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
  save(Results,file='./Wreath1.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./Wreath1.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}

#################################################################
##                            Wreath2                          ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Load in Wreath data
data(wreath)

## Estimate mixture model parameters
MC_wreath <- Mclust(wreath, G = 14, modeNames = 'VVV')

## Setup parameters
# Number of observations to simulation
NN <- 10^7
# Set the number repetitions
Rep <- 100
# Number of components to fit
Groups <- 14
# Set a random seed
set.seed(20190130)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)

# Conduct simulation study
for (ii in 1:Rep) {
  # Simulate data
  Data <- simVVV(MC_wreath$parameters,NN)
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
  save(Results,file='./Wreath2.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./Wreath2.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}