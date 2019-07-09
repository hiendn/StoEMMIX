#################################################################
##                            Poi1                            ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190708_ExpPois_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

## Extract various required variables
# Get the number of subpopulations
g <- 3
## Estimate mixture model parameters
# Proportions
Pi <- c(0.8,0.1,0.1)
# Lambda
Lambda <- c(1,5,12)
# True matrix
True_matrix <- matrix(NA,g,2)
for (ii in 1:g) {
  True_matrix[ii,] <- c(Pi[ii],Lambda[ii])
}

## Setup parameters
# Number of observations to simulation
NN <- 10^6
# Set the number repetitions
Rep <- 100
# Number of components to fit
Groups <- 3
# Set a random seed
set.seed(20190708)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)
Timing <- matrix(NA,100,5)
ARI_results <- matrix(NA,100,9)
SE_results <- matrix(NA,100,9)

# Conduct simulation study
for (rr in 1:Rep) {
  
  # Simulate data
  IDs <- sample(1:3,NN,replace=TRUE,prob=Pi)
  Data <- matrix(NA,nrow=NN,ncol=1)
  Data[IDs==1] <- rpois(sum(IDs==1),Lambda[1])
  Data[IDs==2] <- rpois(sum(IDs==2),Lambda[2])
  Data[IDs==3] <- rpois(sum(IDs==3),Lambda[3])
  
  # Randomly generate labels for initialization
  Samp <- sample(1:Groups,NN,replace = T)
  
  # Initialize parameters
  Pi_int <- table(Samp)/NN
  Lambda_int <- c(mean(Data[Samp==1]),mean(Data[Samp==2]),mean(Data[Samp==3]))
  # Run batch EM algorithm 
  Tick <- proc.time()[3]
  MC <- EMPoisson_pol(data=Data,pi_r=Pi_int,lambda_r = Lambda_int,maxit_r = 10,groups_r = g)
  Timing[rr,1] <- proc.time()[3]-Tick
  # Get likelihood value for batch EM algorithm
  Results[rr,1] <-MC$`reg_log-likelihood`
  # Get parameter estimates
  MC_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    MC_matrix[ii,] <- c(MC$reg_proportions[ii],MC$reg_lambda[ii])
  }
  SE_results[rr,1] <- sum(apply((as.matrix(dist(rbind(MC_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI
  ARI_results[rr,1] <- adjustedRandIndex(IDs,Pois_clust(Data,MC$reg_proportions,MC$reg_lambda))
  
  # Run minibatch algorithm with batch size 10000
  Tick <- proc.time()[3]
  Sto <- stoPoisson_pol(Data, Pi_int, Lambda_int,
                      10*NN/10000,Groups,0.6,1-10^-10,10000)
  Results[rr,2] <- Sto$`reg_log-likelihood`
  Results[rr,3] <- Sto$`pol_log-likelihood`
  Timing[rr,2] <- proc.time()[3]-Tick
  # Get parameter estimates
  Reg_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],Sto$reg_lambda[ii])
  }
  SE_results[rr,2] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  ARI_results[rr,2] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$reg_proportions,Sto$reg_lambda))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],Sto$pol_lambda[ii])
  }
  SE_results[rr,3] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  ARI_results[rr,3] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$pol_proportions,Sto$pol_lambda))
  
  # Run minibatch algorithm with batch size 20000
  Tick <- proc.time()[3]
  Sto <- stoPoisson_pol(Data, Pi_int, Lambda_int,
                      10*NN/20000,Groups,0.6,1-10^-10,20000)
  Results[rr,4] <- Sto$`reg_log-likelihood`
  Results[rr,5] <- Sto$`pol_log-likelihood`
  Timing[rr,3] <- proc.time()[3]-Tick
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],Sto$reg_lambda[ii])
  }
  SE_results[rr,4] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  ARI_results[rr,4] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$reg_proportions,Sto$reg_lambda))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],Sto$pol_lambda[ii])
  }
  SE_results[rr,5] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  ARI_results[rr,5] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$pol_proportions,Sto$pol_lambda))
  
  # Save and print outputs
  save(Results,file='./Poi1.rdata')
  print(c(rr,Results[rr,]))
  save(Timing,file='./Poi1timing.rdata')
  save(ARI_results,file='./Poi1ARI.rdata')
  save(SE_results,file='./Poi1SE.rdata')
  print(c(rr,Timing[rr,]))
  print(c(rr,ARI_results[rr,]))
  print(c(rr,SE_results[rr,]))
  
  # Also sink results to a text file
  sink('./Poi1.txt',append = TRUE)
  cat(rr,Results[rr,],'\n')
  sink()
  sink('./Poi1timing.txt',append = TRUE)
  cat(rr,Timing[rr,],'\n')
  sink()
  sink('./Poi1ARI.txt',append = TRUE)
  cat(rr,ARI_results[rr,],'\n')
  sink()
  sink('./Poi1SE.txt',append = TRUE)
  cat(rr,SE_results[rr,],'\n')
  sink()
}

#################################################################
##                           Poi2                             ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190708_ExpPois_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

## Extract various required variables
# Get the number of subpopulations
g <- 3
## Estimate mixture model parameters
# Proportions
Pi <- c(0.8,0.1,0.1)
# Lambda
Lambda <- c(1,5,12)
# True matrix
True_matrix <- matrix(NA,g,2)
for (ii in 1:g) {
  True_matrix[ii,] <- c(Pi[ii],Lambda[ii])
}

## Setup parameters
# Number of observations to simulation
NN <- 10^7
# Set the number repetitions
Rep <- 100
# Number of components to fit
Groups <- 3
# Set a random seed
set.seed(20190708)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)
Timing <- matrix(NA,100,5)
ARI_results <- matrix(NA,100,9)
SE_results <- matrix(NA,100,9)

# Conduct simulation study
for (rr in 1:Rep) {
  
  # Simulate data
  IDs <- sample(1:3,NN,replace=TRUE,prob=Pi)
  Data <- matrix(NA,nrow=NN,ncol=1)
  Data[IDs==1] <- rpois(sum(IDs==1),Lambda[1])
  Data[IDs==2] <- rpois(sum(IDs==2),Lambda[2])
  Data[IDs==3] <- rpois(sum(IDs==3),Lambda[3])
  
  # Randomly generate labels for initialization
  Samp <- sample(1:Groups,NN,replace = T)
  
  # Initialize parameters
  Pi_int <- table(Samp)/NN
  Lambda_int <- c(mean(Data[Samp==1]),mean(Data[Samp==2]),mean(Data[Samp==3]))
  # Run batch EM algorithm 
  Tick <- proc.time()[3]
  MC <- EMPoisson_pol(data=Data,pi_r=Pi_int,lambda_r = Lambda_int,maxit_r = 10,groups_r = g)
  Timing[rr,1] <- proc.time()[3]-Tick
  # Get likelihood value for batch EM algorithm
  Results[rr,1] <-MC$`reg_log-likelihood`
  # Get parameter estimates
  MC_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    MC_matrix[ii,] <- c(MC$reg_proportions[ii],MC$reg_lambda[ii])
  }
  SE_results[rr,1] <- sum(apply((as.matrix(dist(rbind(MC_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI
  ARI_results[rr,1] <- adjustedRandIndex(IDs,Pois_clust(Data,MC$reg_proportions,MC$reg_lambda))
  
  # Run minibatch algorithm with batch size 10000
  Tick <- proc.time()[3]
  Sto <- stoPoisson_pol(Data, Pi_int, Lambda_int,
                            10*NN/100000,Groups,0.6,1-10^-10,100000)
  Results[rr,2] <- Sto$`reg_log-likelihood`
  Results[rr,3] <- Sto$`pol_log-likelihood`
  Timing[rr,2] <- proc.time()[3]-Tick
  # Get parameter estimates
  Reg_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],Sto$reg_lambda[ii])
  }
  SE_results[rr,2] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  ARI_results[rr,2] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$reg_proportions,Sto$reg_lambda))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],Sto$pol_lambda[ii])
  }
  SE_results[rr,3] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  ARI_results[rr,3] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$pol_proportions,Sto$pol_lambda))
  
  # Run minibatch algorithm with batch size 20000
  Tick <- proc.time()[3]
  Sto <- stoPoisson_pol(Data, Pi_int, Lambda_int,
                            10*NN/200000,Groups,0.6,1-10^-10,200000)
  Results[rr,4] <- Sto$`reg_log-likelihood`
  Results[rr,5] <- Sto$`pol_log-likelihood`
  Timing[rr,3] <- proc.time()[3]-Tick
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],Sto$reg_lambda[ii])
  }
  SE_results[rr,4] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  ARI_results[rr,4] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$reg_proportions,Sto$reg_lambda))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,2)
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],Sto$pol_lambda[ii])
  }
  SE_results[rr,5] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  ARI_results[rr,5] <- adjustedRandIndex(IDs,Pois_clust(Data,Sto$pol_proportions,Sto$pol_lambda))
  
  # Save and print outputs
  save(Results,file='./Poi2.rdata')
  print(c(rr,Results[rr,]))
  save(Timing,file='./Poi2timing.rdata')
  save(ARI_results,file='./Poi2ARI.rdata')
  save(SE_results,file='./Poi2SE.rdata')
  print(c(rr,Timing[rr,]))
  print(c(rr,ARI_results[rr,]))
  print(c(rr,SE_results[rr,]))
  
  # Also sink results to a text file
  sink('./Poi2.txt',append = TRUE)
  cat(rr,Results[rr,],'\n')
  sink()
  sink('./Poi2timing.txt',append = TRUE)
  cat(rr,Timing[rr,],'\n')
  sink()
  sink('./Poi2ARI.txt',append = TRUE)
  cat(rr,ARI_results[rr,],'\n')
  sink()
  sink('./Poi2SE.txt',append = TRUE)
  cat(rr,SE_results[rr,],'\n')
  sink()
}