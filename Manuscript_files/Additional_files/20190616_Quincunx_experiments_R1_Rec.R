#################################################################
##                            Quin1                            ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

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
Mu[7,] <- c(1,1,1)
Mu[8,] <- c(1,-1,-1)
Mu[9,] <- c(1,1,1)
# Covariance matrices
Sigma <- array(NA,c(3,3,9))
for (ii in 1:9) {
  Sigma[,,ii] <- diag(rep(1,3))/64
}
# Set the dimension of the covariance array
dim(Sigma) <- c(d, d, g)
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
Groups <- 9
# Set a random seed
set.seed(20190129)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)
Timing <- matrix(NA,100,5)
ARI_results <- matrix(NA,100,9)
SE_results <- matrix(NA,100,9)

# Conduct simulation study
for (rr in 1:Rep) {
  # Simulate data
  Pre_data <- simdataset(NN,Pi,Mu,Sigma)
  IDs <- Pre_data$id
  Data <- Pre_data$X
  # Randomly generate labels for initialization
  Samp <- sample(1:Groups,NN,replace = T)
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  # Run batch EM algorithm 
  Tick <- proc.time()[3]
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=10))
  Timing[rr,1] <- proc.time()[3]-Tick
  # Get likelihood value for batch EM algorithm
  Results[rr,1] <- MC$loglik
  # Get parameter estimates
  MC_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    MC_matrix[ii,] <- c(MC$parameters$pro[ii],
                        t(MC$parameters$mean)[ii,],
                        MC$parameters$variance$sigma[,,ii][upper.tri(MC$parameters$variance$sigma[,,ii],diag = T)])
  }
  SE_results[rr,1] <- sum(apply((as.matrix(dist(rbind(MC_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI
  ARI_results[rr,1] <- adjustedRandIndex(IDs,apply(MC$z,1,which.max))
  
  # Run minibatch algorithm with batch size 10000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/10000,Groups,0.6,1-10^-10,10000)
  Results[rr,2] <- Sto$`reg_log-likelihood`
  Results[rr,3] <- Sto$`pol_log-likelihood`
  Timing[rr,2] <- proc.time()[3]-Tick
  
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                        t(Sto$reg_means)[ii,],
                        Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,2] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,2] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,3] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,3] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Run minibatch algorithm with batch size 20000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/20000,Groups,0.6,1-10^-10,20000)
  Results[rr,4] <- Sto$`reg_log-likelihood`
  Results[rr,5] <- Sto$`pol_log-likelihood`
  Timing[rr,3] <- proc.time()[3]-Tick
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                         t(Sto$reg_means)[ii,],
                         Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,4] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,4] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,5] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,5] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Run truncated minibatch algorithm with batch size 10000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/10000,Groups,0.6,1-10^-10,10000,
                           1000,1000,1000)
  Results[rr,6] <- Sto$`reg_log-likelihood`
  Results[rr,7] <- Sto$`pol_log-likelihood`
  Timing[rr,4] <- proc.time()[3]-Tick
  
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                         t(Sto$reg_means)[ii,],
                         Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,6] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,6] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,7] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,7] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Run truncated minibatch algorithm with batch size 20000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/20000,Groups,0.6,1-10^-10,20000,
                           1000,1000,1000)
  Results[rr,8] <- Sto$`reg_log-likelihood`
  Results[rr,9] <- Sto$`pol_log-likelihood`
  Timing[rr,5] <- proc.time()[3]-Tick
   
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                         t(Sto$reg_means)[ii,],
                         Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,8] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,8] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,9] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,9] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Save and print outputs
  save(Results,file='./Quin1.rdata')
  print(c(rr,Results[rr,]))
  save(Timing,file='./Quin1timing.rdata')
  save(ARI_results,file='./Quin1ARI.rdata')
  save(SE_results,file='./Quin1SE.rdata')
  print(c(rr,Timing[rr,]))
  print(c(rr,ARI_results[rr,]))
  print(c(rr,SE_results[rr,]))
  
  # Also sink results to a text file
  sink('./Quin1.txt',append = TRUE)
  cat(rr,Results[rr,],'\n')
  sink()
  sink('./Quin1timing.txt',append = TRUE)
  cat(rr,Timing[rr,],'\n')
  sink()
  sink('./Quin1ARI.txt',append = TRUE)
  cat(rr,ARI_results[rr,],'\n')
  sink()
  sink('./Quin1SE.txt',append = TRUE)
  cat(rr,SE_results[rr,],'\n')
  sink()
}

#################################################################
##                           Quin2                             ##
#################################################################

# Load libraries
library('MixSim')
library('mclust')

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

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
Mu[7,] <- c(1,1,1)
Mu[8,] <- c(1,-1,-1)
Mu[9,] <- c(1,1,1)
# Covariance matrices
Sigma <- array(NA,c(3,3,9))
for (ii in 1:9) {
  Sigma[,,ii] <- diag(rep(1,3))/64
}
# Set the dimension of the covariance array
dim(Sigma) <- c(d, d, g)
# True matrix
True_matrix <- matrix(NA,g,1+d+d+choose(d,2))
for (ii in 1:g) {
  True_matrix[ii,] <- c(Pi[ii],Mu[ii,],Sigma[,,ii][upper.tri(Sigma[,,ii],diag = T)])
}

## Setup parameters
# Number of observations to simulation
NN <- 10^7
# Set the number repetitions
Rep <- 100
# Number of components to fit
Groups <- 9
# Set a random seed
# set.seed(20190129)
set.seed(20190703)
# Construct a matrix to store the results
Results <- matrix(NA,100,9)
Timing <- matrix(NA,100,5)
ARI_results <- matrix(NA,100,9)
SE_results <- matrix(NA,100,9)

# Conduct simulation study
for (rr in 79:Rep) {
  # Simulate data
  Pre_data <- simdataset(NN,Pi,Mu,Sigma)
  IDs <- Pre_data$id
  Data <- Pre_data$X
  # Randomly generate labels for initialization
  Samp <- sample(1:Groups,NN,replace = T)
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  # Run batch EM algorithm 
  Tick <- proc.time()[3]
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=10))
  Timing[rr,1] <- proc.time()[3]-Tick
  # Get likelihood value for batch EM algorithm
  Results[rr,1] <- MC$loglik
  # Get parameter estimates
  MC_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    MC_matrix[ii,] <- c(MC$parameters$pro[ii],
                        t(MC$parameters$mean)[ii,],
                        MC$parameters$variance$sigma[,,ii][upper.tri(MC$parameters$variance$sigma[,,ii],diag = T)])
  }
  SE_results[rr,1] <- sum(apply((as.matrix(dist(rbind(MC_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI
  ARI_results[rr,1] <- adjustedRandIndex(IDs,apply(MC$z,1,which.max))
  
  # Run minibatch algorithm with batch size 10000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/100000,Groups,0.6,1-10^-10,100000)
  Results[rr,2] <- Sto$`reg_log-likelihood`
  Results[rr,3] <- Sto$`pol_log-likelihood`
  Timing[rr,2] <- proc.time()[3]-Tick
  
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                         t(Sto$reg_means)[ii,],
                         Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,2] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,2] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,3] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,3] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Run minibatch algorithm with batch size 20000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/200000,Groups,0.6,1-10^-10,200000)
  Results[rr,4] <- Sto$`reg_log-likelihood`
  Results[rr,5] <- Sto$`pol_log-likelihood`
  Timing[rr,3] <- proc.time()[3]-Tick
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                         t(Sto$reg_means)[ii,],
                         Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,4] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,4] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,5] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,5] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Run truncated minibatch algorithm with batch size 10000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/100000,Groups,0.6,1-10^-10,100000,
                           1000,1000,1000)
  Results[rr,6] <- Sto$`reg_log-likelihood`
  Results[rr,7] <- Sto$`pol_log-likelihood`
  Timing[rr,4] <- proc.time()[3]-Tick
  
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                         t(Sto$reg_means)[ii,],
                         Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,6] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,6] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,7] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,7] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Run truncated minibatch algorithm with batch size 20000
  Tick <- proc.time()[3]
  Sto <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                           msEst$parameters$variance$sigma,
                           10*NN/200000,Groups,0.6,1-10^-10,200000,
                           1000,1000,1000)
  Results[rr,8] <- Sto$`reg_log-likelihood`
  Results[rr,9] <- Sto$`pol_log-likelihood`
  Timing[rr,5] <- proc.time()[3]-Tick
  
  # Get parameter estimates (regular)
  Reg_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Reg_matrix[ii,] <- c(Sto$reg_proportions[ii],
                         t(Sto$reg_means)[ii,],
                         Sto$reg_covariances[,,ii][upper.tri(Sto$reg_covariances[,,ii],diag = T)])
  }
  SE_results[rr,8] <- sum(apply((as.matrix(dist(rbind(Reg_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (reg)
  Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
  ARI_results[rr,8] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Get parameter estimates (polyak)
  Pol_matrix <- matrix(NA,g,1+d+d+choose(d,2))
  for (ii in 1:g) {
    Pol_matrix[ii,] <- c(Sto$pol_proportions[ii],
                         t(Sto$pol_means)[ii,],
                         Sto$pol_covariances[,,ii][upper.tri(Sto$pol_covariances[,,ii],diag = T)])
  }
  SE_results[rr,9] <- sum(apply((as.matrix(dist(rbind(Pol_matrix,True_matrix),diag=T,upper=T))[1:g,(g+1):(2*g)])^2,1,min))
  # Get ARI (pol)
  Cluster <- GMM_arma_cluster(t(Data),Sto$pol_proportions,Sto$pol_means,Sto$pol_covariances)
  ARI_results[rr,9] <- adjustedRandIndex(IDs,unlist(Cluster))
  
  # Save and print outputs
  save(Results,file='./Quin2_2.rdata')
  print(c(rr,Results[rr,]))
  save(Timing,file='./Quin2timing_2.rdata')
  save(ARI_results,file='./Quin2ARI_2.rdata')
  save(SE_results,file='./Quin2SE_2.rdata')
  print(c(rr,Timing[rr,]))
  print(c(rr,ARI_results[rr,]))
  print(c(rr,SE_results[rr,]))
  
  # Also sink results to a text file
  sink('./Quin2_2.txt',append = TRUE)
  cat(rr,Results[rr,],'\n')
  sink()
  sink('./Quin2timing_2.txt',append = TRUE)
  cat(rr,Timing[rr,],'\n')
  sink()
  sink('./Quin2ARI_2.txt',append = TRUE)
  cat(rr,ARI_results[rr,],'\n')
  sink()
  sink('./Quin2SE_2.txt',append = TRUE)
  cat(rr,SE_results[rr,],'\n')
  sink()
}