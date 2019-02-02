##################################################################
##                        MNIST PC10 G10                        ##
##################################################################

# Load libraries
library(mclust)

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Set a random seed
set.seed(20190202)

Results <- matrix(NA, 100, 7)
for (ii in 1:100) {
  
  # Set number of PCs
  dPC <- 10
  
  # Set number of groups
  Groups <- 10
  
  # Declare number of epochs
  Epoch <- 10
  
  # Set Data to be PCA of dimensions dPC
  Data <- PCA$scores[,1:dPC]
  
  # Sample starting allocation
  Samp <- sample(1:10,70000,replace = T)
  
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  
  # Run Mclust
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=Epoch))
  
  # Estimate the Mixture model using minibatch with truncation
  Sto_trunc <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                                 msEst$parameters$variance$sigma,
                                 Epoch*dim(Data)[1]/7000,Groups,0.6,1-10^-10,7000,
                                 1000,1000,1000)
  
  # Conduct a K-means for comparison
  KM <- kmeans(Data,centers = Groups,iter.max = Epoch, nstart = 1)
  
  # Obtain clustering outcomes
  Cluster_reg <- GMM_arma_cluster(t(Data),Sto_trunc$reg_proportions,
                                  Sto_trunc$reg_means,
                                  Sto_trunc$reg_covariances)
  Cluster_pol <- GMM_arma_cluster(t(Data),Sto_trunc$pol_proportions,
                                  Sto_trunc$pol_means,
                                  Sto_trunc$pol_covariances)
  
  # Compute ARIs
  Results[ii,1] <- adjustedRandIndex(apply(MC$z,1,which.max), c(train_label,test_label))
  Results[ii,2] <- adjustedRandIndex(Cluster_reg$Cluster,c(train_label,test_label))
  Results[ii,3] <- adjustedRandIndex(Cluster_pol$Cluster,c(train_label,test_label))
  Results[ii,4] <- adjustedRandIndex(KM$cluster,c(train_label,test_label))
  
  # Likelihoods
  Results[ii,5] <- MC$loglik
  Results[ii,6] <- Sto_trunc$`reg_log-likelihood`
  Results[ii,7] <- Sto_trunc$`pol_log-likelihood`
  
  # Save and print outputs
  save(Results,file='./MNIST_PC10_G10.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./MNIST_PC10_G10.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}

##################################################################
##                        MNIST PC20 G10                        ##
##################################################################

# Load libraries
library(mclust)

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Set a random seed
set.seed(20190202)

Results <- matrix(NA, 100, 7)
for (ii in 1:100) {
  
  # Set number of PCs
  dPC <- 20
  
  # Set number of groups
  Groups <- 10
  
  # Declare number of epochs
  Epoch <- 10
  
  # Set Data to be PCA of dimensions dPC
  Data <- PCA$scores[,1:dPC]
  
  # Sample starting allocation
  Samp <- sample(1:10,70000,replace = T)
  
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  
  # Run Mclust
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=Epoch))
  
  # Estimate the Mixture model using minibatch with truncation
  Sto_trunc <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                                 msEst$parameters$variance$sigma,
                                 Epoch*dim(Data)[1]/7000,Groups,0.6,1-10^-10,7000,
                                 1000,1000,1000)
  
  # Conduct a K-means for comparison
  KM <- kmeans(Data,centers = Groups,iter.max = Epoch, nstart = 1)
  
  # Obtain clustering outcomes
  Cluster_reg <- GMM_arma_cluster(t(Data),Sto_trunc$reg_proportions,
                                  Sto_trunc$reg_means,
                                  Sto_trunc$reg_covariances)
  Cluster_pol <- GMM_arma_cluster(t(Data),Sto_trunc$pol_proportions,
                                  Sto_trunc$pol_means,
                                  Sto_trunc$pol_covariances)
  
  # Compute ARIs
  Results[ii,1] <- adjustedRandIndex(apply(MC$z,1,which.max), c(train_label,test_label))
  Results[ii,2] <- adjustedRandIndex(Cluster_reg$Cluster,c(train_label,test_label))
  Results[ii,3] <- adjustedRandIndex(Cluster_pol$Cluster,c(train_label,test_label))
  Results[ii,4] <- adjustedRandIndex(KM$cluster,c(train_label,test_label))
  
  # Likelihoods
  Results[ii,5] <- MC$loglik
  Results[ii,6] <- Sto_trunc$`reg_log-likelihood`
  Results[ii,7] <- Sto_trunc$`pol_log-likelihood`
  
  # Save and print outputs
  save(Results,file='./MNIST_PC20_G10.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./MNIST_PC20_G10.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}

##################################################################
##                        MNIST PC50 G10                        ##
##################################################################

# Load libraries
library(mclust)

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Set a random seed
set.seed(20190202)

Results <- matrix(NA, 100, 7)
for (ii in 1:100) {
  
  # Set number of PCs
  dPC <- 50
  
  # Set number of groups
  Groups <- 10
  
  # Declare number of epochs
  Epoch <- 10
  
  # Set Data to be PCA of dimensions dPC
  Data <- PCA$scores[,1:dPC]
  
  # Sample starting allocation
  Samp <- sample(1:10,70000,replace = T)
  
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  
  # Run Mclust
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=Epoch))
  
  # Estimate the Mixture model using minibatch with truncation
  Sto_trunc <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                                 msEst$parameters$variance$sigma,
                                 Epoch*dim(Data)[1]/7000,Groups,0.6,1-10^-10,7000,
                                 1000,1000,1000)
  
  # Conduct a K-means for comparison
  KM <- kmeans(Data,centers = Groups,iter.max = Epoch, nstart = 1)
  
  # Obtain clustering outcomes
  Cluster_reg <- GMM_arma_cluster(t(Data),Sto_trunc$reg_proportions,
                                  Sto_trunc$reg_means,
                                  Sto_trunc$reg_covariances)
  Cluster_pol <- GMM_arma_cluster(t(Data),Sto_trunc$pol_proportions,
                                  Sto_trunc$pol_means,
                                  Sto_trunc$pol_covariances)
  
  # Compute ARIs
  Results[ii,1] <- adjustedRandIndex(apply(MC$z,1,which.max), c(train_label,test_label))
  Results[ii,2] <- adjustedRandIndex(Cluster_reg$Cluster,c(train_label,test_label))
  Results[ii,3] <- adjustedRandIndex(Cluster_pol$Cluster,c(train_label,test_label))
  Results[ii,4] <- adjustedRandIndex(KM$cluster,c(train_label,test_label))
  
  # Likelihoods
  Results[ii,5] <- MC$loglik
  Results[ii,6] <- Sto_trunc$`reg_log-likelihood`
  Results[ii,7] <- Sto_trunc$`pol_log-likelihood`
  
  # Save and print outputs
  save(Results,file='./MNIST_PC50_G10.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./MNIST_PC50_G10.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}

##################################################################
##                        MNIST PC100 G10                       ##
##################################################################

# Load libraries
library(mclust)

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Set a random seed
set.seed(20190203)

Results <- matrix(NA, 100, 7)
for (ii in 1:100) {
  
  # Set number of PCs
  dPC <- 100
  
  # Set number of groups
  Groups <- 10
  
  # Declare number of epochs
  Epoch <- 10
  
  # Set Data to be PCA of dimensions dPC
  Data <- PCA$scores[,1:dPC]
  
  # Sample starting allocation
  Samp <- sample(1:10,70000,replace = T)
  
  # Initialize parameters
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  
  # Run Mclust
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=Epoch))
  
  # Estimate the Mixture model using minibatch with truncation
  Sto_trunc <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                                 msEst$parameters$variance$sigma,
                                 Epoch*dim(Data)[1]/7000,Groups,0.6,1-10^-10,7000,
                                 1000,1000,1000)
  
  # Conduct a K-means for comparison
  KM <- kmeans(Data,centers = Groups,iter.max = Epoch, nstart = 1)
  
  # Obtain clustering outcomes
  Cluster_reg <- GMM_arma_cluster(t(Data),Sto_trunc$reg_proportions,
                                  Sto_trunc$reg_means,
                                  Sto_trunc$reg_covariances)
  Cluster_pol <- GMM_arma_cluster(t(Data),Sto_trunc$pol_proportions,
                                  Sto_trunc$pol_means,
                                  Sto_trunc$pol_covariances)
  
  # Compute ARIs
  Results[ii,1] <- adjustedRandIndex(apply(MC$z,1,which.max), c(train_label,test_label))
  Results[ii,2] <- adjustedRandIndex(Cluster_reg$Cluster,c(train_label,test_label))
  Results[ii,3] <- adjustedRandIndex(Cluster_pol$Cluster,c(train_label,test_label))
  Results[ii,4] <- adjustedRandIndex(KM$cluster,c(train_label,test_label))
  
  # Likelihoods
  Results[ii,5] <- MC$loglik
  Results[ii,6] <- Sto_trunc$`reg_log-likelihood`
  Results[ii,7] <- Sto_trunc$`pol_log-likelihood`
  
  # Save and print outputs
  save(Results,file='./MNIST_PC100_G10.rdata')
  print(c(ii,Results[ii,]))
  # Also sink results to a text file
  sink('./MNIST_PC100_G10.txt',append = TRUE)
  cat(ii,Results[ii,],'\n')
  sink()
}

