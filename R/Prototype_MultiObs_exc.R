### Load libraries
library(mvnfast)
library(mvtnorm)
library(mclust)
library(EMMIXskew)
library(mixtools)

### Load Iris Data
# Data <- iris
# Data <- Data[,-5]
Data <- as.matrix(Data)

### Gain function
Gain <- function(COUNT) {
  1/COUNT^.6
}

### Initialization parameters
# Groups <- 3
Dim_vec <- dim(Data)
Batch_size <- NN/10000
MAX_num <- (attributes(MC)$info[1]*NN)/Batch_size
Pol_start <- round(MAX_num/2)

### Initialize Pi
Pi_vec <- msEst$parameters$pro

### Initialize Mean
# hcTree <- hcVVV(data = Data)
# cl <- hclass(hcTree, Groups)
msEst <- mstep(modelName = "VVV", data = Data, z = unmap(A$id))
Mean_list <- list()
for (ii in 1:Groups) {
  Mean_list[[ii]] <- msEst$parameters$mean[,ii]
}
### Initialize Covariance
Cov_list <- list()
for (ii in 1:Groups) {
  Cov_list[[ii]] <- msEst$parameters$variance$sigma[,,ii]
}

Refit2 <- mvnormalmixEM(x = Data,
                       lambda = Pi_vec,
                       mu = Mean_list,
                       sigma = Cov_list,
                       k = Groups,
                       maxit = 1)

### Tau functions
Tau_fun <- function(X_vec,Pi_vec,Mean_list,Cov_list) {
  Tau_vec <- c()
  for (ii in 1:Groups) {
    Tau_vec[ii] <- Pi_vec[ii]*dmvnorm(X_vec,
                                      Mean_list[[ii]],
                                      Cov_list[[ii]])
  }
  Tau_vec <- Tau_vec/sum(Tau_vec)
  return(Tau_vec)
}
X_samp <- Data[sample(Dim_vec[1],Batch_size,replace=T),]
Multi_Tau_fun <- function(X_samp,Pi_vec,Mean_list,Cov_list) {
  t(apply(X_samp,1,function(x){Tau_fun(x,Pi_vec,Mean_list,Cov_list)}))
}

T1 <- Pi_vec*Batch_size
T2 <- lapply(1:Groups,function(ii){Mean_list[[ii]]*Pi_vec[ii]})
T3 <- mapply(function(x,y,z){x*z+y%*%t(y)/x},T1,T2,Cov_list)
T3 <- lapply(seq_len(ncol(T3)), function(i) matrix(T3[,i],Dim_vec[2],Dim_vec[2]))
### Algorithm
COUNT <- 0
# Initialize some lists
Run_Pi_list <- list()
Run_Mean_list <- list()
Run_Cov_list <- list()

while(COUNT<MAX_num) {
  
  ### Iterate Count
  COUNT <- COUNT + 1
  
  ### Averaging lists
  Run_Pi_list[[COUNT]] <- Pi_vec
  Run_Mean_list[[COUNT]] <- matrix(unlist(Mean_list),Dim_vec[2],Groups,byrow = F)
  Run_Cov_list[[COUNT]] <- array(unlist(Cov_list),dim=c(Dim_vec[2],Dim_vec[2],Groups))
  
  X_samp <- Data[sample(Dim_vec[1],Batch_size,replace=T),]
  Tau <- Multi_Tau_fun(X_samp,Pi_vec,Mean_list,Cov_list)
  
  T1 <- (1-Gain(COUNT))*T1+Gain(COUNT)*colSums(Tau)
  
  T2 <- lapply(seq_len(Groups),function(ii) {
    (1-Gain(COUNT))*T2[[ii]]+Gain(COUNT)*colSums(Tau[,ii]*X_samp)
  })
  
  T3 <- lapply(seq_len(Groups),function(ii) {
    (1-Gain(COUNT))*T3[[ii]]+Gain(COUNT)*
      Reduce('+',lapply(1:Batch_size,function(jj){Tau[jj,ii]*X_samp[jj,]%*%t(X_samp[jj,])}))
  })
  

  Pi_vec <- T1/Batch_size
  Mean_list <- lapply(1:Groups,function(ii){
    T2[[ii]]/T1[[ii]]
  })
  Cov_list <- lapply(seq_len(Groups),function(ii){
    (T3[[ii]]-T2[[ii]]%*%t(T2[[ii]])/T1[ii])/T1[ii]
  })
  
  print(COUNT)
  # print(Pi_vec)
  # print(Mean_list)
  # print(Cov_list)
}
# 
# Cluster <- sapply(1:Dim_vec[1],function(ii) {
#   Tau_fun(Data[ii,],Mean_list,)
# })

Refit <- mvnormalmixEM(x = Data,
                       lambda = Pi_vec,
                       mu = Mean_list,
                       sigma = Cov_list,
                       k = Groups,
                       maxit = 1)
print(Refit$all.loglik[1])

### Polyak Averaging
Red_Pi_vec <- Reduce('+',Run_Pi_list[-(1:Pol_start)])/(COUNT-Pol_start+1)
Red_Mean_list <- Reduce('+',Run_Mean_list[-(1:Pol_start)])/(COUNT-Pol_start+1)
Red_Mean_list <- lapply(1:Groups,function(ii){Red_Mean_list[,ii]})
Red_Cov_list <- Reduce('+',Run_Cov_list[-(1:Pol_start)])/(COUNT-Pol_start+1)
Red_Cov_list <- lapply(1:Groups,function(ii){Red_Cov_list[,,ii]})
Refit_pol <- mvnormalmixEM(x = Data,
                       lambda = Red_Pi_vec,
                       mu = Red_Mean_list,
                       sigma = Red_Cov_list,
                       k = Groups,
                       maxit = 1)
print(Refit_pol$all.loglik[1])

MC$loglik
print(Refit2$all.loglik[1])