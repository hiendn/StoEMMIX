### Load libraries
library(mvnfast)
library(mvtnorm)

### Load Iris Data
Data <- iris
Data <- Data[,-5]
Data <- as.matrix(Data)

### Gain function
Gain <- function(COUNT) {
  0.1*COUNT^-0.51
}

### Initialization parameters
Groups <- 3
Dim_vec <- dim(Data)
MAX_num <- 3000

### Initialize Pi
Pi_vec <- msEst$parameters$pro

### Initialize Mean
hcTree <- hcVVV(data = Data)
cl <- hclass(hcTree, Groups)
msEst <- mstep(modelName = "VVV", data = Data, z = unmap(cl))
Mean_list <- list()
for (ii in 1:Groups) {
  Mean_list[[ii]] <- msEst$parameters$mean[,ii]
}
### Initialize Covariance
Cov_list <- list()
for (ii in 1:Groups) {
  Cov_list[[ii]] <- msEst$parameters$variance$sigma[,,ii]
}
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

T1 <- Pi_vec
T2 <- lapply(1:Groups,function(ii){Mean_list[[ii]]*Pi_vec[ii]})
T3 <- mapply(function(x,y,z){x*z+y%*%t(y)/x},T1,T2,Cov_list)
T3 <- lapply(seq_len(ncol(T3)), function(i) matrix(T3[,i],Dim_vec[2],Dim_vec[2]))
### Algorithm
COUNT <- 0
while(COUNT<=MAX_num) {
  
  COUNT <- COUNT + 1
  X_vec <- Data[sample(Dim_vec[1],1),]
  Tau <- Tau_fun(X_vec,Pi_vec,Mean_list,Cov_list)
  
  T1 <- (1-Gain(COUNT))*T1+Gain(COUNT)*Tau
  
  T2 <- lapply(seq_len(Groups),function(ii) {
    (1-Gain(COUNT))*T2[[ii]]+Gain(COUNT)*Tau[ii]*X_vec
  })
  
  T3 <- lapply(seq_len(Groups),function(ii) {
    (1-Gain(COUNT))*T3[[ii]]+Gain(COUNT)*Tau[ii]*X_vec%*%t(X_vec)
  })
  

  Pi_vec <- T1
  Mean_list <- lapply(1:Groups,function(ii){
    T2[[ii]]/T1[[ii]]
  })
  Cov_list <- lapply(seq_len(Groups),function(ii){
    (T3[[ii]]-T2[[ii]]%*%t(T2[[ii]])/T1[ii])/T1[ii]
  })
  
  print(COUNT)
  print(Pi_vec)
  print(Mean_list)
  print(Cov_list)
}
# 
# Cluster <- sapply(1:Dim_vec[1],function(ii) {
#   Tau_fun(Data[ii,],Mean_list,)
# })