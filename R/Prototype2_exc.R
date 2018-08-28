### Load libraries
library(mvnfast)
library(mvtnorm)

### Load Iris Data
Data <- iris
Data <- Data[,-5]
Data <- as.matrix(Data)

### Gain function
Gain <- function(COUNT) {
  0.1^COUNT
}

### Initialization parameters
Groups <- 3
Dim_vec <- dim(Data)
MAX_num <- 1000

### Initialize Pi
Pi_vec <- rep(1/Groups,Groups)
Pi_vec_new <- Pi_vec

### Initialize Mu
Kmeans <- kmeans(Data,Groups)
Mean_list <- list()
for (ii in 1:Groups) {
  Mean_list[[ii]] <- Kmeans$centers[ii,]
}
Mean_list_new <- Mean_list

### Initialize Covariance
Cov_list <- list()
for (ii in 1:Groups) {
  Cov_list[[ii]] <- cov(Data)
}
Cov_list_new <- Cov_list

Tau_fun <- function(X_vec,Mean_list,Cov_list) {
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
T2 <- Mean_list
T3 <- mapply(function(x,y,z){x*z+y%*%t(y)/x},T1,T2,Cov_list)
T3 <- lapply(seq_len(ncol(T3)), function(i) matrix(T3[,i],Dim_vec[2],Dim_vec[2]))
### Algorithm
COUNT <- 1
while(COUNT<=MAX_num) {
  COUNT <- COUNT + 1
  X_vec <- Data[sample(1),]
  Tau <- Tau_fun(X_vec,Mean_list,Cov_list)
  T1 <- (1-Gain(COUNT))*T1+Gain(COUNT)*Tau
  T2 <- lapply(seq_len(Groups),function(ii) {
    (1-Gain(COUNT))*T2[[ii]]+Gain(COUNT)*Tau[ii]*X_vec
  })
  T3 <- lapply(seq_len(Groups),function(ii) {
    (1-Gain(COUNT))*T3[[ii]]+Gain(COUNT)*Tau[ii]*X_vec%*%t(X_vec)
  })
  print(T1)
  print(T2)
  print(T3)
}
Mean_list <- T2
Pi_vec <- T1
Cov_list <- lapply(seq_len(Groups),function(ii){
  (T3[[ii]]-T2[[ii]]%*%t(T2[[ii]])/T1[ii])/T1[ii]
})

Pi_vec
Mean_list
Cov_list