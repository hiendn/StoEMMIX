### Load libraries
library(mvnfast)

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
MAX_num <- 10000

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

### Algorithm
COUNT <- 0
while(COUNT<=MAX_num) {
  COUNT <- COUNT + 1

  X_vec <- Data[sample(1),]

  Tau_vec <- c()
  for (ii in 1:Groups) {
    Tau_vec[ii] <- Pi_vec[ii]*dmvnorm(X_vec,
                                   Mean_list[[ii]],
                                   Cov_list[[ii]])
  }
  Tau_vec <- Tau_vec/sum(Tau_vec)
  Pi_vec_new <- (1-Gain(COUNT))*Pi_vec+Gain(COUNT)*Tau_vec

  for (ii in 1:Groups) {
    Mean_list_new[[ii]] <- (1-Gain(COUNT))*Pi_vec[ii]*Mean_list[[ii]] +
      Gain(COUNT)*Tau_vec[ii]*X_vec
    Mean_list_new[[ii]] <- Mean_list_new[[ii]]/Pi_vec_new[ii]
  }

  for (ii in 1:Groups) {
    Cov_list_new[[ii]] <- (1-Gain(COUNT))*Pi_vec[ii]*(Cov_list[[ii]]+t(t(Mean_list[[ii]]))%*%t(Mean_list[[ii]]))+
      Gain(COUNT)*Tau_vec[ii]*t(t(X_vec))%*%t(X_vec)
    Cov_list_new[[ii]] <- Cov_list_new[[ii]]/Pi_vec_new[ii] - t(t(Mean_list_new[[ii]]))%*%t(Mean_list_new[[ii]])
  }

  Pi_vec <- Pi_vec_new
  Mean_list <- Mean_list_new
  Cov_list <- Cov_list_new

  print(Pi_vec)
  print(Mean_list)
  print(Cov_list)
}
