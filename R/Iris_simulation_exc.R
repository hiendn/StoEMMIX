library(MixSim)
data("iris", package = "datasets")
p <- ncol(iris) - 1
id <- as.integer(iris[, 5])
K <- max(id)
# estimate mixture parameters
Pi <- prop.table(tabulate(id))
Mu <- t(sapply(1:K, function(k){ colMeans(iris[id == k, -5]) }))
S <- sapply(1:K, function(k){ var(iris[id == k, -5]) })
dim(S) <- c(p, p, K)

NN <- 1000000
Groups <- 3
Results <- matrix(NA,20,5)
for (ii in 1:20) {
  Data <- simdataset(NN,Pi,Mu,S)$X
  Samp <- sample(1:Groups,NN,replace = T)
  msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
  MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=c(10,10)))
  Results[ii,1] <- MC$loglik
  print(Results)
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/1000,Groups,0.6,1,1000)
  Results[ii,2] <- Sto$`reg_log-likelihood`
  Results[ii,3] <- Sto$`pol_log-likelihood`
  print(Results)
  Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/2000,Groups,0.6,1,2000)
  Results[ii,4] <- Sto$`reg_log-likelihood`
  Results[ii,5] <- Sto$`pol_log-likelihood`
  print(Results)
}


NN <- 1000000
Groups <- 3
Results <- matrix(NA,20,7)
for (ii in 1:20) {
  rm(.Random.seed)
  Data <- simdataset(NN,Pi,Mu,S)$X
  Rand <- round(100000*runif(1))
  set.seed(Rand)
  # Samp <- Data[sample.int(NN,Groups),-1]
  KK <- kmeans(Data,centers=Groups,iter.max=1,nstart=1)
  Samp <- KK$centers
  set.seed(Rand)
  KK <- kmeans(Data,centers=Groups,iter.max=10,nstart=1)
  Results[ii,1] <- KK$tot.withinss
  print(Results)
  Sto <- softkmeans_pol(t(Data),t(Samp),
                        10*NN/1000,Groups,0.6,1,1000)
  Results[ii,2] <- sum(softkmeans_ss(t(Data),Sto$reg_means,softkmeans_clust(t(Data),Sto$reg_means)$allocations)$Within_SS)
  Results[ii,3] <- sum(softkmeans_ss(t(Data),Sto$pol_means,softkmeans_clust(t(Data),Sto$pol_means)$allocations)$Within_SS)
  print(Results)
  Sto <- softkmeans_pol(t(Data),t(Samp),
                        10*NN/2000,Groups,0.6,1,2000)
  Results[ii,4] <- sum(softkmeans_ss(t(Data),Sto$reg_means,softkmeans_clust(t(Data),Sto$reg_means)$allocations)$Within_SS)
  Results[ii,5] <- sum(softkmeans_ss(t(Data),Sto$pol_means,softkmeans_clust(t(Data),Sto$pol_means)$allocations)$Within_SS)
  print(Results)
  # set.seed(Rand)
  # MB <- MiniBatchKmeans(Data,clusters=Groups,batch_size = 1000,num_init = 1,max_iters = 10*NN/1000, early_stop_iter = 10*NN/1000,tol=1e-16,
  #                       CENTROIDS = KK$centers)
  # Results[ii,6] <- sum(softkmeans_ss(t(Data),t(MB$centroids),predict_KMeans(Data,MB$centroids)-1)$Within_SS)
  # set.seed(Rand)
  # MB <- MiniBatchKmeans(Data,clusters=Groups,batch_size = 2000,num_init = 1,max_iters = 10*NN/2000, early_stop_iter = 10*NN/2000,tol=1e-16,
  #                       CENTROIDS = KK$centers)
  # Results[ii,7] <- sum(softkmeans_ss(t(Data),t(MB$centroids),predict_KMeans(Data,MB$centroids)-1)$Within_SS)
  # print(Results)
}
