library(mclust)
library(colorspace)

data(wreath)
MC2 <- Mclust(wreath, G = 14, modeNames = 'VVV')

SIM <- simVVV(MC$parameters,10000)

plot(SIM[,c(2,3)],
     col=rainbow_hcl(14)[SIM[,1]],
     xlab='',ylab='',axes = F)


plot(wreath,col=rainbow_hcl(14)[MC$classification])

par(mar=c(1,1,1,1))
plot(MC,what='density',xlab='',ylab='',axes = F,main='')
points(SIM[,c(2,3)],
       col=rainbow_hcl(14)[SIM[,1]])

NN <- 1000000
Groups <- 14
Results <- matrix(NA,20,5)
for (ii in 1:20) {
  Data <- simVVV(MC2$parameters,NN)
  Samp <- sample(1:Groups,NN,replace = T)
  msEst <- mstep(modelName = "VVV", data = Data[,-1], z = unmap(Samp))
  MC <- em('VVV', data=Data[,-1], parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=c(10,10)))
  Results[ii,1] <- MC$loglik
  print(Results)
  Sto <- stoEMMIX_pol(t(Data[,-1]), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/1000,Groups,0.6,1,1000)
  Results[ii,2] <- Sto$`reg_log-likelihood`
  Results[ii,3] <- Sto$`pol_log-likelihood`
  print(Results)
  Sto <- stoEMMIX_pol(t(Data[,-1]), msEst$parameters$pro, msEst$parameters$mean,
                      msEst$parameters$variance$sigma,
                      10*NN/2000,Groups,0.6,1,2000)
  Results[ii,4] <- Sto$`reg_log-likelihood`
  Results[ii,5] <- Sto$`pol_log-likelihood`
  print(Results)
}

NN <- 1000000
Groups <- 14
Results <- matrix(NA,20,7)
for (ii in 1:20) {
  rm(.Random.seed)
  Data <- simVVV(MC2$parameters,NN)
  Rand <- round(100000*runif(1))
  set.seed(Rand)
  # Samp <- Data[sample.int(NN,Groups),-1]
  KK <- kmeans(Data[,-1],centers=Groups,iter.max=1,nstart=1)
  Samp <- KK$centers
  set.seed(Rand)
  KK <- kmeans(Data[,-1],centers=Groups,iter.max=10,nstart=1)
  Results[ii,1] <- KK$tot.withinss
  print(Results)
  Sto <- softkmeans_pol(t(Data[,-1]),t(Samp),
                      10*NN/1000,Groups,0.6,1,1000)
  Results[ii,2] <- sum(softkmeans_ss(t(Data[,-1]),Sto$reg_means,softkmeans_clust(t(Data[,-1]),Sto$reg_means)$allocations)$Within_SS)
  Results[ii,3] <- sum(softkmeans_ss(t(Data[,-1]),Sto$pol_means,softkmeans_clust(t(Data[,-1]),Sto$pol_means)$allocations)$Within_SS)
  print(Results)
  Sto <- softkmeans_pol(t(Data[,-1]),t(Samp),
                        10*NN/2000,Groups,0.6,1,2000)
  Results[ii,4] <- sum(softkmeans_ss(t(Data[,-1]),Sto$reg_means,softkmeans_clust(t(Data[,-1]),Sto$reg_means)$allocations)$Within_SS)
  Results[ii,5] <- sum(softkmeans_ss(t(Data[,-1]),Sto$pol_means,softkmeans_clust(t(Data[,-1]),Sto$pol_means)$allocations)$Within_SS)
  print(Results)
  # set.seed(Rand)
  # MB <- MiniBatchKmeans(Data[,-1],clusters=Groups,batch_size = 1000,num_init = 1,max_iters = 10*NN/1000, early_stop_iter = 10*NN/1000,tol=1e-16,
  #                       CENTROIDS = KK$centers)
  # Results[ii,6] <- sum(softkmeans_ss(t(Data[,-1]),t(MB$centroids),predict_KMeans(Data[,-1],MB$centroids)-1)$Within_SS)
  # set.seed(Rand)
  # MB <- MiniBatchKmeans(Data[,-1],clusters=Groups,batch_size = 2000,num_init = 1,max_iters = 10*NN/2000, early_stop_iter = 10*NN/2000,tol=1e-16,
  #                       CENTROIDS = KK$centers)
  # Results[ii,7] <- sum(softkmeans_ss(t(Data[,-1]),t(MB$centroids),predict_KMeans(Data[,-1],MB$centroids)-1)$Within_SS)
  # print(Results)
  }

library(reshape2)
colnames(Results) <- c('k-means','soft1r','soft1p','soft2r','soft2p','mini1','mini2')
DF <- melt(Results)
boxplot(value ~ Var2, data = DF, lwd = 2, ylab = 'TSS')
stripchart(value ~ Var2, vertical = TRUE, data = DF, 
           method = "jitter", add = TRUE, pch = 20, col = rainbow_hcl(20)[DF$Var1])
