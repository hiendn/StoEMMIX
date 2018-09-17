library(MixSim)
library(mclust)
Groups <- 5
NN <- 10000
Q <- MixSim(MaxOmega = 0.5, K = Groups, p = 10,PiLow=1/(2*Groups))
A <- simdataset(n = NN, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.noise = 0)
Data <- A$X
# hcTree <- hcVVV(data = Data)
# cl <- hclass(hcTree, Groups)
Samp <- sample(1:Groups,NN,replace = T)
msEst <- mstep(modelName = "VVV", data = Data, z = unmap(Samp))
Dim_vec <- dim(Data)

MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=c(10,10)))

MC$parameters
attributes(MC)$info
MC$loglik
# 
# REFIT <- em('VVV',data=Data,parameters =list(
#   pro = Pi_vec,
#   mean = matrix(unlist(Mean_list),Dim_vec[2],Groups),
#   variance = Cov_list
# ))
# 