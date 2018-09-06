library(MixSim)
library(mclust)
Groups <- 5
NN <- 1000000
Q <- MixSim(MaxOmega = 0.5, K = Groups, p = 3,PiLow=1/(2*Groups))
A <- simdataset(n = NN, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, n.noise = 0)
Data <- A$X
# hcTree <- hcVVV(data = Data)
# cl <- hclass(hcTree, Groups)
msEst <- mstep(modelName = "VVV", data = Data, z = unmap(sample(1:Groups,NN,replace = T)))
Dim_vec <- dim(Data)

MC <- em('VVV', data=Data, parameters = msEst$parameters)

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
