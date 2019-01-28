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

MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=1))
MC$loglik
MC <- em('VVV', data=Data, parameters = msEst$parameters, control = emControl(eps=0,tol=0,itmax=10))
MC$loglik
Sto <- stoEMMIX_pol(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                    msEst$parameters$variance$sigma,
                    1000,5,0.6,1,100)
Sto$`reg_log-likelihood`
Sto$`pol_log-likelihood`
Stot <- stoEMMIX_poltrunc(t(Data), msEst$parameters$pro, msEst$parameters$mean,
                         msEst$parameters$variance$sigma,
                         1000,5,0.6,1,100,1000,1000,1000)
Stot$`reg_log-likelihood`
Stot$`pol_log-likelihood`
# 
# REFIT <- em('VVV',data=Data,parameters =list(
#   pro = Pi_vec,
#   mean = matrix(unlist(Mean_list),Dim_vec[2],Groups),
#   variance = Cov_list
# ))
# 