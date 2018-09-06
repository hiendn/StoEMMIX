hcTree <- hcVVV(data = Data)
cl <- hclass(hcTree, Groups)
msEst <- mstep(modelName = "VVV", data = Data, z = unmap(cl))

MC <- em('VVV', data=Data, parameters = msEst$parameters)
MC$parameters

REFIT <- em('VVV',data=Data,parameters =list(
  pro = Pi_vec,
  mean = matrix(unlist(Mean_list),Dim_vec[2],Groups),
  variance = Cov_list
))

