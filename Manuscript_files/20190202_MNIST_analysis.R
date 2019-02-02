#################################################################
##                            MNIST                            ##
#################################################################

# Load source files for minibatch EM algorithms
source('https://raw.githubusercontent.com/hiendn/StoEMMIX/master/Manuscript_files/20190128_main_functions.R')

# Set memory limit
Sys.setenv('R_MAX_VSIZE'=10000000000000)

# Set number of PCs
dPC <- 20

# Declare number of epochs
Epoch <- 10

# Set Data to be PCA of dimensions dPC
Data <- PCA$scores[,1:dPC]

# Set a random seed
set.seed(20190202)

# Conduct a K-means initialization
KM <- kmeans(Data,centers = 10,iter.max = 100, nstart = 10)

# Extract K-means clusters
id <- KM$cluster
# Get maximum number of clusters
K <- max(id)

# Estimate mixture parameters from K-means initialization
Pi <- prop.table(tabulate(id))
Mu <- t(sapply(1:K, function(k){ colMeans(Data[id == k,]) }))
S <- sapply(1:K, function(k){ var(Data[id == k,]) })
dim(S) <- c(dim(Data)[2], dim(Data)[2], K)

Sto <- stoEMMIX_poltrunc(t(Data), Pi, t(Mu), S,
                         10*dim(Data)[1]/7000,K,0.6,1-10^-10,7000,
                         1000,1000,1000)

