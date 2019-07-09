##################################################################
##                      Wreath Data Figure                      ##
##################################################################

# Load libraries
library(mclust)
library(colorspace)

# Load wreath data
data(wreath)

# Set a random seed (in case one is needed)
set.seed(20190129)

# Fit the Mclust model
MC <- Mclust(wreath, G = 14, modeNames = 'VVV')

# Open a plot device
pdf(file='./Wreath.pdf',width=6,height=6,paper='special') 

# Plot the data colored by the subpopulations
plot(wreath,col=rainbow_hcl(14)[MC$classification],
     xlab=expression(y[1]),ylab=expression(y[2]), lwd = 2)
grid()

# Plot subpopulation means
points(MC$parameters$mean[1,],MC$parameters$mean[2,],
       pch=4,cex=2,lwd=3)

# Close plot device
dev.off()

#################################################################
##                        Iris1 Results                        ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Iris1.rdata')

# Rename the columns of the data
colnames(Results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(Results)

# Open a plot device
pdf(file='./Iris1.pdf',width=12,height=4,paper='special') 

# Construct boxplot
boxplot(value ~ Var2, data = DF, lwd = 1, ylab = 'log-likelihood',range=0)

# Insert a grid
grid()

# Plot points over boxplot
stripchart(value ~ Var2, vertical = TRUE, data = DF, 
           method = "jitter", add = TRUE, pch = 20, 
           col = rainbow_hcl(dim(Results)[2]+1,alpha=0.5)[DF$Var1])

# Close plot device
dev.off()

#################################################################
##                        Iris1 Results timing                 ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Iris1timing.rdata')

# Rename the columns of the data
colnames(Timing) <- c('EM','N=n/10','N=n/5','N=n/10,T','N=n/5,T')

# Melt the data
DF <- melt(Timing)

# Open a plot device
pdf(file='./Iris1timing.pdf',width=12,height=4,paper='special') 

# Construct boxplot
boxplot(value ~ Var2, data = DF, lwd = 1, ylab = 'time (s)',range=0)

# Insert a grid
grid()

# Plot points over boxplot
stripchart(value ~ Var2, vertical = TRUE, data = DF, 
           method = "jitter", add = TRUE, pch = 20, 
           col = rainbow_hcl(dim(Timing)[2]+1,alpha=0.5)[DF$Var1])

# Close plot device
dev.off()

#################################################################
##                        Iris1 Results SE                     ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Iris1SE.rdata')

# Rename the columns of the data
colnames(SE_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(SE_results)

# Open a plot device
pdf(file='./Iris1SE.pdf',width=12,height=4,paper='special') 

# Construct boxplot
boxplot(value ~ Var2, data = DF, lwd = 1, ylab = 'SE',range=0)

# Insert a grid
grid()

# Plot points over boxplot
stripchart(value ~ Var2, vertical = TRUE, data = DF, 
           method = "jitter", add = TRUE, pch = 20, 
           col = rainbow_hcl(dim(SE_results)[2]+1,alpha=0.5)[DF$Var1])

# Close plot device
dev.off()

#################################################################
##                        Iris1 Results ARI                    ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Iris1ARI.rdata')

# Rename the columns of the data
colnames(ARI_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(ARI_results)

# Open a plot device
pdf(file='./Iris1ARI.pdf',width=12,height=4,paper='special') 

# Construct boxplot
boxplot(value ~ Var2, data = DF, lwd = 1, ylab = 'ARI',range=0)

# Insert a grid
grid()

# Plot points over boxplot
stripchart(value ~ Var2, vertical = TRUE, data = DF, 
           method = "jitter", add = TRUE, pch = 20, 
           col = rainbow_hcl(dim(ARI_results)[2]+1,alpha=0.5)[DF$Var1])

# Close plot device
dev.off()
