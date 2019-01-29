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

