#################################################################
##                        Wreath1 Results                        ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath1R1.rdata')

# Rename the columns of the data
colnames(Results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(Results)

# Open a plot device
pdf(file='./Wreath1.pdf',width=12,height=4,paper='special') 

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
##                        Wreath1 Results timing                 ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath1timing.rdata')

# Rename the columns of the data
colnames(Timing) <- c('EM','N=n/10','N=n/5','N=n/10,T','N=n/5,T')

# Melt the data
DF <- melt(Timing)

# Open a plot device
pdf(file='./Wreath1timing.pdf',width=12,height=4,paper='special') 

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
##                        Wreath1 Results SE                     ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath1SE.rdata')

# Rename the columns of the data
colnames(SE_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(SE_results)

# Open a plot device
pdf(file='./Wreath1SE.pdf',width=12,height=4,paper='special') 

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
##                        Wreath1 Results ARI                    ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath1ARI.rdata')

# Rename the columns of the data
colnames(ARI_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(ARI_results)

# Open a plot device
pdf(file='./Wreath1ARI.pdf',width=12,height=4,paper='special') 

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

#################################################################
##                        Wreath2 Results                        ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath2R1.rdata')

# Rename the columns of the data
colnames(Results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(Results)

# Open a plot device
pdf(file='./Wreath2.pdf',width=12,height=4,paper='special') 

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
##                        Wreath2 Results timing                 ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath2timing.rdata')

# Rename the columns of the data
colnames(Timing) <- c('EM','N=n/10','N=n/5','N=n/10,T','N=n/5,T')

# Melt the data
DF <- melt(Timing)

# Open a plot device
pdf(file='./Wreath2timing.pdf',width=12,height=4,paper='special') 

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
##                        Wreath2 Results SE                     ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath2SE.rdata')

# Rename the columns of the data
colnames(SE_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(SE_results)

# Open a plot device
pdf(file='./Wreath2SE.pdf',width=12,height=4,paper='special') 

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
##                        Wreath2 Results ARI                    ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Wreath2ARI.rdata')

# Rename the columns of the data
colnames(ARI_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(ARI_results)

# Open a plot device
pdf(file='./Wreath2ARI.pdf',width=12,height=4,paper='special') 

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

