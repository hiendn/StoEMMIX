#################################################################
##                        Quin1 Results                        ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin1.rdata')

# Rename the columns of the data
colnames(Results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(Results)

# Open a plot device
pdf(file='./Quin1.pdf',width=12,height=4,paper='special') 

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
##                        Quin1 Results timing                 ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin1timing.rdata')

# Rename the columns of the data
colnames(Timing) <- c('EM','N=n/10','N=n/5','N=n/10,T','N=n/5,T')

# Melt the data
DF <- melt(Timing)

# Open a plot device
pdf(file='./Quin1timing.pdf',width=12,height=4,paper='special') 

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
##                        Quin1 Results SE                     ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin1SE.rdata')

# Rename the columns of the data
colnames(SE_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(SE_results)

# Open a plot device
pdf(file='./Quin1SE.pdf',width=12,height=4,paper='special') 

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
##                        Quin1 Results ARI                    ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin1ARI.rdata')

# Rename the columns of the data
colnames(ARI_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(ARI_results)

# Open a plot device
pdf(file='./Quin1ARI.pdf',width=12,height=4,paper='special') 

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
##                        Quin2 Results                        ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin2.rdata')

# Rename the columns of the data
colnames(Results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(Results)

# Open a plot device
pdf(file='./Quin2.pdf',width=12,height=4,paper='special') 

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
##                        Quin2 Results timing                 ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin2timing.rdata')

# Rename the columns of the data
colnames(Timing) <- c('EM','N=n/10','N=n/5','N=n/10,T','N=n/5,T')

# Melt the data
DF <- melt(Timing)

# Open a plot device
pdf(file='./Quin2timing.pdf',width=12,height=4,paper='special') 

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
##                        Quin2 Results SE                     ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin2SE.rdata')

# Rename the columns of the data
colnames(SE_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(SE_results)

# Open a plot device
pdf(file='./Quin2SE.pdf',width=12,height=4,paper='special') 

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
##                        Quin2 Results ARI                    ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Quin2ARI.rdata')

# Rename the columns of the data
colnames(ARI_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P','N=n/10,T','N=n/10,PT','N=n/5,T','N=n/5,PT')

# Melt the data
DF <- melt(ARI_results)

# Open a plot device
pdf(file='./Quin2ARI.pdf',width=12,height=4,paper='special') 

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

