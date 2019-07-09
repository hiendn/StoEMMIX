#################################################################
##                        Exp1 Results                        ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp1.rdata')
Results <- Results[,-c(6:9)]

# Rename the columns of the data
colnames(Results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P')

# Melt the data
DF <- melt(Results)

# Open a plot device
pdf(file='./Exp1.pdf',width=12,height=4,paper='special') 

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
##                        Exp1 Results timing                 ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp1timing.rdata')
Timing <- Timing[,-c(4:5)]

# Rename the columns of the data
colnames(Timing) <- c('EM','N=n/10','N=n/5')

# Melt the data
DF <- melt(Timing)

# Open a plot device
pdf(file='./Exp1timing.pdf',width=12,height=4,paper='special') 

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
##                        Exp1 Results SE                     ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp1SE.rdata')
SE_results <- SE_results[,-c(6:9)]


# Rename the columns of the data
colnames(SE_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P')

# Melt the data
DF <- melt(SE_results)

# Open a plot device
pdf(file='./Exp1SE.pdf',width=12,height=4,paper='special') 

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
##                        Exp1 Results ARI                    ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp1ARI.rdata')
ARI_results <- ARI_results[,-c(6:9)]

# Rename the columns of the data
colnames(ARI_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P')

# Melt the data
DF <- melt(ARI_results)

# Open a plot device
pdf(file='./Exp1ARI.pdf',width=12,height=4,paper='special') 

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
##                        Exp2 Results                        ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp2.rdata')
Results <- Results[,-c(6:9)]

# Rename the columns of the data
colnames(Results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P')

# Melt the data
DF <- melt(Results)

# Open a plot device
pdf(file='./Exp2.pdf',width=12,height=4,paper='special') 

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
##                        Exp2 Results timing                 ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp2timing.rdata')
Timing <- Timing[,-c(4:5)]


# Rename the columns of the data
colnames(Timing) <- c('EM','N=n/10','N=n/5')

# Melt the data
DF <- melt(Timing)

# Open a plot device
pdf(file='./Exp2timing.pdf',width=12,height=4,paper='special') 

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
##                        Exp2 Results SE                     ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp2SE.rdata')
SE_results <- SE_results[,-c(6:9)]


# Rename the columns of the data
colnames(SE_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P')

# Melt the data
DF <- melt(SE_results)

# Open a plot device
pdf(file='./Exp2SE.pdf',width=12,height=4,paper='special') 

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
##                        Exp2 Results ARI                    ##
#################################################################

# Load libraries
library(reshape2)
library(colorspace)

# Load result data
load('Exp2ARI.rdata')
ARI_results <- ARI_results[,-c(6:9)]


# Rename the columns of the data
colnames(ARI_results) <- c('EM','N=n/10','N=n/10,P','N=n/5','N=n/5,P')

# Melt the data
DF <- melt(ARI_results)

# Open a plot device
pdf(file='./Exp2ARI.pdf',width=12,height=4,paper='special') 

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

