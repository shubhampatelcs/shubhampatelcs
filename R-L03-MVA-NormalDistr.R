##############################################################################
##############################################################################
#
# This file contains the codes for transforming different types of variables 
# to Obtain a normally distributed data
#
##############################################################################
##############################################################################

# install.packages("EnvStats") # Uncommnet this line the first time you run the program
library(EnvStats) # Needed for the Boxcox transformation

data("USJudgeRatings")
attach(USJudgeRatings)
head(USJudgeRatings)

##############################################################################
# Skewed Right Data
##############################################################################

x = CONT

dev.off() # resets all graphical parameters 
# Histograms of Original and Transformed Data
par(mfrow=c(2,3))
hist(x,breaks=8,main="Original Data")
hist(sqrt(x),breaks=8,main="Square Root Transformation")
# If the min value the data is less than 1, the logx/(1/x) transformation will not
# work well. To avoid issues we add 1 to move it to be bigger than 1. 
if(min(x)<1)
{
  x1 = x + 1
}
hist(log(x1),breaks=8,main="Log Transformation")
hist(1/x1,breaks=8,main="Inverse Transformation")

# Box-Cox Transformation
# Estimate value of lambda
BCL = boxcox(x, optimize = TRUE)
lambda = BCL$lambda # best guess of lambda
x2 = (x^lambda-1)/lambda # transformed data

hist(x2,breaks=8,main="Box-Cox Transformation")

dev.off() # resets all graphical parameters 
# QQPlots of Original and Transformed Data
par(mfrow=c(2,3))
qqnorm(x,main="Original Data")
qqline(x,lwd=1.5,col="red")
qqnorm(sqrt(x),main="Square Root Transformation")
qqline(sqrt(x),lwd=1.5,col="red")
qqnorm(log(x1),main="Log Transformation")
qqline(log(x1),lwd=1.5,col="red")
qqnorm(1/x1,main="Inverse Transformation")
qqline(1/x1,lwd=1.5,col="red")
qqnorm(x2,main="Box-Cox Transformation")
qqline(x2,lwd=1.5,col="red")

##############################################################################
# Skewed Left Data
##############################################################################
x.org = PHYS

# Data are skewed left so we need to "Mirror Transform" the data first
x = max(x.org+1)-x.org 

dev.off() # resets all graphical parameters 
# Histograms of Original and Transformed Data
par(mfrow=c(2,3))
hist(x.org,breaks=8,main="Original Data")
hist(x,breaks=8,main="Mirror Transformed Data")
hist(sqrt(x),breaks=8,main="Square Root Transformation")
# If the min value the data is less than 1, the logx/(1/x) transformation will not
# work well. To avoid issues we add 1 to move it to be bigger than 1. 
if(min(x)<1)
{
  x = x + 1
}
hist(log(x),breaks=8,main="Log Transformation")
hist(1/x,breaks=8,main="Inverse Transformation")

# Box-Cox Transformation
# Estimate value of lambda
BCL = boxcox(x, optimize = TRUE)
lambda = BCL$lambda # best guess of lambda
x1 = (x^lambda-1)/lambda # transformed data

hist(x1,breaks=8,main="Box-Cox Transformation")

dev.off() # resets all graphical parameters 
# QQPlots of Original and Transformed Data
par(mfrow=c(2,3))
qqnorm(x.org,main="Original Data")
qqline(x.org,lwd=1.5,col="red")
qqnorm(x,main="Mirror Transformed Data")
qqline(x,lwd=1.5,col="red")
qqnorm(sqrt(x),main="Square Root Transformation")
qqline(sqrt(x),lwd=1.5,col="red")
# If the min value the data is less than 1, the logx/(1/x) transformation will not
# work well. To avoid issues we add 1 to move it to be bigger than 1. 
if(min(x)<1)
{
  x = x + 1
}
qqnorm(log(x),main="Log Transformation")
qqline(log(x),lwd=1.5,col="red")
qqnorm(1/x,main="Inverse Transformation")
qqline(1/x,lwd=1.5,col="red")
qqnorm(x1,main="Box-Cox Transformation")
qqline(x1,lwd=1.5,col="red")

##############################################################################
# Very Skewed Right Data
##############################################################################
set.seed(113)
x = rgamma(200, shape = 1, scale = 1)
x = x[x>.05]

dev.off() # resets all graphical parameters 
# Histograms of Original and Transformed Data
par(mfrow=c(2,3))
hist(x,breaks=8,main="Original Data")
hist(sqrt(x),breaks=8,main="Square Root Transformation")

# If the min value the data is less than 1, the logx/(1/x) transformation will not
# work well. To avoid issues we add 1 to move it to be bigger than 1. 
if(min(x)<1)
{
  x = x + 1
}

hist(log(x),breaks=8,main="Log Transformation")
hist(1/x,breaks=8,main="Inverse Transformation")

# Box-Cox Transformation
# Estimate value of lambda
BCL = boxcox(x, optimize = TRUE)
lambda = BCL$lambda # best guess of lambda
x1 = (x^lambda-1)/lambda # transformed data

hist(x1,breaks=8,main="Box-Cox Transformation")

dev.off() # resets all graphical parameters 
# QQPlots of Original and Transformed Data
par(mfrow=c(2,3))
qqnorm(x,main="Original Data")
qqline(x,lwd=1.5,col="red")
qqnorm(sqrt(x),main="Square Root Transformation")
qqline(sqrt(x),lwd=1.5,col="red")

# If the min value the data is less than 1, the logx/(1/x) transformation will not
# work well. To avoid issues we add 1 to move it to be bigger than 1. 
if(min(x)<1)
{
  x = x + 1
}

qqnorm(log(x),main="Log Transformation")
qqline(log(x),lwd=1.5,col="red")
qqnorm(1/x,main="Inverse Transformation")
qqline(1/x,lwd=1.5,col="red")
qqnorm(x1,main="Box-Cox Transformation")
qqline(x1,lwd=1.5,col="red")


###############################################################
###############################################################
# Checking Multivariate Normality
###############################################################
###############################################################

###############################################################
# Reading Data in R
###############################################################
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')
# Read data
SData = read.xlsx(xlsxFile= '../Datasets/StiffnessData.xlsx')
head(SData,5) # Displays firs 5 rows of data

###############################################################
###############################################################
# Univariate Analysis
###############################################################
###############################################################

###############################################################
# Plot Histograms
###############################################################
# Long way 

dev.off() # resets all graphical parameters 

par(mfrow = c(2,2))

hist(SData[, 1], probability = T,
     xlab = paste("x", 1, sep=""),
     main = paste0("Histogram of x", 1))

hist(SData[, 2], probability = T,
     xlab = paste("x", 2, sep=""),
     main = paste0("Histogram of x", 1))

hist(SData[, 3], probability = T,
     xlab = paste("x", 3, sep=""),
     main = paste0("Histogram of x", 1))

hist(SData[, 4], probability = T,
     xlab = paste("x", 4, sep=""),
     main = paste0("Histogram of x", 1))

###############################################################
# Plot Histograms: Looping
###############################################################

dev.off() # resets all graphical parameters 
par(mfrow = c(2,2))

for(i in 1:4){
  hist(SData[, i], probability = T,
       xlab = paste("x", i, sep=""),
       main = paste0("Histogram of x", i))
}

###############################################################
# Plot Histograms: using apply function
###############################################################

apply(SData[,1:4],2,hist) # apply given functon (hist) to each col (2)


###############################################################
# Plot QQPlots: Looping
###############################################################

dev.off() # resets all graphical parameters 
par(mfrow = c(2,2))

for(i in 1:4){
  qqnorm(SData[, i],
         main = paste0("Q-Q plot of x", i), pch=19, cex=1.5)
qqline(SData[, i],col = "#006CD1")
  }

###############################################################
# Shapiro-Wilk
###############################################################

# Shapiro-Wilks test for each column using the apply function
# See ?apply for more details
apply(SData[,1:4], 2, shapiro.test)


###############################################################
###############################################################
# Multivariate Analysis
###############################################################
###############################################################

###############################################################
# Scatter matrix
###############################################################

# install.packages("psych")
library(psych) 

dev.off() # resets all graphical parameters 
pairs.panels(SData[,1:4], main = "Advanced Scatter Matrix for Paper Quality Dataset",
             smooth=FALSE,
             ellipses=FALSE) 

###############################################################
# Scatter plots with ellipses to identify outliers
###############################################################

library(car)

dev.off() # resets all graphical parameters 
par(mfrow = c(2,3))

for(i in 1:3){
  for(j in (i+1):4){
  x1 <- SData[, i]
x2 <- SData[, j]
dataEllipse(x1, x2,
            xlim = c(700, 3500), ylim = c(700, 3000),
            pch=19, col = c("#994F00", "#006CD1"), lty=2,
            ellipse.label=c(0.5, 0.95), levels = c(0.5, 0.95),
            fill=TRUE, fill.alpha=0.1,
            xlab = paste0("X",i),
            ylab = paste0("X",j)
            )
  }
}

###############################################################
# Chi-sq QQPlot to identify outliers: Version 1
###############################################################
# install.packages("heplots")
library(heplots)

dev.off() # resets all graphical parameters 
cqplot(SData[,1:4], id.n=3) # id the 3 most extreme values
dev.off()

###############################################################
# Chi-sq QQPlot to identify outliers: Version 2
###############################################################
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot

dev.off() # resets all graphical parameters 
mvn(SData[, 1:4], mvnTest = "royston", multivariatePlot = "qq")

###############################################################
# Estimated bivariate joint pdf
###############################################################
dev.off() # resets all graphical parameters 
# Henze-Zirkler test and construct perspective plot
mvn(SData[, 1:2], mvnTest = "hz", multivariatePlot = "persp")

###############################################################
# Estimated bivariate contours of joint pdf
###############################################################
dev.off() # resets all graphical parameters 
# Henze-Zirkler test and contour plot
mvn(SData[, 1:2], mvnTest = "hz", multivariatePlot = "contour")







