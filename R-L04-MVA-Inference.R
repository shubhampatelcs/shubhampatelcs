# This File contains codes for Performing Inference in MV

MVN.Simultaneous.CIs.f = function(X,method="Hotelling",conf.level = 0.95,ContrastMAT=NULL)
{
  # Contrast matrix should be given in such way that COLUMNS correspond to Contrasts
  # Computing the sample data information
  NOBS = nrow(X)
  NVAR = ncol(X)
  XBAR = apply(X,2,mean)
  SCOV = cov(X)
  
  # if no contrast matrix is provided, we compute all CIs for the marginal means
  # This is achieved by setting the Contrast Matrix to be the Identity
  if(is.null(ContrastMAT))
  {
    ContrastMAT=diag(NVAR)
  }else
  {
    # check if it is a single contrast/vector and convert it to a matrix
    if(is.null(dim(ContrastMAT)))
      ContrastMAT=matrix(ContrastMAT,ncol=1)
  }
  NCONTR = ncol(ContrastMAT)
  
  # Compute the covariance matrix of the contrasts
  CCOVMAT = t(ContrastMAT)%*%SCOV%*%ContrastMAT
  
  # Checking if method provided is one of the available ones
  METHODS = c("Hotelling","Chi-Sq","Bonferroni") 
  if(!method%in%METHODS)
  {
    stop(paste("Method Unknown:", method))
  }else
  {
    MRANK = match(method,METHODS)[1] # taking only first element if multiple methods provided
  }
  # Hotelling/Scheffe method
  CUTOFFS = vector(length = length(METHODS))
  CUTOFFS[1] = sqrt((NOBS-1)*NVAR/(NOBS-NVAR)* qf(1-conf.level,NVAR,NOBS-NVAR, lower.tail=FALSE))
  # Asymptotic/Chi-sq Method
  CUTOFFS[2] = sqrt(qchisq(1-conf.level,NVAR,lower.tail=FALSE))
  # Bonferroni Method
  CUTOFFS[3] = qt((1-conf.level)/(2*NCONTR),NOBS-1,lower.tail=FALSE)
  # Critical Value
  CRITVAL = CUTOFFS[MRANK]
  
  # Compute the Intervals
  SES = sqrt(1/NOBS * diag(CCOVMAT)) 
  PES = t(ContrastMAT)%*% XBAR
  LCBS =  PES - CRITVAL*SES
  UCBS =  PES + CRITVAL*SES
  CIs = data.frame(LCBS,UCBS)
  RES = list(method,conf.level,ContrastMAT,CIs)
  names(RES) = c("method","conf.level","Contrasts","CIs")
  return(RES)
}

###############################################################
###############################################################
# Multivariate Inference about mean vecotrs: Ex4.5
###############################################################
###############################################################
# Entering Data
X = array(dim=c(3,2))
X[1,] = c(6,9)
X[2,] = c(10,6)
X[3,] = c(8,3)

mu0 = c(9,5)
n = 3
# Computing numerical summaries
Xbar = colMeans(X)
S = var(X)

# T2 Hotteling Test Stat

IS = solve(S) 
T2 = n*t(Xbar-mu0)%*%IS%*%(Xbar-mu0)

# F critical value
FCV = qf(0.05,2,1, lower.tail=FALSE)

# Threshold 
Ta2 = (3-1)*2/(3-2)*FCV 

# S eigenvalues
EI = eigen(S)
EI$values

# Half lenght of eigen axes

sqrt(EI$values) * sqrt(Ta2/n)

###############################################################
###############################################################
# Example 4.7 Lectures (5.2 in JW)
###############################################################
###############################################################
# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
SData = read.xlsx(xlsxFile= '../Datasets/SweatData.xlsx')
head(SData,5) # Displays firs 5 rows of data
# COLUMNS: Both Numeric and Categorical. No mixing WITHIN 
# each column

###############################################################
# Numerical Summaries in R
###############################################################

XBar = apply(SData,2,mean) # Variable (column) means
# apply(SData,2,sd) # Variable (column) sds

S = cov(SData) # covariance Matrix

###############################################################
# Chi-sq QQPlot to Check Normality
###############################################################
# install.packages("heplots")
library(heplots)
D2 = cqplot(SData[,1:3], id.n=3) # id the 3 most extreme values

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot
mvn(MData, mvnTest = "royston")
# Henze-Zirkler test and construct perspective plot
mvn(MData, mvnTest = "hz")

###############################################################
###############################################################
# HT P-val: Hotelling's T^2 test
###############################################################
###############################################################
library(DescTools)
mu0 = c(6,50,10)
HotellingsT2Test(SData[,1:3], mu = mu0, test = "f")

# Note: To get the T^2_obs in the lecture we need to  
# multiply the value of the T.2 given from function as 
# T^2_obs = (n-1)*p*T.2/(n-p) = (n-1)*p/(n-p) * 5.7307 

###############################################################
# HT RR: Hotelling's T^2 test Manual computations
###############################################################
# Computing the T^2_obs manually
n=20
p = 3
IS = solve(S)
T2Obs = n*t(XBar-mu0)%*%IS%*%(XBar-mu0)
# F critical value
FCV = qf(0.05,p,n-p, lower.tail=FALSE)
# Threshold 
Ta2 = (n-1)*p/(n-p)*FCV 

###############################################################
# Simultaneous CIs: Using R packages
###############################################################
library(mvdalab)
MVcis(SData, level = 0.95)
# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(SData,method='Hotelling',conf.level = 0.95)

# Simultaneous CIs: Manual Computation
LowBound = XBar - sqrt(Ta2) * sqrt(1/n * diag(S))
UpperBound = XBar + sqrt(Ta2) * sqrt(1/n * diag(S))

###############################################################
###############################################################
# Chi-Sq analyses
###############################################################
###############################################################

# P-value, Need the T Obs computed earlier
pchisq(T2Obs,df = 3, lower.tail=FALSE)
# Critical value
Chi2a = qchisq(.05,3,lower.tail=FALSE)
# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(SData,method='Chi-Sq',conf.level = 0.95)

# Simultaneous CIs: Manual Computation
LowBound = XBar - sqrt(Chi2a) * sqrt(1/n * diag(S))
UpperBound = XBar + sqrt(Chi2a) * sqrt(1/n * diag(S))

###############################################################
###############################################################
# Bonferroni analyses
###############################################################
###############################################################
# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(SData,method='Bonferroni',conf.level = 0.95)
# Simultaneous CIs: Manual Computation
# Critical value
Ta = qt(.05/(2*3),n-1,lower.tail=FALSE)
LowBound = XBar - Ta * sqrt(1/n * diag(S))
UpperBound = XBar + Ta * sqrt(1/n * diag(S))

###############################################################
# Eigen values and Confidence region
###############################################################
# install.packages("MVQuickGraphs")
library(MVQuickGraphs)
# S eigenvalues
EI = eigen(S)
EI$values

# F critical value
FCV = qf(0.05,p,n-p, lower.tail=FALSE)
# Threshold 
Ta2 = (n-1)*p/(n-p)*FCV 

# Half lenght of eigen axes
sqrt(EI$values) * sqrt(Ta2/n)

# Confidence ellipses

confidenceEllipse(X.mean = XBar,
                  eig = EI,
                  n = n, p = p,
                  center=TRUE,
                  axes = FALSE,
                  alpha = 0.05
)
# Above function doesn't creat title and labels so we mus add them
title(main="A 95% Confidence Ellipse for mu = [mu1,mu2]",
      xlab="mu2 (Door Open)",ylab="mu1 (Door Closed)")






###############################################################
###############################################################
# Example 4.8 Lectures (5.3 in JW)
###############################################################
###############################################################
# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
MData = read.xlsx(xlsxFile= '../Datasets/MicrowaveData.xlsx')
head(MData,5) # Displays firs 5 rows of data
# Data MData in a Dataframe object, it allows for mixed type of 
# COLUMNS: Both Numeric and Categorical. No mixing WITHIN 
# each column

# Remove 1st column
MData = MData[, -1]

###############################################################
# Chi-sq QQPlot to Check Normality
###############################################################
library(heplots)
cqplot(MData, id.n=3) # id the 3 most extreme values

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test
mvn(MData, mvnTest = "royston")
# Henze-Zirkler test 
mvn(MData, mvnTest = "hz")

###############################################################
# Tranform the data
###############################################################

MData[,1] = (MData[,1])^(.25)
MData[,2] = (MData[,2])^(.25)

###############################################################
# Numerical Summaries in R
###############################################################

XBar = apply(MData,2,mean) # Variable (column) means
# apply(SData,2,sd) # Variable (column) sds

S = cov(MData) # covariance Matrix

###############################################################
# Chi-sq QQPlot to Check Normality
###############################################################
# install.packages("heplots")
library(heplots)
cqplot(MData, id.n=3) # id the 3 most extreme values

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test
mvn(MData, mvnTest = "royston")
# Henze-Zirkler test 
mvn(MData, mvnTest = "hz")

###############################################################
###############################################################
# HT P-val: Hotelling's T^2 test
###############################################################
###############################################################
library(DescTools)
mu0 = c(.5,.5)
HotellingsT2Test(MData, mu = mu0, test = "f")

# Note: To get the T^2_obs in the lecture we need to  
# multiply the value of the T.2 given from function as 
# T^2_obs = (n-1)*p*T.2/(n-p) = (n-1)*p/(n-p) * 5.7307 

###############################################################
# HT RR: Hotelling's T^2 test Manual computations
###############################################################

# Computing the T^2_obs manually
n=42
p = 2
IS = solve(S)
T2Obs = n*t(XBar-mu0)%*%IS%*%(XBar-mu0)
# F critical value
FCV = qf(0.05,p,n-p, lower.tail=FALSE)
# Threshold 
Ta2 = (n-1)*p/(n-p)*FCV 

###############################################################
# Simultaneous CIs: Using R packages
###############################################################
library(mvdalab)
MVcis(MData, level = 0.95)
# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(MData,method='Hotelling',conf.level = 0.95)

# Simultaneous CIs: Manual Computation
LowBound = XBar - sqrt(Ta2) * sqrt(1/n * diag(S))
UpperBound = XBar + sqrt(Ta2) * sqrt(1/n * diag(S))

###############################################################
# Confidence interval for linear combinations: 
###############################################################
# Contrast Mat with the linear combinations as columns
ContMat = matrix(nrow=p,ncol=3) 
ContMat[,1] = c(1,1) # X1+X2
ContMat[,2] = c(-1,1) # X2-X1
ContMat[,3] = c(.5,.5) # (X1+X2)/2

# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(MData,method='Hotelling',conf.level = 0.95,ContrastMAT = ContMat)

# Manual Computations 
AXBar = t(ContMat)%*%XBar 
ASA = t(ContMat)%*%S%*%ContMat

LowBound = AXBar - sqrt(Ta2) * sqrt(1/n * diag(ASA))
UpperBound = AXBar + sqrt(Ta2) * sqrt(1/n * diag(ASA))

###############################################################
###############################################################
# Chi-Sq analyses
###############################################################
###############################################################

# P-value, Need the T Obs computed earlier
pchisq(T2Obs,df = 2, lower.tail=FALSE)
# Critical value
Chi2a = qchisq(.05,2,lower.tail=FALSE)
# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(MData,method='Chi-Sq',conf.level = 0.95)

# Simultaneous CIs: Manual Computation
LowBound = XBar - sqrt(Chi2a) * sqrt(1/n * diag(S))
UpperBound = XBar + sqrt(Chi2a) * sqrt(1/n * diag(S))

###############################################################
# Confidence interval for linear combinations: 
###############################################################
# Contrast Mat with the linear combinations as columns
ContMat = matrix(nrow=p,ncol=3) 
ContMat[,1] = c(1,1) # X1+X2
ContMat[,2] = c(-1,1) # X2-X1
ContMat[,3] = c(.5,.5) # (X1+X2)/2

# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(MData,method='Chi-Sq',conf.level = 0.95,ContrastMAT = ContMat)

# Manual Computations 
AXBar = t(ContMat)%*%XBar 
ASA = t(ContMat)%*%S%*%ContMat

LowBound = AXBar - sqrt(Chi2a) * sqrt(1/n * diag(ASA))
UpperBound = AXBar + sqrt(Chi2a) * sqrt(1/n * diag(ASA))

###############################################################
###############################################################
# Bonferroni analyses
###############################################################
###############################################################
# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(MData,method='Bonferroni',conf.level = 0.95)
# Simultaneous CIs: Manual Computation
# Critical value
Ta = qt(.05/(2*2),n-1,lower.tail=FALSE)
LowBound = XBar - Ta * sqrt(1/n * diag(S))
UpperBound = XBar + Ta * sqrt(1/n * diag(S))

###############################################################
# Confidence interval for linear combinations: 
###############################################################
# Contrast Mat with the linear combinations as columns
ContMat = matrix(nrow=p,ncol=3) 
ContMat[,1] = c(1,1) # X1+X2
ContMat[,2] = c(-1,1) # X2-X1
ContMat[,3] = c(.5,.5) # (X1+X2)/2

# Simultaneous CIs: Using instructor function
MVN.Simultaneous.CIs.f(MData,method='Bonferroni',conf.level = 0.95,ContrastMAT = ContMat)

# Manual Computations 
AXBar = t(ContMat)%*%XBar 
ASA = t(ContMat)%*%S%*%ContMat

# Need to adjust Critical Value cause more CI's 
Tac = qt(.05/(2*3),n-1,lower.tail=FALSE)

LowBound = AXBar - Tac * sqrt(1/n * diag(ASA))
UpperBound = AXBar + Tac * sqrt(1/n * diag(ASA))

###############################################################
###############################################################
# eigenvalues analyses
###############################################################
###############################################################

###############################################################
# Eigen values and Confidence region
###############################################################
# install.packages("MVQuickGraphs")
library(MVQuickGraphs)
# S eigenvalues
EI = eigen(S)
EI$values

# F critical value
FCV = qf(0.05,p,n-p, lower.tail=FALSE)
# Threshold 
Ta2 = (n-1)*p/(n-p)*FCV 

# Half lenght of eigen axes
sqrt(EI$values) * sqrt(Ta2/n)

# Confidence ellipses

confidenceEllipse(X.mean = XBar,
                  eig = EI,
                  n = n, p = p,
                  center=TRUE,
                  axes = FALSE,
                  alpha = 0.05
)
# Above function doesn't creat title and labels so we mus add them
title(main="A 95% Confidence Ellipse for mu = [mu1,mu2]",
      xlab="mu2 (Door Open)",ylab="mu1 (Door Closed)")












