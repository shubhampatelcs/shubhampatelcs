
MVN.HT.CIs.f = function(X,conf.level = 0.95,alpha=.05,mu0=rep(0,ncol(X)),ContrastMAT=NULL,SigDig=3)
{
  # Contrast matrix should be given in such way that COLUMNS correspond to Contrasts
  # Computing the sample data information
  NOBS = nrow(X)
  NVAR = ncol(X)
  XBAR = apply(X,2,mean)
  SCOV = cov(X)
  
  # HT info:
  ISCOV = solve(SCOV)
  TSQOBS = NOBS*t(XBAR-mu0)%*%ISCOV%*%(XBAR-mu0)
  
  HTMETHODS = c("Hotelling","Chi-Sq") 
  HTCUTOFFS = vector(length = length(HTMETHODS))
  # Hotelling/Scheffe method
  HTCUTOFFS[1] = qf(alpha,NVAR,NOBS-NVAR, lower.tail=FALSE)
  # Asymptotic/Chi-sq Method
  HTCUTOFFS[2] = qchisq(alpha,NVAR,lower.tail=FALSE)
  
  HTPVALS = vector(length = length(HTMETHODS)) 
  HTPVALS[1] = pf((NOBS-NVAR)/((NOBS-1)*NVAR)*TSQOBS,NVAR,NOBS-NVAR, lower.tail=FALSE)
  HTPVALS[2] = pchisq(TSQOBS,df = NVAR, lower.tail=FALSE)
  
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
  NMETH =  length(METHODS)
  # Hotelling/Scheffe method
  # Critical Value
  CUTOFFS = vector(length = length(METHODS))
  CUTOFFS[1] = sqrt((NOBS-1)*NVAR/(NOBS-NVAR)* qf(1-conf.level,NVAR,NOBS-NVAR, lower.tail=FALSE))
  # Asymptotic/Chi-sq Method
  CUTOFFS[2] = sqrt(qchisq(1-conf.level,NVAR,lower.tail=FALSE))
  # Bonferroni Method
  CUTOFFS[3] = qt((1-conf.level)/(2*NCONTR),NOBS-1,lower.tail=FALSE)
  
  # Compute the Intervals
  CIs = vector(mode="list", length = NMETH) 
  # lapply(1:NMETH, matrix, data= NA, nrow=NCONTR, ncol=2)
  
  SES = sqrt(1/NOBS * diag(CCOVMAT)) 
  PES = t(ContrastMAT)%*% XBAR
  for(m in 1:NMETH)
  {
    CRITVAL = CUTOFFS[m]
    LCBS =  PES - CRITVAL*SES
    UCBS =  PES + CRITVAL*SES
    CIs[[m]] = round(data.frame(LCBS,UCBS),SigDig)
  }
  CUTOFFS[1] = qf(1-conf.level,NVAR,NOBS-NVAR, lower.tail=FALSE)
  CUTOFFS[2] = qchisq(1-conf.level,NVAR,lower.tail=FALSE)
  names(CIs) = METHODS
  RES = list(NOBS,round(XBAR,SigDig),round(SCOV,SigDig),round(ISCOV,SigDig),mu0,
             round(TSQOBS,3),HTPVALS,alpha,round(HTCUTOFFS,3),conf.level,ContrastMAT,
             round(PES,SigDig),round(diag(CCOVMAT),SigDig),
             round(CUTOFFS,3),CIs)
  names(RES) = c( "N","XBar","S","IS","Mu0","Tobs","PVals","Alpha","HT CritVal","CI conf.level",
                  "Contrasts","CMeans","CVars","CI Crit. Vals","CIs")
  return(RES)
}
# MVN.HT.CIs.f(CData,conf.level = 0.95,alpha=.05,mu0=mu0,ContrastMAT=NULL)

MVN.HT.CIs.NumSum.f = function(NOBS,XBAR,SCOV,conf.level = 0.95,alpha=.05,mu0=rep(0,ncol(X)),ContrastMAT=NULL,SigDig=3)
{
  # Contrast matrix should be given in such way that COLUMNS correspond to Contrasts
  # Computing the sample data information
  NVAR = ncol(SCOV)

  # HT info:
  ISCOV = solve(SCOV)
  TSQOBS = NOBS*t(XBAR-mu0)%*%ISCOV%*%(XBAR-mu0)
  
  HTMETHODS = c("Hotelling","Chi-Sq") 
  HTCUTOFFS = vector(length = length(HTMETHODS))
  # Hotelling/Scheffe method
  HTCUTOFFS[1] = qf(alpha,NVAR,NOBS-NVAR, lower.tail=FALSE)
  # Asymptotic/Chi-sq Method
  HTCUTOFFS[2] = qchisq(alpha,NVAR,lower.tail=FALSE)
  
  HTPVALS = vector(length = length(HTMETHODS)) 
  HTPVALS[1] = pf((NOBS-NVAR)/((NOBS-1)*NVAR)*TSQOBS,NVAR,NOBS-NVAR, lower.tail=FALSE)
  HTPVALS[2] = pchisq(TSQOBS,df = NVAR, lower.tail=FALSE)
  
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
  NMETH =  length(METHODS)
  # Hotelling/Scheffe method
  # Critical Value
  CUTOFFS = vector(length = length(METHODS))
  CUTOFFS[1] = sqrt((NOBS-1)*NVAR/(NOBS-NVAR)* qf(1-conf.level,NVAR,NOBS-NVAR, lower.tail=FALSE))
  # Asymptotic/Chi-sq Method
  CUTOFFS[2] = sqrt(qchisq(1-conf.level,NVAR,lower.tail=FALSE))
  # Bonferroni Method
  CUTOFFS[3] = qt((1-conf.level)/(2*NCONTR),NOBS-1,lower.tail=FALSE)
  
  # Compute the Intervals
  CIs = vector(mode="list", length = NMETH) 
  # lapply(1:NMETH, matrix, data= NA, nrow=NCONTR, ncol=2)
  
  SES = sqrt(1/NOBS * diag(CCOVMAT)) 
  PES = t(ContrastMAT)%*% XBAR
  for(m in 1:NMETH)
  {
    CRITVAL = CUTOFFS[m]
    LCBS =  PES - CRITVAL*SES
    UCBS =  PES + CRITVAL*SES
    CIs[[m]] = round(data.frame(LCBS,UCBS),SigDig)
  }
  CUTOFFS[1] = qf(1-conf.level,NVAR,NOBS-NVAR, lower.tail=FALSE)
  CUTOFFS[2] = qchisq(1-conf.level,NVAR,lower.tail=FALSE)
  names(CIs) = METHODS
  RES = list(NOBS,round(XBAR,SigDig),round(SCOV,SigDig),round(ISCOV,SigDig),mu0,
             round(TSQOBS,3),HTPVALS,alpha,round(HTCUTOFFS,3),conf.level,ContrastMAT,
             round(PES,SigDig),round(diag(CCOVMAT),SigDig),
             round(CUTOFFS,3),CIs)
  names(RES) = c( "N","XBar","S","IS","Mu0","Tobs","PVals","Alpha","HT CritVal","CI conf.level",
                  "Contrasts","CMeans","CVars","CI Crit. Vals","CIs")
  return(RES)
}

MVN2Sample.HT.CIs.f = function(X,Y,conf.level = 0.95,alpha=.05,mu0=rep(0,ncol(X)),
                               ContrastMAT=NULL,SigDig=3,var.eq=TRUE)
{
  N1 = nrow(X)
  X1Bar = colMeans(X)
  S1 = cov(X)
  N2 = nrow(Y)
  X2Bar = colMeans(Y)
  S2 = cov(Y)
  
  OUTHTCI = MVN2Sample.HT.CIs.NumSum.f(NOBS1=N1,XBAR1=X1Bar,SCOV1=S1,NOBS2=N2,XBAR2=X2Bar,SCOV2=S2,
                                       conf.level = conf.level,
                                       alpha=alpha,mu0=mu0,ContrastMAT=ContrastMAT,SigDig=SigDig,var.eq=var.eq)
  return(OUTHTCI)
}

MVN2Sample.HT.CIs.NumSum.f = function(NOBS1,XBAR1,SCOV1,NOBS2,XBAR2,SCOV2,conf.level = 0.95,
                                      alpha=.05,mu0=rep(0,length(XBAR1)),ContrastMAT=NULL,SigDig=3,var.eq=TRUE)
{
  # Contrast matrix should be given in such way that COLUMNS correspond to Contrasts
  # Computing the sample data information
  NVAR = ncol(SCOV1)
  
  # HT info:
  XBAR = XBAR1-XBAR2
  
  if(var.eq)
  {
    SPOOLED = (NOBS1-1)/(NOBS1+NOBS2-2) * SCOV1 + (NOBS2-1)/(NOBS1+NOBS2-2) * SCOV2
    SCOV=(1/NOBS1+1/NOBS2)*SPOOLED
  }else
  {
    SPOOLED = matrix(0,ncol=NVAR,nrow = NVAR)
    SCOV= 1/NOBS1*SCOV1+1/NOBS2*SCOV2
  }
  ISCOV = solve(SCOV)
  TSQOBS = t(XBAR-mu0)%*%ISCOV%*%(XBAR-mu0)
  
  HTMETHODS = c("Hotelling","Chi-Sq") 
  HTCUTOFFS = vector(length = length(HTMETHODS))
  # Hotelling/Scheffe method
  HTCUTOFFS[1] = qf(alpha,NVAR,NOBS1+NOBS2-1-NVAR, lower.tail=FALSE)
  # Asymptotic/Chi-sq Method
  HTCUTOFFS[2] =qchisq(alpha,NVAR,lower.tail=FALSE)
  
  HTPVALS = vector(length = length(HTMETHODS)) 
  HTPVALS[1] = pf((NOBS1+NOBS2-1-NVAR)/((NOBS1+NOBS2-2)*NVAR)*TSQOBS,NVAR,NOBS1+NOBS2-1-NVAR, lower.tail=FALSE)
  HTPVALS[2] = pchisq(TSQOBS,df = NVAR, lower.tail=FALSE)
  
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
  NMETH =  length(METHODS)
  # Hotelling/Scheffe method
  # Critical Value
  CUTOFFS = vector(length = length(METHODS))
  CUTOFFS[1] = sqrt(((NOBS1+NOBS2-2)*NVAR)/(NOBS1+NOBS2-1-NVAR)* qf(1-conf.level,NVAR,NOBS1+NOBS2-1-NVAR, lower.tail=FALSE))
  # Asymptotic/Chi-sq Method
  CUTOFFS[2] = sqrt(qchisq(1-conf.level,NVAR,lower.tail=FALSE))
  # Bonferroni Method
  CUTOFFS[3] = qt((1-conf.level)/(2*NCONTR),NOBS1+NOBS2-2,lower.tail=FALSE)
  
  # Compute the Intervals
  CIs = vector(mode="list", length = NMETH) 
  # lapply(1:NMETH, matrix, data= NA, nrow=NCONTR, ncol=2)
  
  SES = sqrt(diag(CCOVMAT)) 
  PES = t(ContrastMAT)%*% XBAR
  for(m in 1:NMETH)
  {
    CRITVAL = CUTOFFS[m]
    LCBS =  PES - CRITVAL*SES
    UCBS =  PES + CRITVAL*SES
    CIs[[m]] = round(data.frame(LCBS,UCBS),SigDig)
  }
  CUTOFFS[1] = qf(1-conf.level,NVAR,NOBS1+NOBS2-1-NVAR, lower.tail=FALSE)
  CUTOFFS[2] = qchisq(1-conf.level,NVAR,lower.tail=FALSE)
  names(CIs) = METHODS
  ST1 = list(NOBS1,round(XBAR1,SigDig),round(SCOV1,SigDig))
  names(ST1) = c("N","XBar","S")
  ST2 = list(NOBS2,round(XBAR2,SigDig),round(SCOV2,SigDig))
  names(ST2) = c("N","XBar","S")
  STDIFF = list(var.eq, round(XBAR,SigDig),round(SPOOLED,SigDig),round(SCOV,SigDig),round(ISCOV,SigDig))
  names(STDIFF) = c("Eq.Var","XBarDiff","SPooled","SE","ISE")
  
  RES = list(ST1,ST2,STDIFF,mu0,
             round(TSQOBS,3),HTPVALS,alpha,round(HTCUTOFFS,3),conf.level,ContrastMAT,
             round(PES,SigDig),round(diag(CCOVMAT),SigDig),
             round(CUTOFFS,3),CIs)
  names(RES) = c( "G1","G2","DIFF","Mu0","Tobs","PVals","Alpha","HT CritVal","CI conf.level",
                  "Contrasts","CMeans","CVars","CI Crit. Vals","CIs")
  return(RES)
}

###############################################################
###############################################################
# Example 5.1 Lectures (6.1 in JW)
###############################################################
###############################################################
# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
WPData = read.xlsx(xlsxFile= 'WasteWaterData.xlsx')
head(WData,5) # Displays firs 5 rows of data

#Compute Differences between the two Agencies for each Variable
# Commercial-State
DData = WData[,c(1,2)] - WData[,c(3,4)]

###############################################################
# Numerical Summaries in R
###############################################################
n=11
p = 2

XBar = apply(DData,2,mean) # Variable (column) means
# apply(DData,2,sd) # Variable (column) sds

S = cov(DData) # covariance Matrix

###############################################################
# Chi-sq QQPlot to Check Normality
###############################################################
# install.packages("heplots")
library(heplots)

cqplot(DData, id.n=3) # id the 3 most extreme values

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot
mvn(DData, mvnTest = "royston")
# Henze-Zirkler test and construct perspective plot
mvn(DData, mvnTest = "hz")

###############################################################
# HT P-val: Hotelling's T^2 test
###############################################################
library(DescTools)

mu0 = c(0,0)
HotellingsT2Test(DData, mu = mu0, test = "f")

# Note: To get the T^2_obs in the lecture we nee to  
# multiply the value of the T.2 given from function as 
T.2 = 6.1377
T2Obs = (n-1)*p*T.2/(n-p)

###############################################################
# HT RR: Hotelling's T^2 test Manual computations
###############################################################

# Computing the T^2_obs manually
IS = solve(S)

T2 = n*t(XBar-mu0)%*%IS%*%(XBar-mu0)

# F critical value
FCV = qf(0.05,p,n-p, lower.tail=FALSE)

# Threshold 
Ta2 = (n-1)*p/(n-p)*FCV 

###############################################################
# HT And Simultaneous CI using Instructor Function
###############################################################
CISigl = 0.05
RES = MVN.HT.CIs.f(DData,conf.level = 1-CISigL,alpha=.05,mu0=MuO,ContrastMAT=NULL,SigDig = NSigDig)
RES
###############################################################
# Linear Combo
###############################################################
ContCISigL = 0.10

nc = 3
# Contrast Mat with the linear combinations as columns
ContMat = matrix(nrow=p,ncol=nc) 
ContMat[,1] = 1/2*rep(1,p) # 
ContMat[,2] = c(1,-1) # 
ContMat[,3] = c(-1,1) # 

RESC = MVN.HT.CIs.f(DData,conf.level = 1-ContCISigL,alpha=.05,mu0=MuO,ContrastMAT=ContMat,SigDig = NSigDig)

###############################################################
# EigenValues and Axes Half lengths
###############################################################

# S eigenvalues
EI = eigen(S)
EI$values

# Half lenght of eigen axes

sqrt(EI$values) * sqrt(Ta2/n)



###############################################################
###############################################################
# Example 5.1 Lectures (6.1 in JW)
###############################################################
###############################################################
# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
DAData = read.xlsx(xlsxFile= 'DogAnesthesiaData.xlsx')
head(DAData,5) # Displays firs 5 rows of data

#Compute Differences between the two Agencies for each Variable
# Commercial-State
# DAData = WData[,c(1,2)] - WData[,c(3,4)]

###############################################################
# Numerical Summaries in R
###############################################################
n=18
p = 4

XBar = apply(DAData,2,mean) # Variable (column) means
# apply(DAData,2,sd) # Variable (column) sds

S = cov(DAData) # covariance Matrix

###############################################################
# Chi-sq QQPlot to Check Normality
###############################################################
# install.packages("heplots")
library(heplots)

cqplot(DAData, id.n=3) # id the 3 most extreme values

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot
mvn(DAData, mvnTest = "royston")
# Henze-Zirkler test and construct perspective plot
mvn(DAData, mvnTest = "hz")

###############################################################
# HT P-val: Hotelling's T^2 test
###############################################################
library(DescTools)

# Must take differences to perform HT.
DDAData = DAData[,2:4] - DAData[,1]
mu0 = c(0,0,0)
HotellingsT2Test(DDAData, mu = mu0, test = "f")

# Note: To get the T^2_obs in the lecture we nee to  
# multiply the value of the T.2 given from function as 
T.2 = 6.1377
T2Obs = (n-1)*p*T.2/(n-p)

###############################################################
# HT RR: Hotelling's T^2 test Manual computations
###############################################################

# Computing the T^2_obs manually
IS = solve(S)

T2 = n*t(XBar-mu0)%*%IS%*%(XBar-mu0)

# F critical value
FCV = qf(0.05,p,n-p, lower.tail=FALSE)

# Threshold 
Ta2 = (n-1)*p/(n-p)*FCV 

###############################################################
# HT And Simultaneous CI using Instructor Function
###############################################################
CISigl = 0.05
RES = MVN.HT.CIs.f(DAData,conf.level = 1-CISigL,alpha=.05,mu0=MuO,ContrastMAT=NULL,SigDig = NSigDig)
RES
###############################################################
# Linear Combo
###############################################################
ContCISigL = 0.10

nc = 3
# Contrast Mat with the linear combinations as columns
ContMat = matrix(nrow=p,ncol=nc) 
ContMat[,1] = c(-1,-1,1,1) # 
ContMat[,2] = c(1,-1,1,-1) # 
ContMat[,3] = c(1,-1,-1,1) # 

RESC = MVN.HT.CIs.f(DAData,conf.level = 1-ContCISigL,alpha=.05,mu0=MuO,ContrastMAT=ContMat,SigDig = NSigDig)

###############################################################
# EigenValues and Axes Half lengths
###############################################################

# S eigenvalues
EI = eigen(S)
EI$values

# Half lenght of eigen axes

sqrt(EI$values) * sqrt(Ta2/n)




###############################################################
###############################################################
# Example 5.1 Lectures (6.1 in JW)
###############################################################
###############################################################

###############################################################
###############################################################
# Raw Data
###############################################################
###############################################################
# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')
# Read data
WPData = read.xlsx(xlsxFile= 'WasteWaterData.xlsx')
head(WData,5) # Displays firs 5 rows of data

#Compute Differences between the two Agencies for each Variable
# Commercial-State
DData = WData[,c(1,2)] - WData[,c(3,4)]
NSigDig = 3
p = 2

# HT/CI Info
CISigL = .05 
ContCISigL = .10

colMeans(DData)
MuO= c(0,0)

# Contrast Info
VarNames = paste0("X_",1:p)
LinCombNames = c("Y","Z","W")
ContExpres = "\\frac{D_1+D_2}{2};X_1-X_2;X_2-X_1"
p = ncol(DData)
nc = length(LinCombNames)
# Contrast Mat with the linear combinations as columns
ContMat = matrix(nrow=p,ncol=nc) 
ContMat[,1] = 1/2*rep(1,p) # 
ContMat[,2] = c(1,-1) # 
ContMat[,3] = c(-1,1) # 

###############################################################
# Chi-sq QQPlot to Check Normality Untransformed data
###############################################################
# install.packages("heplots")
library(heplots)

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot
NPvalsUntr = c(.5,.5)
NPvalsUntr = c(mvn(DData, mvnTest = "royston")$multivariateNormality$`p value`,mvn(DData, mvnTest = "hz")$multivariateNormality$`p value`)
NPvalsUntr

ChSqQQPlotDataUntr =  round(ChiSq.QQPlot.Data.f(DData),4)

###############################################################
# Plot QQPlots: Looping
###############################################################

dev.off() # resets all graphical parameters 
par(mfrow = c(2,2))
for(i in 1:p){
  qqnorm(DData[, i],
         main = paste0("Q-Q plot of x", i), pch=19, cex=1.5)
  qqline(DData[, i],col = "#006CD1")
}
cqplot(DData, id.n=3)


# Box-Cox Transformation
# Estimate value of lambda
library("EnvStats")
TLambdas = rep(1,p)
DData1 = DData
par(mfrow = c(2,2))
for(i in 1:p){
  x = DData[, i]
  BCL = boxcox(x, optimize = TRUE)
  lambda = BCL$lambda # best guess of lambda
  TLambdas[i] = lambda
  x = (x^lambda-1)/lambda # transformed data
  qqnorm(x,
         main = paste0("Q-Q plot of x", i), pch=19, cex=1.5)
  qqline(x,col = "#006CD1")
  DData1[,i] = x
}
TLambdas=round(TLambdas,4)
cqplot(DData1, id.n=3)
abline(a=0,b=1,col="blue")

# Perform Royston's test Henze-Zirkler 
NPvalsUntr = c(.5,.5)
NPvalsUntr = c(mvn(DData, mvnTest = "royston")$multivariateNormality$`p value`,mvn(DData, mvnTest = "hz")$multivariateNormality$`p value`)
NPvalsUntr

XBarUntr = round(apply(DData,2,mean),NSigDig) # Variable (column) means
SUntr = round(cov(DData),NSigDig) # covariance Matrix

DData = DData1

###############################################################
# Chi-sq QQPlot to Check Normality Transformed Data
###############################################################

# install.packages("heplots")
library(heplots)
D2 = cqplot(DData, id.n=3) # id the 3 most extreme values

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot
NPvals = c(.5,.5)
NPvals = c(mvn(DData, mvnTest = "royston")$multivariateNormality$`p value`,mvn(DData, mvnTest = "hz")$multivariateNormality$`p value`)
NPvals

ChSqQQPlotData =  round(ChiSq.QQPlot.Data.f(DData),NSigDig)

RES = MVN.HT.CIs.f(DData,conf.level = 1-CISigL,alpha=.05,mu0=MuO,ContrastMAT=NULL,SigDig = NSigDig)

###############################################################
# Linear Combo
###############################################################

RESC = MVN.HT.CIs.f(DData,conf.level = 1-ContCISigL,alpha=.05,mu0=MuO,ContrastMAT=ContMat,SigDig = NSigDig)



###############################################################
###############################################################
# Example 5.2 Lectures (6.2 in JW)
###############################################################
###############################################################

###############################################################
###############################################################
# Raw Data
###############################################################
###############################################################
# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')
# Read data
DAData = read.xlsx(xlsxFile= 'DogAnesthesiaData.xlsx')
head(DAData,5) # Displays firs 5 rows of data

DDAData = DAData[,2:4] - DAData[,1]
MuO = c(0,0,0)

#Compute Differences between the two Agencies for each Variable
# Commercial-State
# DAData = WData[,c(1,2)] - WData[,c(3,4)]
NSigDig = 1
p = 4

# HT/CI Info
CISigL = .10 
ContCISigL = .20

colMeans(DAData)
MuO= c(350,400,475,500)

# Contrast Info
VarNames = paste0("X_",1:p)
LinCombNames = c("H","C","I")
ContExpres = "X_3+X_4 - (X_1+X_2);X_1+X_3 - (X_2+X_4);X_1+X_4 - (X_2+X_3)"
p = ncol(DAData)
nc = length(LinCombNames)
# Contrast Mat with the linear combinations as columns
ContMat = matrix(nrow=p,ncol=nc) 
ContMat[,1] = c(-1,-1,1,1) # 
ContMat[,2] = c(1,-1,1,-1) # 
ContMat[,3] = c(1,-1,-1,1) # 

###############################################################
# Chi-sq QQPlot to Check Normality Untransformed data
###############################################################
# install.packages("heplots")
library(heplots)

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot
NPvalsUntr = c(.5,.5)
NPvalsUntr = c(mvn(DAData, mvnTest = "royston")$multivariateNormality$`p value`,mvn(DAData, mvnTest = "hz")$multivariateNormality$`p value`)
NPvalsUntr

ChSqQQPlotDataUntr =  round(ChiSq.QQPlot.Data.f(DAData),4)

###############################################################
# Plot QQPlots: Looping
###############################################################

dev.off() # resets all graphical parameters 
par(mfrow = c(2,2))
for(i in 1:p){
  qqnorm(DAData[, i],
         main = paste0("Q-Q plot of x", i), pch=19, cex=1.5)
  qqline(DAData[, i],col = "#006CD1")
}
cqplot(DAData, id.n=3)


# Box-Cox Transformation
# Estimate value of lambda
library("EnvStats")
TLambdas = rep(1,p)
DAData1 = DAData
par(mfrow = c(2,2))
for(i in 1:p){
  x = DAData[, i]
  BCL = boxcox(x, optimize = TRUE)
  lambda = BCL$lambda # best guess of lambda
  TLambdas[i] = lambda
  x = (x^lambda-1)/lambda # transformed data
  qqnorm(x,
         main = paste0("Q-Q plot of x", i), pch=19, cex=1.5)
  qqline(x,col = "#006CD1")
  DAData1[,i] = x
}
TLambdas=round(TLambdas,4)
cqplot(DAData1, id.n=3)
abline(a=0,b=1,col="blue")

# Perform Royston's test Henze-Zirkler 
NPvalsUntr = c(.5,.5)
NPvalsUntr = c(mvn(DAData, mvnTest = "royston")$multivariateNormality$`p value`,mvn(DAData, mvnTest = "hz")$multivariateNormality$`p value`)
NPvalsUntr

XBarUntr = round(apply(DAData,2,mean),NSigDig) # Variable (column) means
SUntr = round(cov(DAData),NSigDig) # covariance Matrix

DAData = DAData1

###############################################################
# Chi-sq QQPlot to Check Normality Transformed Data
###############################################################

# install.packages("heplots")
library(heplots)
cqplot(DAData, id.n=3) # id the 3 most extreme values

# Running test stat for MVN
# install.packages("MVN")
library(MVN)
# Perform Royston's test and create a chi-square plot
NPvals = c(.5,.5)
NPvals = c(mvn(DAData, mvnTest = "royston")$multivariateNormality$`p value`,mvn(DAData, mvnTest = "hz")$multivariateNormality$`p value`)
NPvals

ChSqQQPlotData =  round(ChiSq.QQPlot.Data.f(DAData),NSigDig)

RES = MVN.HT.CIs.f(DAData,conf.level = 1-CISigL,alpha=.05,mu0=MuO,ContrastMAT=NULL,SigDig = NSigDig)

###############################################################
# Linear Combo
###############################################################

RESC = MVN.HT.CIs.f(DAData,conf.level = 1-ContCISigL,alpha=.05,mu0=MuO,ContrastMAT=ContMat,SigDig = NSigDig)



###############################################################
###############################################################
###############################################################
# 2 Sample MVN Hotelling's
###############################################################
###############################################################
###############################################################

###############################################################
# Using Numerical Summaries
###############################################################

# EX 6.3
p = 2
VarNames = paste0("X_",1:p)

N1 = 45 
XB1 = c(204.4,556.6)
S1 = matrix(c(13825.3,23823.4,23823.4,73107.4),ncol=2)
N2 = 55
XB2 = c(130.0,355.0)
S2 = matrix(c(8632.0,19616.7,19616.7,55964.5),ncol=2)

# EX6.4 
p = 2
VarNames = paste0("X_",1:p)

N1 = 50 
XB1 = c(8.3,4.1)
S1 = matrix(c(2,1,1,6),ncol=2)
N2 = 50
XB2 = c(10.2,3.9)
S2 = matrix(c(2,1,1,4),ncol=2)

CISigL = .05 
RES = MVN2Sample.HT.CIs.NumSum.f(N1,XB1,S1,N2,XB2,S2,conf.level = 1-CISigL,var.eq = TRUE)

###############################################################
# generating fake chi-sq data to populate the chi-sq QQPlot
###############################################################
# ChSqQQPlotDataUntr =  round(ChiSq.QQPlot.Data.f(DAData),4)

D2 = sort(rchisq(N1,df = p))
Q = qchisq((1:N1-.5)/N1,df = p)
plot(Q,D2)
ChSqQQPlotDataA =  round(cbind(Q,D2),4)

D2 = sort(rchisq(N2,df = p))
Q = qchisq((1:N2-.5)/N2,df = p)
plot(Q,D2)
ChSqQQPlotDataB =  round(cbind(Q,D2),4)

# Running test stat for MVN
NPvalsA = c(.5,.5)
NPvalsB = c(.5,.5)

###############################################################
# Raw Data
###############################################################
# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')
# Read data
LData = read.xlsx(xlsxFile= 'GenusLizardsData.xlsx')
head(LData,5) # Displays firs 5 rows of data

# Separate the two geni into two objects X and Y
attach(LData) # attaching the data so we can use columns names
X = LData[Genus==0,1:2]
Y = LData[Genus==1,1:2]

p = ncol(X)
VarNames = paste0("X_",1:p)

###############################################################
# Chi-sq QQPlot to Check Normality Untransformed data
###############################################################
# install.packages("heplots")
library(heplots)
# install.packages("MVN")
library(MVN)

###############################################################
# Plot QQPlots, Royston and HZ tests
###############################################################

# Group A
dev.off() # resets all graphical parameters 
par(mfrow = c(2,2))
for(i in 1:p){
  qqnorm(X[, i],
         main = paste0("Q-Q plot of x", i), pch=19, cex=1.5)
  qqline(X[, i],col = "#006CD1")
}
# MV Chi2 QQPlot
cqplot(X, id.n=3)

# Testing MVN
NPvalsA = c(.5,.5)
NPvalsA = c(mvn(X, mvnTest = "royston")$multivariateNormality$`p value`,mvn(X, mvnTest = "hz")$multivariateNormality$`p value`)
NPvalsA

ChSqQQPlotDataA =  round(ChiSq.QQPlot.Data.f(X),4)

# Group B
par(mfrow = c(2,2))
for(i in 1:p){
  qqnorm(Y[, i],
         main = paste0("Q-Q plot of x", i), pch=19, cex=1.5)
  qqline(Y[, i],col = "#006CD1")
}
cqplot(Y, id.n=3)

NPvalsB = c(.5,.5)
NPvalsB = c(mvn(Y, mvnTest = "royston")$multivariateNormality$`p value`,mvn(Y, mvnTest = "hz")$multivariateNormality$`p value`)
NPvalsB

ChSqQQPlotDataB =  round(ChiSq.QQPlot.Data.f(Y),4)

####################################################
# HT
####################################################

library(DescTools)

# Using the build in function to compute p-vals
HotellingsT2Test(X,Y, test = "f") # F Test
HotellingsT2Test(X,Y, test = "chi") # Chi-Sq

# Using instructors function that also provides CIs
RES = MVN2Sample.HT.CIs.f(X,Y,conf.level = 0.95,alpha=.05,var.eq=TRUE)
RES 


###############################################################
###############################################################
###############################################################
# MANOVA
###############################################################
###############################################################
###############################################################

###############################################################
# MANOVA: Hand calculations
###############################################################
G1 = matrix(c(3,1,4,0,5,2,8,5),ncol=2,byrow = TRUE)
G2 = matrix(c(9,7,3,5,7,8,5,4),ncol=2,byrow = TRUE)
G3 = matrix(c(10,5,9,6,8,5,9,4),ncol=2,byrow = TRUE)

XBar1 = colMeans(G1)
S1 = cov(G1)

XBar2 = colMeans(G2)
S2 = cov(G2)

XBar3 = colMeans(G3)
S3 = cov(G3)

# Overall mean 
XBarT = colMeans(rbind(G1,G2,G3))

B1 = (XBar1-XBarT)%*%t((XBar1-XBarT))
B2 = (XBar2-XBarT)%*%t((XBar2-XBarT))
B3 = (XBar3-XBarT)%*%t((XBar3-XBarT))

B = 4*B1+4*B2+4*B3 # matrix(c(34.68,17.32,17.32,34.68),ncol=2,byrow = TRUE)
W = 3*S1+3*S2+3*S3
BW = B + W

Lambda = det(W)/det(BW)

p =2
g = 3 
n = 12
a = .1

df = n - g
a/(p*g*(g-1))

ta = qt(a/(p*g*(g-1)),df = n-g,lower.tail = FALSE)

XBar1 - XBar3 - ta *sqrt((1/4+1/4)*diag(W)/(n-g))
XBar1 - XBar3 + ta *sqrt((1/4+1/4)*diag(W)/(n-g))

###############################################################
###############################################################
# Example 5.3 Lectures (6.24 in JW)
###############################################################
###############################################################
library('openxlsx')
# Read data
SData = read.xlsx(xlsxFile= 'EgyptianSkullData.xlsx')
head(SData,5) # Displays firs 5 rows of data

SData[,1] = as.factor(SData[,1]) 
###############################################################
# MANOVA: R Computations
###############################################################
# MANOVA with car package
library(car)

###############################################################
# Numerical summaries by Group
###############################################################

# Computing Means and Var/Cov Mats per groups
VMeans = statList(SData[,-1],SData[,1],FUN = colMeans)
VMat = statList(SData[,-1],SData[,1],FUN = var)
Ns = table(SData[,1]) # Sample sizes
p = 4
g = 3
n = nrow(SData)

# Manually computing the W matrix
W = (Ns[1]-1)*VMat[[1]] + (Ns[2]-1)*VMat[[2]] + (Ns[3]-1)*VMat[[3]]

###############################################################
# Chi-sq QQPlot to Check Normality Untransformed data
###############################################################
library(heplots)
library(MVN)

###############################################################
# Plot QQPlots, Royston and HZ tests
###############################################################

# Split data by group
SDataByGr = split(SData,SData[,1])

# Checking Multivariate normality for each group
dev.off() # resets all graphical parameters 
par(mfrow = c(2,2))
for(i in 1:g){
  X = SDataByGr[[i]][,-1]
  cqplot(X, id.n=3)
   print(mvn(X, mvnTest = "royston"))
   print(mvn(X, mvnTest = "hz"))
}

# MANOVA with built in function
# For some reason id doesn't like the dataframe
# So we need to split the data
Y = as.matrix(SData[,-1])
Gr = as.factor(SData[,1])

# First fit a linear regression
LM.res = lm(Y~Gr)

# Call Manova, This command gives all Tests
# Wilk's lambda, “Pillai”, “Hotelling-Lawley” and “Roy” 
SUM = summary(Manova(LM.res),"Wilks")
SUM

# The W matrix using the buld it functions
W = SUM$SSPE
W 

# Checking which of the groups differ from eachother
install.packages("biotools")
library(biotools)
mvpaircomp(LM.res,factor1 = "Gr",test = "Wilks", adjust = "bonferroni")

# Checking each variable for differences across the groups
summary.aov(LM.res)

# Simultaneous CIs, Critical value for Bonferroni
k = p*g*(g-1)/2
CISigl = .05
ta = qt(CISigl/(2*k),df=n-g,lower.tail = FALSE)

# Variable and group of interest information 
i = 2 # Group A
j = 3 # Group B
v = 1

LL = VMeans[[i]] - VMeans[[j]] - ta * sqrt((1/Ns[i]+1/Ns[j])*diag(W)/(n-g))
UL = VMeans[[i]] - VMeans[[j]] + ta * sqrt((1/Ns[i]+1/Ns[j])*diag(W)/(n-g))

cbind(LL,UL)

