# Load necessary libraries
# install.packages("corrr")
library('corrr')

# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
SPData = read.xlsx(xlsxFile= 'StockPriceData.xlsx')
head(SPData,5) # Displays firs 5 rows of data
# You can get quick info of the data
# str(SPDataTemp)

# Check for missing values. Rows with missing vals should be removed
colSums(is.na(SPData))

# Keep only complete cases
SPData = SPData[complete.cases(SPData),]
head(SPData,5) # Displays firs 5 rows of data

# Visualize variables to check variances
dev.off()
boxplot(SPData, main="Distributions of Vars in Stock Data")

# Compute numerical summaries for each variable
library(tableone)
SPDataNumSum = CreateTableOne(data= SPData,includeNA = FALSE)
print(SPDataNumSum) 

# Standardize data
SPData = as.data.frame(scale(SPData))

# Check transformation
SPDataNumSum = CreateTableOne(data= SPData,includeNA = FALSE)
print(SPDataNumSum) 

# Compute Correlation Matrix 
R = cor(SPData)
round(R,3)

# Correlations Between the vars in the study
library("corrplot")
dev.off()
corrplot(R, is.corr=TRUE, col= c("#C994C7","#DF65B0","#E7298A"),tl.col="black")

# Computer the Loading and Variable Specific/Unique Matrices
library(psych)

# Using PCs 
Eig.Res = eigen(R)
# L Matrix for m=2, the minus is so it matches the ones in the book
L = round(Eig.Res$vectors[,1:2]%*% diag(Eig.Res$values[1:2]^.5),3)
PSI = diag(diag(R - L%*%t(L)))
# Error Matrix 
EPS = R - L%*%t(L) - PSI 

# Communalities: row sum of squares
H = round(rowSums(L^2),3)

# Using the mle method
FA.Res = fa(SPData, nfactors = 2, rotate = "none", fm = "ml")
L1 = FA.Res$loadings[,1:2]
PSI1 = diag(FA.Res$uniquenesses)
# Error Matrix 
EPS1 = R - L1%*%t(L1) - PSI1 

# Communalities: row sum of squares
H1 = FA.Res$communality

# Comparison of Rotations
# install.packages("GPArotation") # needed for the Quartimax Roation
library(GPArotation)
FA.Res.None = fa(SPData, nfactors = 2, rotate = "none", fm = "ml")
FA.Res.Varimax = fa(SPData, nfactors = 2, rotate = "varimax", fm = "ml")
FA.Res.Quartimax = fa(SPData, nfactors = 2, rotate = "quartimax", fm = "ml")

par(mfrow = c(1, 3))
# The argument cut = 0.3 specifies that any loading less than 0.3 (in absolute value) 
# will not be drawn, simple = F draws all the loadings (not just the largest 
# loading per variable).

fa.diagram(FA.Res.None, cut = 0.3, simple = F, main = "No Rotation")
fa.diagram(FA.Res.Varimax, cut = 0.3, simple = F, main = "Varimax Rotation")
fa.diagram(FA.Res.Quartimax, cut = 0.3, simple = F, main = "Quartimax Rotation")

####################################################################################
####################################################################################
####################################################################################
# Factor Selection Analysis using Covariance/Correlation Matrix
####################################################################################
####################################################################################
####################################################################################

# Read data
DData = read.xlsx(xlsxFile= 'DecathlonData.xlsx',colNames = TRUE, rowNames = TRUE)
head(DData,5) # Displays firs 5 rows of data
# You can get quick info of the data
# str(SPDataTemp)

# Compute Correlation Matrix 
R = as.matrix(DData)

# Correlations Between the vars in the study
library("corrplot")
corrplot(R, is.corr=TRUE, col= c("#C994C7","#DF65B0","#E7298A"),tl.col="black")

# Test to find the appropriate number of factors: 
pvals = c()
for(k in 1:5)
  pvals[k] = factanal(covmat = R,n.obs = 280,factors = k)$PVAL

pvals

# Varimax rotation to the resulting model
FA.Res4 = fa(R, nfactors = 4, rotate = "varimax", fm = "ml")
L4 = FA.Res4$loadings[,1:4]
PSI4 = diag(FA.Res4$uniquenesses)
# Error Matrix 
EPS4 = R - L4%*%t(L4) - PSI4
round(EPS4,3)
# Communalities: row sum of squares
H4 = FA.Res4$communality

# Visualize factors
fa.diagram(FA.Res4, cut = 0.32, simple = F, main = "Decathlon Data")














