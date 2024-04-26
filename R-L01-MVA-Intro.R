# This File contains codes for Performing ANOVA

# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

###############################################################
###############################################################
# Reading Data in R
###############################################################
###############################################################

# Read data
HData = read.xlsx(xlsxFile= '../Datasets/HealthData.xlsx')
head(HData,5) # Displays firs 5 rows of data
# Data sto#005AB5 in a Dataframe object, it allows for mixed type of 
# COLUMNS: Both Numeric and Categorical. No mixing WITHIN 
# each column

# Accessing part of the data
HData[2,1] # element in position (row,col) = (2,1) 
HData[1:3,] # Accessing 3 first rows of data
HData[,1:3] # Accessing 3 first columns of data
HData[,c(2,4)] # Accessing columns 2 and 4 of data

# We can access columns using their names 
HData$age

# Attach function -lets you work with the column names directly
# no need to use the "$" symbol
# Attach can only be used for datafames
attach(HData)
age
# [1] 20 25 32 24 28 35

###############################################################
###############################################################
# Numerical Summaries in R
###############################################################
###############################################################

Xbar = apply(HData,2,mean) # Variable (column) means
apply(HData,2,sd) # Variable (column) sds

S = cov(HData) # covariance Matrix
R = cor(HData) # correlation Matrix

###############################################################
##############################################################
# Scatter Plots in R
###############################################################
###############################################################

# Read data
SDData = read.xlsx(xlsxFile= '../Datasets/ScatterDotplot.xlsx')
head(SDData) # Check if data read correctly

###############################################################
# Scatterplot
###############################################################
dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 
x1 = SDData[,1]
x2 = SDData[,2]

plot(x = x1,y = x2, 
     main="Scatter Plot of Var 1 vs. Var 2",
     xlab = "Var 1",
     ylab = "Var 2",
     pch = 19, # Symbol type
     cex = 1, # Symbol size
)

###############################################################
###############################################################
# Scatter Plot and Outliers
###############################################################
###############################################################

# Read data
BPData = read.xlsx(xlsxFile= '../Datasets/BloodPressu#005AB5ata.xlsx')
head(BPData) # Check if data read correctly

###############################################################
# Histograms
###############################################################
dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 


x = BPData[,1]
range(x)
HEndPoints = seq(65,95,5)
hist(x,breaks = HEndPoints,right = FALSE,
     main="Diastolic Blood Pressure",
     xlab = "DBP",
     ylab = "No. People",
)


dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 
x = BPData[,2]
range(x)
HEndPoints = seq(120,170,10)
hist(x,breaks = HEndPoints,right = FALSE,
     main="Systolic Blood Pressure",
     xlab = "SBP",
     ylab = "No. People",
)

###############################################################
# Scatterplot
###############################################################
dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 

plot(x = BPData[,1],y = BPData[,2], 
     main="Scatter Plot of SBP vs. DPB",
     xlab = "SBP",
     ylab = "DBP",
     pch = 19, # Symbol type
     cex = 4, # Symbol size
)

###############################################################
###############################################################

###############################################################
###############################################################
# Scatter Matrix Plots in R
###############################################################
###############################################################

# Read data
SMData = read.xlsx(xlsxFile= '../Datasets/PaperQuality.xlsx')
head(SMData) # Check if data read correctly

###############################################################
# Basic ScatterMatrix
###############################################################
dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 

# Basic
pairs(SMData[,2:4], main = "Scatter Plot Matrix for Paper Quality  Dataset")

# No lower panels
pairs(SMData[,2:4], main = "Scatter Matrix for Paper Quality  Dataset",
      lower.panel=NULL)


###############################################################
# ScatterMatrix w/ Marginal Histograms and Correlations
###############################################################
# install.packages("psych")
library(psych) 

dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 

pairs.panels(SMData[,2:4], main = "Advanced Scatter Matrix for Paper Quality Dataset",
             smooth=FALSE,
             ellipses=FALSE) 

###############################################################
###############################################################

###############################################################
###############################################################
# 3D ScatterPlot in R
###############################################################
###############################################################
# install.packages("scatterplot3d") # Install
library("scatterplot3d") # load

# Read data
S3DData = read.xlsx(xlsxFile= '../Datasets/Lizards.xlsx')
head(S3DData) # Check if data read correctly

###############################################################
# 3D Scatterplot: Static
###############################################################

dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 

scatterplot3d(S3DData[,2:4],
              main="3D Scatter Plot For Lizard Data",
              xlab = "Mass",
              ylab = "SVL",
              zlab = "HLS", 
              pch = 19, # Symbol type
              cex.symbols = 2, # Symbol size
              angle = 55) # the x-y axis angle. Changes perception

###############################################################
# 3D ScatterPlot in R: Interactive
###############################################################
#install.packages("rgl")
library(rgl)
plot3d(S3DData[,2:4], type="s", size=3,
       main="3D Scatter Plot For Lizard Data",
       xlab = "Mass",
       ylab = "SVL",
       zlab = "HLS")

###############################################################
# Cool Iris data: 3d and colors
###############################################################
attach(iris) # data existing within R
head(iris)

plot3d(iris[,1:3], type="s", size=3,
       main="3D Scatter Plot For Lizard Data",
       xlab = "Sepal.Length",
       ylab = "Sepal.Width",
       zlab = "Petal.Length")

plot3d(iris[,1:3], type="s", size=3,
       main="3D Scatter Plot For Lizard Data",
       xlab = "Sepal.Length",
       ylab = "Sepal.Width",
       zlab = "Petal.Length",
       col=as.integer(iris[,5]) # use different color by species
)

###############################################################
###############################################################
# Star Plots in R
###############################################################
###############################################################

# Read data
SData = read.xlsx(xlsxFile= '../Datasets/USairpollution.xlsx')
head(SData) # Check if data read correctly

###############################################################
# Star Plots
###############################################################

dev.off() # resets all graphical parameters, 
# it will give an error if you have not created a graph before this one 

palette(rainbow(7, s = 0.6, v = 0.75)) # Define colors for each of 12 vars
stars(SData[,2:8],draw.segments = TRUE,
              main="Star Plots For Pollution Data",
              len = 0.8, key.loc = c(14, 1.75),
              labels=SData[,1], # names of points
              ncol=6
           )


