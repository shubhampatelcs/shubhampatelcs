# This File contains codes for Performing Inference PCA Analysis in R
# Load necessary libraries
install.packages("corrr")
library('corrr')

install.packages("FactoMineR")
install.packages("factoextra")

library(ggplot2)

library("FactoMineR") 
library("factoextra")

# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
MPDataTemp = read.xlsx(xlsxFile= 'ProteinData.xlsx')
head(MPDataTemp,5) # Displays firs 5 rows of data
# You can get quick info of the data
# str(MPDataTemp)

# Extract only numerical variables for the analysis
MPData = MPDataTemp[,2:10]

# Labeling every Observation (row) with the country name 
# Useful for graphing purposes 
rownames(MPData) = MPDataTemp[,1] 
head(MPData,5) # Displays firs 5 rows of data

# Check for missing values. Rows with missing vals should be removed
colSums(is.na(MPData))

# Keep only complete cases
MPData = MPData[complete.cases(MPData),]
head(MPData,5) # Displays firs 5 rows of data

# Visualize variables to check variances
boxplot(MPData, main="Distributions of Vars in Protein Data")

# Compute numerical summaries for each varable
library(tableone)
MPDataNumSum = CreateTableOne(data= MPData,includeNA = FALSE)
print(MPDataNumSum) 

# Check original bi-variate scatterplots to see if obs (countries) can be grouped 
# Looking at regions where each country belongs to

pairs(MPData,main="Protein", pch=19)

VNames = names(MPData)
par(mfrow=c(3,4))
k=0
for(i in 1:8)
  for(j in (i+1):9)
  {
    plot(MPData[,i],MPData[,j],xlab = VNames[i],ylab = VNames[j], pch=19, col=as.numeric(Region)+1)
    text(MPData[,i],MPData[,j], labels=rownames(MPData), cex= 0.7, pos=3)
  }

# Perform PCA analysis on standardized data
PCA.Res = PCA(MPData, graph = FALSE,ncp = 9,scale.unit = TRUE)

# Getting the eigenvectors/loadings 
sweep(PCA.Res$var$coord,2,sqrt(PCA.Res$eig[1:ncol(PCA.Res$var$coord),1]),FUN="/")

# Alternative way of Computing loadings
PCA.Res1  = prcomp(MPData,scale. = TRUE)
round(PCA.Res1$rotation,2) 

# See what information the PCA.Res object 
PCA.Res 
str(PCA.Res)

# Summary of PCA analysis
summary(PCA.Res)

# Scree plot to visualize variance explained by each principal component
fviz_eig(PCA.Res, addlabels = TRUE, ylim = c(0, 50))

# Variable contribution to each PC:
C1 = PCA.Res$var$contrib
round(C1,1) 

# Quality of representation: Contribution of Vars to first 2 PCs
fviz_pca_var(PCA.Res, col.var = "cos2",
gradient.cols =  c("black","#00AFBB", "#E7B800"),
repel = TRUE, # Avoid text overlapping (slow if many points),
)+ labs(title ="Correlation Plot", x = "PC1", y = "PC2")

# Alternative way to see the contributions of each var to each PC
library("corrplot")
corrplot(PCA.Res$var$contrib, is.corr=FALSE, col= c("#C994C7","#DF65B0","#E7298A"))

# Contributions of variables to PC1
# Axes gives the PC component
# top gives number of variables included in graph
fviz_contrib(PCA.Res, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(PCA.Res, choice = "var", axes = 2, top = 10)
# Can plot contributions to multiple PCs at once
fviz_contrib(PCA.Res, choice = "var", axes = 1:2, top = 10)

# Extract scores for each observation, i.e., coordinates on new axes
PCA.Res$ind$coord

# Grouping observations by geographic region 
Country = c('Albania','Austria','Belgium','Bulgaria','Czechoslovakia','Denmark',
            'E Germany','Finland','France','Greece','Hungary','Ireland','Italy',
            'Netherlands','Norway','Poland','Portugal','Romania','Spain','Sweden',
            'Switzerland','UK','USSR','W Germany','Yugoslavia')
Region = factor(c("Balkans","West Europe","West Europe","Balkans","East Europe","Scandinavia",
                  "East Europe","Scandinavia","West Europe","Mediteranean","East Europe","West Europe","Mediteranean",
                  "Scandinavia","Scandinavia","East Europe","Mediteranean","Balkans","Mediteranean","Scandinavia",
                  "West Europe","West Europe","East Europe","West Europe","Balkans"))

PCA.Region = factor(c("Balkans","West Europe","West Europe","Balkans","East Europe","West Europe",
                      "East Europe","West Europe","West Europe","Mediteranean","East Europe","West Europe","Mediteranean",
                      "West Europe","West Europe","East Europe","Mediteranean","Balkans","Mediteranean","West Europe",
                      "West Europe","West Europe","East Europe","West Europe","Balkans"))

# Visualize the observations on the first 2 PCs 
fviz_pca_ind(PCA.Res,repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(PCA.Res, habillage=Region,
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(PCA.Res, habillage=PCA.Region,
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(PCA.Res, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

# Check original bi-variate scatterplots to see if obs (countries) can be grouped 
# Looking at regions where each country belongs to
  
pairs(MPData,main="Protein", pch=19)
pairs(MPData,main="Protein", pch=19, col=as.numeric(Region)+1)
pairs(MPData,main="Protein", pch=19, col=as.numeric(PCA.Region)+1)

VNames = names(MPData)
par(mfrow=c(3,4))
k=0
for(i in 1:8)
  for(j in (i+1):9)
  {
  plot(MPData[,i],MPData[,j],xlab = VNames[i],ylab = VNames[j], pch=19, col=as.numeric(Region)+1)
  text(MPData[,i],MPData[,j], labels=rownames(MPData), cex= 0.7, pos=3)
}


