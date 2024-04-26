# This File contains codes for Performing Hierarchical CA
################################################################################
################################################################################
# Toy Example with Two variables 
################################################################################
################################################################################

# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
EData = read.xlsx(xlsxFile= 'CAExampleData.xlsx')
head(EData,5) # Displays firs 5 rows of data
# You can get quick info of the data
# str(MPDataTemp)

# Check for missing values. Rows with missing vals should be removed
colSums(is.na(EData))

# Keep only complete cases
EData = EData[complete.cases(EData),]
head(EData,5) # Displays firs 5 rows of data

# Visualize variables to check variances
dev.off()
boxplot(EData, main="Distributions of Vars in Example Data")

# Compute numerical summaries for each varable
library(tableone)
MPDataNumSum = CreateTableOne(data= EData,includeNA = FALSE)
print(MPDataNumSum) 

# Visualize data to check potential groupings: 
# A) Scatterplot only if data bivariate
plot(EData[,1],EData[,2], 
     main="Scatterplot of Example Data", xlab = expression(x[1]), ylab = expression(x[2]), 
     panel.first = grid(),pch=19, xlim = c(0, 25), ylim = c(0, 25))
text(EData[,1],EData[,2], labels=1:nrow(EData), cex= 0.7, pos=3)

# B) PCA plot PC1 vs. PC2
# Perform PCA analysis on standardized data
PCA.Res = PCA(EData, graph = FALSE,ncp = 9,scale.unit = TRUE)
# Plot points on the fits two PCs
dev.off()
fviz_pca_ind(PCA.Res,repel = TRUE # Avoid text overlapping (slow if many points)
)

################################################################################
################################################################################
# Classification Analysis
################################################################################
################################################################################

# Standardize data
EData.STD = as.data.frame(scale(EData))

# Compute the matrix with the distances. You can choose the distance 
# method of you choice, but usually we use the Euclidean,
Dist.Mat = dist(x = EData.STD, method = "euclidean")

# Run clustering algorithm. We can use any of the methods we learned in class
# Here I use the Nearest Neighbor (single linkage)
EClust.NN.Res = hclust(d = Dist.Mat, method = "single")

# See what info was created from the function
names(EClust.NN.Res)
#History of the clustering; minuses mean original variables;
#  no minus means the number of a past cluster
EClust.NN.Res$merge 

# Distances of points/clusters merged
EClust.NN.Res$height

cbind(EClust.NN.Res$merge,EClust.NN.Res$height)
  
# Visualizing the classes and the merging history
# We can visualize the classification effect of particular 
# distance thresholds will have by adding the relevant horizontal lines
dev.off()
plot(EClust.NN.Res, main="Hierarchical Tree Diagram \n (Dendrogram) of Example Data")
abline(h = c(.4,1), lty = "dashed", lwd = 2)

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.NN.Res, main="4 Class Dendrogram of Example Data")
EClust.NN.Res4 = rect.hclust(EClust.NN.Res, k = 4, border = "red")

# Seeing the classes
EClust.NN.Res4

# Cluster memberships
# We can impose specific distance thresholds to obtain the corresponding clusters
EClust.NN.ResD1 = cutree(tree = EClust.NN.Res, h = 1)
EClust.NN.ResD1

# We can impose specific number of clusters
EClust.NN.Res4 = cutree(tree = EClust.NN.Res, k = 4)
EClust.NN.Res4

# Compare resulting clusters of the above splits
data.frame(EData,EClust.NN.ResD1, EClust.NN.Res4)
  
# Scatterplots with the specific clusters 

# A) Scatterplot only if data bivariate
dev.off()
plot(EData[,1],EData[,2], 
     main="Cluster Scatterplot of Example Data", xlab = expression(x[1]), ylab = expression(x[2]), 
     panel.first = grid(),pch=19, xlim = c(0, 25), ylim = c(0, 25),col = EClust.NN.Res4)
text(EData[,1],EData[,2], labels=1:nrow(EData), cex= 0.7, pos=3)

# B) PCA plot PC1 vs. PC2
# Plot points on the fits two PCs
dev.off()
fviz_pca_ind(PCA.Res,habillage=factor(EClust.NN.Res4), repel = TRUE # Avoid text overlapping (slow if many points)
)


################################################################################
# Using other agglomerative methods:
################################################################################

################################################################################
# Furthest Neighborhood (Complete)
EClust.FN.Res = hclust(d = Dist.Mat, method = "complete")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.FN.Res, main="4 Class Dendrogram: Furthest Neighbor")
EClust.FN.Res4 = rect.hclust(EClust.FN.Res, k = 4, border = "red")

################################################################################
# Centroid
# CAUTION: Centroid needs SQUARE DISTANCES 
EClust.C.Res = hclust(d = Dist.Mat^2, method = "centroid")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.C.Res, main="4 Class Dendrogram: Centroid")
EClust.C.Res4 = rect.hclust(EClust.C.Res, k = 4, border = "red")

################################################################################
# Average
EClust.Av.Res = hclust(d = Dist.Mat, method = "average")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.Av.Res, main="4 Class Dendrogram: Average")
EClust.Av.Res4 = rect.hclust(EClust.Av.Res, k = 4, border = "red")

################################################################################
# Ward's Method/Distance
EClust.W.Res = hclust(d = Dist.Mat, method = "ward.D")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.W.Res, main="4 Class Dendrogram: Ward's")
EClust.W.Res4 = rect.hclust(EClust.W.Res, k = 4, border = "red")





################################################################################
################################################################################
# Protein Data
################################################################################
################################################################################

# Reading Data in R
# install.packages('openxlsx',dependencies=TRUE)
library('openxlsx')

# Read data
PDataTemp = read.xlsx(xlsxFile= 'ProteinData.xlsx')
head(PDataTemp,5) # Displays firs 5 rows of data
# You can get quick info of the data
# str(MPDataTemp)

# Extract only numerical variables for the analysis
PData = PDataTemp[,2:10]

# Labeling every Observation (row) with the country name 
# Useful for graphing purposes 
rownames(PData) = PDataTemp[,1] 
head(PData,5) # Displays firs 5 rows of data

# Check for missing values. Rows with missing vals should be removed
colSums(is.na(PData))

# Visualize variables to check variances
dev.off()
boxplot(PData, main="Distributions of Vars in Protein Data")

# Compute numerical summaries for each varable
library(tableone)
MPDataNumSum = CreateTableOne(data= PData,includeNA = FALSE)
print(MPDataNumSum) 

# Visualize data to check potential groupings: 
# A) Scatterplot only if data bivariate
# Data have multiple vars so no need for bivariate scatterplot

# B) PCA plot PC1 vs. PC2
# Perform PCA analysis on standardized data
PCA.Res = PCA(PData, graph = FALSE,ncp = 9,scale.unit = TRUE)
# Plot points on the fits two PCs
dev.off()
postscript(outfilename,horizontal = TRUE)
fviz_pca_ind(PCA.Res,repel = TRUE # Avoid text overlapping (slow if many points)
)

################################################################################
################################################################################
# Classification Analysis
################################################################################
################################################################################

# Standardize data
PData.STD = as.data.frame(scale(PData))

# Compute the matrix with the distances. You can choose the distance 
# method of you choice, but usually we use the Euclidean,
Dist.Mat = dist(x = PData.STD, method = "euclidean")

# Run clustering algorithm. We can use any of the methods we learned in class
# Here I use the Nearest Neighbor (single linkage)
EClust.NN.Res = hclust(d = Dist.Mat, method = "single")

# See what info was created from the function
names(EClust.NN.Res)
#History of the clustering; minuses mean original variables;
#  no minus means the number of a past cluster
EClust.NN.Res$merge 

# Distances of points/clusters merged
EClust.NN.Res$height

cbind(EClust.NN.Res$merge,EClust.NN.Res$height)

# Visualizing the classes and the merging history
# We can visualize the classification effect of particular 
# distance thresholds will have by adding the relevant horizontal lines
dev.off()
plot(EClust.NN.Res, main="Hierarchical Tree Diagram \n (Dendrogram) of Protein Data")
abline(h = 2.5, lty = "dashed", lwd = 2)

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.NN.Res, main="5 Class Dendrogram of Protein Data")
EClust.NN.Res5 = rect.hclust(EClust.NN.Res, k = 5, border = "red")

# We can impose specific number of clusters
EClust.NN.Res5 = cutree(tree = EClust.NN.Res, k = 5)
EClust.NN.Res5

# Scatterplots with the specific clusters 

# A) Scatterplot only if data bivariate
# Not applicable 

# B) PCA plot PC1 vs. PC2
# Plot points on the fits two PCs
dev.off()
postscript(outfilename,horizontal = TRUE)
fviz_pca_ind(PCA.Res,habillage=factor(EClust.NN.Res5), repel = TRUE # Avoid text overlapping (slow if many points)
)

################################################################################
# Using other agglomerative methods:
################################################################################

################################################################################
# Furthest Neighborhood (Complete)
EClust.FN.Res = hclust(d = Dist.Mat, method = "complete")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.FN.Res, main="Protein Data 5 Class Dendrogram: Furthest Neighbor")
EClust.FN.Res5 = rect.hclust(EClust.FN.Res, k = 5, border = "red")

################################################################################
# Centroid
# CAUTION: Centroid needs SQUARE DISTANCES 
EClust.C.Res = hclust(d = Dist.Mat^2, method = "centroid")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.C.Res, main="Protein Data 5 Class Dendrogram: Centroid")
EClust.C.Res5 = rect.hclust(EClust.C.Res, k = 5, border = "red")

################################################################################
# Average
EClust.Av.Res = hclust(d = Dist.Mat, method = "average")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.Av.Res, main="Protein Data 5 Class Dendrogram: Average")
EClust.Av.Res5 = rect.hclust(EClust.Av.Res, k = 5, border = "red")

################################################################################
# Furthest Neighborhood (Complete)
EClust.W.Res = hclust(d = Dist.Mat, method = "ward.D")

# We can impose specific number of clusters in the graph
dev.off()
plot(EClust.W.Res, main="Protein Data 5 Class Dendrogram: Ward's")
EClust.W.Res5 = rect.hclust(EClust.W.Res, k = 5, border = "red")


################################################################################
################################################################################
################################################################################
################################################################################
# Non-Hierarchical Classification Analysis: K-Means
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
# USArrests Data
################################################################################
################################################################################
# Data exist in R already We just need to evoke them
data("USArrests") 
AData = USArrests
head(AData,5) # Displays firs 5 rows of data

# Check for missing values. Rows with missing vals should be removed
colSums(is.na(AData))

# Keep only complete cases
AData = AData[complete.cases(AData),]
head(AData,5) # Displays firs 5 rows of data

# Visualize variables to check variances
dev.off()
boxplot(AData, main="Distributions of Vars in USArrests Data")

# Compute numerical summaries for each varable
library(tableone)
MADataNumSum = CreateTableOne(data= AData,includeNA = FALSE)
print(MADataNumSum) 

# Visualize data to check potential groupings: 
# A) Scatterplot only if data bivariate
# Data have multiple vars so no need for bivariate scatterplot

# B) PCA plot PC1 vs. PC2
# Perform PCA analysis on standardized data
dev.off()
PCA.Res = PCA(AData, graph = FALSE,ncp = 9,scale.unit = TRUE)
fviz_pca_ind(PCA.Res,labelsize = 5,repel = TRUE)+ # Avoid text overlapping (slow if many points)
  theme(text = element_text(size = 45),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 30))


################################################################################
################################################################################
# Classification Analysis to get an idea of the number of clusters to use
################################################################################
################################################################################

# Standardize data
AData.STD = as.data.frame(scale(AData))

# Compute the matrix with the distances. You can choose the distance 
# method of you choice, but usually we use the Euclidean,
Dist.Mat = dist(x = AData.STD, method = "euclidean")

# Run clustering algorithm. We can use any of the methods we learned in class
# Here I use the Nearest Neighbor (single linkage)
AClust.NN.Res = hclust(d = Dist.Mat, method = "single")

# Visualize Dentrogram with 6 Clusters
dev.off()
plot(AClust.NN.Res, main="6 Class Dendrogram of US Arrest Data")
rect.hclust(AClust.NN.Res, k = 6, border = "red")


################################################################################
################################################################################
# K-means clustering using Standardized data
################################################################################
################################################################################

set.seed(13) # Save seed to reproduce analyses
# Scree/Elbow method for kmeans using W =SSW/SST
dev.off()
fviz_nbclust(AData.STD, kmeans, method = "wss",k.max = 10)

# Silhouette Graph
dev.off()
fviz_nbclust(AData.STD, kmeans, method = "silhouette",k.max = 10)

# Running the K-Means Algorith with 5 clusters
AKMClust.Res5 = kmeans(x = AData.STD, centers = 5, nstart = 15)
AKMClust.Res5

# Visualize data on PC1, PC2 plane
dev.off()
fviz_cluster(AKMClust.Res5, data = AData.STD,labelsize = 15,repel = TRUE)

# Running the K-Means Algorith with 6 clusters
AKMClust.Res6 = kmeans(x = AData.STD, centers = 6, nstart = 15)
AKMClust.Res6

# Visualize data on PC1, PC2 plane
dev.off()
fviz_cluster(AKMClust.Res6, data = AData.STD,labelsize = 15,repel = TRUE)

################################################################################
################################################################################
# PAM clustering using Standardized data
################################################################################
################################################################################

# Scree/Elbow method for pam using W =SSW/SST
set.seed(13)

dev.off()
fviz_nbclust(AData.STD, pam, method = "wss",k.max = 10)
 
# Silhouette Graph
dev.off()
fviz_nbclust(AData.STD, pam, method = "silhouette",k.max = 10)

# Running the PAM Algorithm with 5 clusters
set.seed(8912)
APAMClust.Res5 = pam(x = AData.STD, k = 5, nstart = 15)
APAMClust.Res5

dev.off()
fviz_cluster(APAMClust.Res5, data = AData.STD,labelsize = 15,repel = TRUE)

# Running the K-Means Algorithm with 6 clusters
library(cluster)
set.seed(8912)
APAMClust.Res6 =  pam(x = AData.STD, k = 6, nstart = 15)
APAMClust.Res6

dev.off()
fviz_cluster(APAMClust.Res6, data = AData.STD,labelsize = 15,repel = TRUE)











