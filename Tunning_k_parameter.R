
# Data preparation
# The data set iris is used. We start by excluding the species column and scaling the data using the function scale():

# Load the data

setwd("//crc/profiles/RedirectFolders/hichem/Desktop/EMS-big data paper-March2018/R codes-EMS paper/")

source("main_for_3studies.R")

source("mapping.R")

r = read_data("MuskegonData.mat", 6:11, 14, 4, 5, 0, 1) 
r2 = read_data("Boston_dataset123.mat", 1:6, 11, 7, 8, 1, 2) 
r3 = read_data("change_and_no_change_sewi.mat", 6:22, 23, 4, 5, 0, 1) 

dataset1 = rbind( cbind(r$ch, rep(1, nrow(r$ch)) ), cbind(r$no_ch, rep(0, nrow(r$no_ch) ))) ### 
dataset2 = rbind( cbind(r2$ch, rep(1, nrow(r2$ch)) ), cbind(r2$no_ch, rep(0, nrow(r2$no_ch) ))) ### 
dataset3 = rbind( cbind(r3$ch, rep(1, nrow(r3$ch)) ), cbind(r3$no_ch, rep(0, nrow(r3$no_ch) ))) ### 

dataset1N <- fn.normalize(dataset1, 6:11, 14) 
dataset2N <- fn.normalize(dataset2, 1:6, 11) 
dataset3N <- fn.normalize(dataset3, 6:22, 23) 

k.max <- 10
data <- dataset1N

# res = tunning_k(data, k.max) 

set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(data, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:k.max

# extract wss for 2-15 clusters
# install.packages("tidyverse")
install.packages("purrr")
library(purrr)
library(tidyverse)

wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")

### 


# Install factoextra package as follow:

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")

# The remaining packages can be installed using the code below:

pkgs <- c("cluster",  "NbClust")
install.packages(pkgs)
install.packages("factoextra")

# Load packages:

library(factoextra)
library(cluster)
library(NbClust)

tunning_k <-function(data, k.max) {

# Average silhouette method
# Average silhouette method for k-means clustering
# The R code below determine the optimal number of clusters K for k-means clustering:

library(cluster)
# k.max <- 10
# data <- dataset1N
sil <- rep(0, k.max)

xdata <- data[ sample(nrow(data), size = 20000) ,]

# Compute the average silhouette width for 
# k = 2 to k = 10
for(i in 2:k.max){
  km.res <- kmeans(xdata, centers = i, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(xdata)) # rdist
  sil[i] <- mean(ss[, 3])
}

# Plot the  average silhouette width
plot(1:k.max, sil, type = "b", pch = 19, 
     frame = FALSE, xlab = "Number of clusters k")
abline(v = which.max(sil), lty = 2)

# require(cluster)
# library(cluster)
# fviz_nbclust(data, kmeans, method = "silhouette")
return(c(kop=which.max(sil), val=sil[which.max(sil)]))
}




