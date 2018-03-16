library(ggfortify)
library(plyr)
library(tidyverse)
library(stringr)

# These are all elements
toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w5_all.txt")
col. <- mapvalues(toydf$X1, from = c("c#","e#","s#"), to = c(1:3))
toydf <- toydf %>% select(-X1,-X2)
autoplot(prcomp(toydf), col = "white", loadings = FALSE, loadings.label = FALSE, loadings.label.col = "black", loadings.label.vjust = 0.8, loadings.label.hjust = 0.6) + geom_point(col = alpha(col., 0.5), size=1.5, pch=19)

# These are the rows where mammalian CNEs are and their surrogates are. Encode them with a different color
# toydf <- toydf[c(-7,-8,-18,-19),]
col.[7:8] <- 4
col.[18:19] <- 5
autoplot(prcomp(toydf), col = "white", loadings = FALSE, loadings.label = FALSE, loadings.label.col = "black", loadings.label.vjust = 0.8, loadings.label.hjust = 0.6) + geom_point(col = alpha(col., 0.5), size=1.5, pch=19)

# Do k-means clustering. Decide on the optimum number of k ommiting mammalian CNEs that are the outliers (Figure 5)
par(mfrow=c(2,2))
toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w2_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
toydf <- toydf[c(-7,-8,-18,-19),]
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 2)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w3_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
toydf <- toydf[c(-7,-8,-18,-19),]
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 3)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w4_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
toydf <- toydf[c(-7,-8,-18,-19),]
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 4)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w5_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
toydf <- toydf[c(-7,-8,-18,-19),]
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 5)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

# Do k-means clustering. Decide on the optimum number of k ommiting having all datasets (Supplementary Figure)
par(mfrow=c(2,2))
toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w2_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 2)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w3_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 3)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w4_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 4)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")

toydf <- readr::read_table("~/Dropbox/Papers/CNEs&GenomicSignatures/all-datasets_classification-clustering/GGS_ratio_w5_all.txt")
toydf$X1 <- str_replace(toydf$X1, "c#", "CNE")
toydf$X1 <- str_replace(toydf$X1, "e#", "EXON")
toydf$X1 <- str_replace(toydf$X1, "s#", "SURROGATE")
labels <- toydf$X1
toydf <- toydf %>% select(-X1, -X2) # Remove first two columns

wss <- 0
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(toydf, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
plot(1:15, wss, cex.lab = 1.2, cex.axis = 1.2, font = 2, cex.main = 1.2, font.lab = 2, font.main = 2, las = 1, type = "b", cex = 1.2, lwd = 1.2,
     main = "Optimal Number of clusters for k-means (n = 5)",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")
