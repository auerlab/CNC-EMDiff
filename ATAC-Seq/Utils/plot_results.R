#!/usr/bin/env Rscript

library(ComplexHeatmap)

####### Read in results

ne2vsce4 <- read.table("Sleuth-Prelim/ne2vsce4_v2.txt", header=TRUE, sep="\t")
ne6vsce14 <- read.table("Sleuth-Prelim/ne6vsce14_v2.txt", header=TRUE, sep="\t")

ce14vsce0 <- read.table("Sleuth-Prelim/ce14vsce0.txt", header=TRUE, sep="\t")
ce4vsce0 <- read.table("Sleuth-Prelim/ce4vsce0.txt", header=TRUE, sep="\t")
ce14vsce4 <- read.table("Sleuth-Prelim/ce14vsce4.txt", header=TRUE, sep="\t")

ne2vsne0 <- read.table("Sleuth-Prelim/ne2vsne0.txt", header=TRUE, sep="\t")
ne6vsne0 <- read.table("Sleuth-Prelim/ne6vsne0.txt", header=TRUE, sep="\t")
ne6vsne2 <- read.table("Sleuth-Prelim/ne6vsne2.txt", header=TRUE, sep="\t")

### Start with Chondro plots

ce <- rbind(ce14vsce0, ce4vsce0, ce14vsce4)

### Remove duplicate rows
ce.v2 <- ce[which(duplicated(ce$target_id)==FALSE),]

### Normalize all columns. Take log first

input.v1 <- log(as.matrix(ce.v2[, 6:14]) + 0.01)
col.mean <- apply(input.v1, 2, mean)
col.sd <- apply(input.v1, 2, sd)

input.v2 <- input.v1[, c(1,4,7,2,5,8,3,6,9)]
for(i in 1:9){
input.v2[,i] <- (input.v1[,i] - col.mean[i])/col.sd[i]
}


K = 6
kclus <- kmeans(input.v2, K)
split <- factor(paste0("Cluster ", kclus$cluster), levels = paste0("Cluster ", 1:K))


h <- Heatmap(input.v2, km = 6, cluster_columns=FALSE,
split=split)

h <- Heatmap(input.v2, cluster_columns=FALSE)


pdf("heat_temp.pdf")
draw(h)
dev.off()
