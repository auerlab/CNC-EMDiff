#!/usr/bin/env Rscript

library(ComplexHeatmap)
library(tcltk)

##########################################################################
#   Pause until user presses enter
##########################################################################

pause <- function()
{
    # Change to TRUE for interactive run
    cat("Press Enter to continue...")
    invisible(b <- scan("stdin", character(), nlines=1, quiet=TRUE))
}

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

# Merge rows of all chondro counts
ce <- rbind(ce14vsce0, ce4vsce0, ce14vsce4)
head(ce)
dim(ce)

### Remove duplicate rows
ce.v2 <- ce[which(duplicated(ce$target_id)==FALSE),]
dim(ce.v2)

### Normalize all columns. Take log first

# Extract columns with counts.  +0.01 to avoid zeros?
input.v1 <- log(as.matrix(ce.v2[, 6:14]) + 0.01)
col.mean <- apply(input.v1, 2, mean)
col.sd <- apply(input.v1, 2, sd)

print("col.mean:")
print(col.mean)

print("col.sd")
print(col.sd)

# Need to create input.v2 before the loop below can alter its values
input.v2 <- input.v1[, c(1,4,7,2,5,8,3,6,9)]
print("input.v2:")
head(input.v2)

# Normalize by dividing by mean/sd (subtracting from log is dividing)
for(i in 1:9) {
    input.v2[,i] <- (input.v1[,i] - col.mean[i])/col.sd[i]
}
print("Normalized input.v2:")
head(input.v2)

### Compute k-means

# Determine K (the number of clusters) by trial and error
K = 6
print("Computing K-means...")
kclus <- kmeans(input.v2, K)
print("kclus:")
head(kclus)

# ??
split <- factor(paste0("Cluster ", kclus$cluster),
		levels = paste0("Cluster ", 1:K))
head(split)

print("Generating heatmap...")
h <- Heatmap(input.v2, km = 6, cluster_columns=FALSE, split=split)
# h <- Heatmap(input.v2, cluster_columns=FALSE)

# Force plot to use display even though we're in Rscript
# X11() for typical Unix (BSD, Linux, etc)
# windows() for Windows
# quartz() for mac
# pdf("heat_temp.pdf")
X11()
draw(h)
pause()
dev.off()
