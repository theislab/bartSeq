#################################################################################################################################
#
# Project: BART-Seq
# Author: Fatma Uzbas
#
#############################################   PROTECTION GROUP IDENTIFICATION   ###############################################
#
# This script counts the protection groups that were added to the 5' of the BART-Seq barcodes in all possible NNN combinations.
# 
# The script uses the RData file "results_motifs.RData", which was generated from the next-generation sequencing run deposited
# to the NCBI-GEO with the sample ID: GSM2877059.
#
#
# Explanation:
# (Only the data from Library 3 is relevant for this analysis)
# 
# 1 - The list of protection groups was subsetted, either left or right barcodes
# 2 - An empty list of data frames was created
# 3 - Each iteration of the for loop processes one mutation/wt (17 in total), subsets the trinucleotide protection groups, and creates a data frame
# 4 - Then, the read counts assigned to each trinucleotide protection group was summed across the 17 data frames 
#
#################################################################################################################################




# Working directory
WD = "D:/Users/fatma.uzbas/Desktop/Protection Groups"
setwd(WD)

# Load the RData
load("D:/Users/fatma.uzbas/Desktop/Protection Groups/results_motifs.RData")


# install.packages("gplots")
# install.packages("plyr")
# install.packages("xlsx")


library(gplots)
library(plyr)
library(xlsx)



# Left barcodes 

lefts = results$library3$protection.groups$left

dflistL = list()  

for (mut in names(lefts)) {
    df = lefts[[mut]][,which(nchar(colnames(lefts[[mut]])) == 3)]
    dflistL[[mut]] = df
}

sumLall = aaply(laply(dflistL, as.matrix), c(2,3), sum)

sumL = sumLall[1:10,]   # Subset only used barcodes
rownames(sumL) = c("A","B","C","D","E","F","G","H","I","J")   # Rename with forward barcode IDs



# Right barcodes

rights = results$library3$protection.groups$right

dflistR = list()

for (mut in names(rights)) {
  df = rights[[mut]][,which(nchar(colnames(rights[[mut]])) == 3)]
  dflistR[[mut]] = df
}

sumRall = aaply(laply(dflistR, as.matrix), c(2,3), sum)

sumR = sumRall[1:8,]   # Subset only used barcodes



# Save the dataframe as excel file
write.xlsx2(as.data.frame(rbind(sumL,sumR)), "NGS4_Lib3 - Protection groups - SumAllMuts.xlsx", col.names=TRUE, row.names=TRUE, append=TRUE)



# Create a heatmap of barcodes and protection groups

pdf("NGS4_Lib3_Fatma.pdf", width=18, height=4)

scaledL = log10(sumL[nrow(sumL):1, ])
heatmap.2(scaledL, main="Forward", scale="row", trace="none", dendrogram="none", Rowv=FALSE, lwid=c(1,11), lhei=c(2,7), key.title=NA, key.xlab=NA, key.ylab=NA,
          density.info="none", Colv=FALSE, col=colorRampPalette(c("white","orangered"))(256))

scaledR = log10(sumR[nrow(sumR):1, ])
heatmap.2(scaledR, main="Reverse", scale="row", trace="none", dendrogram="none", Rowv=FALSE, lwid=c(1,11), lhei=c(2,6), key.title=NA, key.xlab=NA, key.ylab=NA,
          density.info="none", Colv=FALSE, col=colorRampPalette(c("white","orangered"))(256))

dev.off()


