#################################################################################################################################
#
#  Project: BART-Seq
#  Author: Fatma Uzbas
#
################################   NORMALIZE THE BART-SEQ DATA   ###################################################
#
#  This script normalizes a sample BART-Seq read count matrix using scaling factors calculated based on the spike-ins
# 
#  The main steps are:
#  1 - Adjusting the outlier spike-in reads and bad barcode-spike-in combinations
#  2 - Calculating the scaling factors: 2^Mean(log2(1+spike-inN)))
#  3 - Normalizing the data
#
###################################################################################################################

install.packages("readxl")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("gplots")
install.packages("GGally")
install.packages("xlsx")

library(readxl)
library(reshape2)
library(ggplot2)
library(gplots)  # For heatmaps
library(GGally)  # For correlograms
library(xlsx)  # For exporting the data to excel

###################################################################################################################

setwd("...your_working_directory.../normalize_readcounts/")

ExcelFile = "...your_working_directory.../normalize_readcounts/BARTseq_sample_count_matrix.xlsx"   # Sample data
NGS = "SampleData"   # Name of the experiment
lib = "Lib1"   # Library ID

# Read in the original data

dat0 = read_excel(ExcelFile, sheet = lib, col_names = TRUE, col_types = NULL, na = "", skip = 0)

genes = c("B2M","CCND1","CCNE1","CER1","DNMT3B","GAPDH","LIN28A","NANOG","POU5F1","SOX2","ZFP42","RNA1","RNA2","RNA6","RNA8")
spikes = c("RNA1","RNA2","RNA6","RNA8")

# Replace any extremely high spike-ins reads, if there are any (none in this case)

dat = dat0

#################################################################################################################
################################   RAW DATA   ###################################################################
#################################################################################################################

# Heatmaps of the raw data as pdf (organized as plates)

divideR = c(0,16,32,48)
divideC = c(0,24,48)

pdf(paste0(NGS,"_",lib," Heatmaps (Raw).pdf"), width=14,height=10)
for (gene in genes) {
  try(heatmap.2(log2(acast(dat0, well.y~well.x, value.var=gene)+1), Rowv=F, Colv=F, dendrogram='none',
                trace="none", density.info="density", key.title=NA, cexRow=1.2, cexCol=1.2, adjCol=c(0.9,0.5),
                main=paste0("\n\n\n",NGS," - ",lib,"\n",gene, "\nRaw"), rowsep=divideR, colsep=divideC, sepwidth=c(0.0001,0.0001), sepcolor="black",
                labRow = acast(dat0, well.y~well.x, value.var="left")[,1], 
                labCol = acast(dat0, well.y~well.x, value.var="right")[1,],
                col=colorRampPalette(c("white","Azure","PowderBlue","SteelBlue","MidnightBlue"))(300), lwid=c(2,10), lhei=c(2,10),
                cellnote=acast(dat0, well.y~well.x, value.var="cell.count"), notecol="gray30", notecex=0.7))
}
dev.off()


#################################################################################################################
################################   SCALING FACTOR: RNA_X  (RAW DATA)  ###########################################
#################################################################################################################

# Check how the bad barcode-spike-in combinations influence the scaling factors

# Calculate the scaling factor RNA_X

datclean = dat
datclean$RNA_X = 2^rowMeans(log2(datclean[spikes]+1))-1


# Boxplots of the raw spike-ins and the scaling factor (black)

pdf(paste0(NGS,"_",lib," Spike-in BoxPlots (Raw).pdf"), width=10, height=3)
for (bc in c("left","right")){
  
  box.melt = melt(datclean[c(bc,spikes,"RNA_X")], id.vars=c(bc))
  box.melt$value = log2(box.melt$value+1)
  
  plot(
    ggplot(box.melt, aes_string(x=bc, y="value")) +
      geom_jitter(data=subset(box.melt, variable != "RNA_X"), alpha=0.8, size=0.3, aes_string(color="variable")) +
      geom_jitter(data=subset(box.melt, variable == "RNA_X"), alpha=1, size=0.4) +
      geom_boxplot(data=subset(box.melt, variable == "RNA_X"), outlier.shape=NA, alpha=0, lwd=0.3) +
      labs(title=paste0(bc, " spike-ins (raw)"), x="", y="Raw reads (log2)") +
      theme_classic() +
      ylim(0,NA) +
      theme(plot.title = element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
      guides(colour=guide_legend(override.aes=list(size=4)), fill=guide_legend(title="Spike-in")))
  }
dev.off()


#################################################################################################################
################################   ADJUST THE BAD COMBINATIONS   ################################################
#################################################################################################################

# Calculate tha bad barcode-spike-in combinations

datclean = dat

spikes.all = data.frame(dat[c("left","right","well.y","well.x",spikes)])    # Subset the relevant columns of the data table
melted.spikes.all = melt(spikes.all,id.vars=c("left","right","well.y","well.x"))    # Melt the four spike-ins in a single variable

# Per spike-in, determine the barcodes with extremely low values:
# Using the data between 0.05 and 0.95 quantiles, apply t.test between the chosen barcode and rest of the barcodes
# (script for the t-test is a modified version of the script provided by Ayelet Alpert)

AllBcs = c(unique(as.character(melted.spikes.all$left)), unique(as.character(melted.spikes.all$right)))

resPvs = do.call('rbind',lapply(AllBcs, function(Bc){
  res = do.call('rbind',lapply(spikes, function(spike){
    RNA1Vals = as.matrix(melted.spikes.all[(melted.spikes.all$left == Bc | melted.spikes.all$right == Bc) & melted.spikes.all$variable == spike,"value"])
    RNA1Vals = as.matrix(RNA1Vals[(quantile(RNA1Vals,0.05,na.rm=T) < RNA1Vals) & (RNA1Vals < quantile(RNA1Vals,0.95,na.rm=T))])
    RNA2Vals = as.matrix(melted.spikes.all[(melted.spikes.all$left != Bc & melted.spikes.all$right != Bc) & melted.spikes.all$variable == spike,"value"])
    RNA2Vals = as.matrix(RNA2Vals[(quantile(RNA2Vals,0.05,na.rm=T) < RNA2Vals) & (RNA2Vals < quantile(RNA2Vals,0.95,na.rm=T))])
    return(data.frame(SP = spike, Bc = Bc, pv = t.test(RNA1Vals, RNA2Vals)[[3]]))
  }))
  return(res)
}))

# Empirically set a threshold (referring the data)

thresh = 1e-60

# Determine barcode-spike-in combination with very low read counts

flagPerSpike = resPvs[resPvs$pv <= thresh, ]    

# Set the flagged combinations as NA

for (i in 1:nrow(flagPerSpike)) {    
  datclean[[as.character(flagPerSpike[i,"SP"])]][which(datclean$left == as.character(flagPerSpike[i,"Bc"]) | datclean$right == as.character(flagPerSpike[i,"Bc"]))] = NA
}
summary(flagPerSpike)

# Remove the barcodes with bad combination with genes: L28 (DNMT3B), L44 (LIN28A), R23 (POU5F1)
# Remove L24 and L47 to prevent over-correction (since they are very low overall) 

datclean = subset(datclean, !(left %in% c("L24","L47","L28","L44")))
datclean = subset(datclean, !(right %in% c("R23")))

# Per spike-in, replace the NAs with the median of the rest

datclean$RNA1[is.na(datclean$RNA1)] = median(datclean$RNA1, na.rm=T)
datclean$RNA2[is.na(datclean$RNA2)] = median(datclean$RNA2, na.rm=T)
datclean$RNA6[is.na(datclean$RNA6)] = median(datclean$RNA6, na.rm=T)
datclean$RNA8[is.na(datclean$RNA8)] = median(datclean$RNA8, na.rm=T)


# Boxplots of the adjusted spike-ins and the scaling factor (black)

datclean$RNA_X = 2^rowMeans(log2(datclean[spikes]+1))-1

pdf(paste0(NGS,"_",lib," Spike-in BoxPlots (Adjusted).pdf"), width=10, height=3)
for (bc in c("left","right")){
  
  box.melt = melt(datclean[c(bc,spikes,"RNA_X")], id.vars=c(bc))
  box.melt$value = log2(box.melt$value+1)
  
  plot(
    ggplot(box.melt, aes_string(x=bc, y="value")) +
      geom_jitter(data=subset(box.melt, variable != "RNA_X"), alpha=0.8, size=0.3, aes_string(color="variable")) +
      geom_jitter(data=subset(box.melt, variable == "RNA_X"), alpha=1, size=0.4) +
      geom_boxplot(data=subset(box.melt, variable == "RNA_X"), outlier.shape=NA, alpha=0, lwd=0.3) +
      labs(title=paste0(bc, " spike-ins (adjusted)"), x="", y="Adjusted reads (log2)") +
      theme_classic() +
      ylim(0,NA) +
      theme(plot.title = element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
      guides(colour=guide_legend(override.aes=list(size=4)), fill=guide_legend(title="Spike-in")))
}
dev.off()


#################################################################################################################
################################   NORMALIZE THE DATA  ##########################################################
#################################################################################################################

datnormal = datclean

# Remove outlier scaling factors to prevent over-correction

min = median(datnormal$RNA_X, na.rm=T)/10
max = median(datnormal$RNA_X, na.rm=T)*10
datnormal = subset(datnormal, (min < RNA_X) & (RNA_X < max))

# To preserve gene magnitudes, center the scaling factors around 1

factor = datnormal$RNA_X / median(datnormal$RNA_X, na.rm=T)

# Normalize

datnormal[c(genes,"RNA_X")] = datnormal[c(genes,"RNA_X")] / factor


# Heatmaps of the normalized data as pdf (organized as plates)

pdf(paste0(NGS,"_",lib," Heatmaps (Normalized).pdf"), width=14,height=10)
for (gene in genes) {
  
  try(heatmap.2(log2(acast(datnormal, well.y~well.x, value.var=gene)+1), Rowv=F, Colv=F, dendrogram='none',
                trace="none", density.info="density", key.title=NA, cexRow=1.2, cexCol=1.2, adjCol=c(0.9,0.5),
                main=paste0("\n\n\n",NGS," - ",lib,"\n",gene, "\nNormalized"),
                labRow = acast(datnormal, well.y~well.x, value.var="left")[,1], 
                labCol = acast(datnormal, well.y~well.x, value.var="right")[1,],
                col=colorRampPalette(c("white","Azure","PowderBlue","SteelBlue","MidnightBlue"))(300), lwid=c(2,10), lhei=c(2,10),
                cellnote=acast(datnormal, well.y~well.x, value.var="cell.count"), notecol="gray30", notecex=0.7))
}
dev.off()


# Correlograms of the raw and normalized spike-ins

pdf(paste0(NGS,"_",lib," Correlation Plots Before-After.pdf"), width=10, height=10)
  ggpairs(spikes.all, columns = which(colnames(spikes.all) %in% spikes), title = "Raw Reads", axisLabels = "show")+theme(plot.title = element_text(hjust = 0.5))
  ggpairs(datnormal, columns = which(colnames(datnormal) %in% spikes), title = "Normalized Reads", axisLabels = "show")+theme(plot.title = element_text(hjust = 0.5))
dev.off()


# Boxplots of the normalized spike-ins and the scaling factor (black)

pdf(paste0(NGS,"_",lib," Spike-in BoxPlots (Normalized).pdf"), width=10, height=3)
for (bc in c("left","right")){
  
  box.melt = melt(datnormal[c(bc,spikes,"RNA_X")], id.vars=c(bc))
  box.melt$value = log2(box.melt$value+1)
  
  plot(
    ggplot(box.melt, aes_string(x=bc, y="value")) +
      geom_jitter(data=subset(box.melt, variable != "RNA_X"), alpha=0.8, size=0.3, aes_string(color="variable")) +
      geom_jitter(data=subset(box.melt, variable == "RNA_X"), alpha=1, size=0.4) +
      geom_boxplot(data=subset(box.melt, variable == "RNA_X"), outlier.shape=NA, alpha=0, lwd=0.3) +
      labs(title=paste0(bc, " spike-ins (normalized)"), x="", y="Normalized reads (log2)") +
      theme_classic() +
      ylim(0,NA) +
      theme(plot.title = element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
      guides(colour=guide_legend(override.aes=list(size=4)), fill=guide_legend(title="Spike-in")))
}
dev.off()


#################################################################################################################
################################   WRITE TO EXCEL THE ADJUSTED & NORMALIZED DATA  ###############################
#################################################################################################################


write.xlsx2(as.data.frame(datclean[setdiff(colnames(datclean),"RNA_X")]), paste0(NGS,"_Adjusted.xlsx"), sheetName = lib, col.names=TRUE, row.names=FALSE, append=TRUE)
write.xlsx2(as.data.frame(datnormal[setdiff(colnames(datnormal),"RNA_X")]), paste0(NGS,"_Normalized.xlsx"), sheetName = lib, col.names=TRUE, row.names=FALSE, append=TRUE)

