
#################################################################################################################################
#
# Project: BART-Seq
# Author: Fatma Uzbas
#
##################################################   SIMPLE DIVISION METHOD   ###################################################
#
# This script corrects the barcode-primer effect observed in spike-in read counts, and normalizes the data using spike-ins.
#
# Correction factors are calculated by dividing the glm predictions with interactions (mp1) to the glm predictions without
# interactions (mp0)
#
# mp2: Variable + left + right + cells + L:V + R:V
# mp1: Variable + left + right + L:V + R:V
# mp0: Variable + left + right
#
# The script has overall three main steps:
# 1 - Filters the data for outliers
# 2 - Models and corrects for the barcode-primer combination effect
# 3 - Normalizes the data
#
#################################################################################################################################


# install.packages("readxl")
# install.packages("reshape2")
# install.packages("MASS")
# install.packages("gplots")
# install.packages("corrplot")
# install.packages("GGally")
# install.packages("xlsx")
# install.packages("scales")
# install.packages("ggplot2")
# install.packages("EnvStats")


library(readxl)
library(reshape2)
library(MASS)
library(gplots)
library(corrplot)  # For correlograms
library(GGally)  # For correlograms
library(xlsx)
library(scales)
library(ggplot2)
library(EnvStats) # For Geometric Mean


WD = "D:/Users/fatma.uzbas/Desktop/GitHub"
ExcelFile = "NGS12_Lib1.xlsx"   # The read count table
NGS = "NGS12"   # Name of the experiment
lib = "Lib1"   # Name of the library & tab in excel file


genes = c("B2M", "CCND1", "CCNE1", "CER1", "DNMT3B", "GAPDH", "LIN28A", "NANOG", "POU5F1", "SOX2", "ZFP42", "RNA1", "RNA2", "RNA6", "RNA8")
spikes = c("RNA1","RNA2","RNA6", "RNA8") 
setwd(WD)


# The original data

dat = read_excel(ExcelFile, sheet = lib, col_names = TRUE, col_types = NULL, na = "", skip = 0)


#################################################################################################################################
###############################################   FILTER THE DATA FOR OUTLIERS   ################################################
#################################################################################################################################



##############################################################################################
#-----------------   Clean up the combinations with very low read counts   ------------------#
##############################################################################################


# Set the cell count of bulk RNA wells as "bulk"

dat2 = dat
dat2$cells[which(dat2$type == "bulk RNA")] = "bulk"



# Take the spike-in part of the original data

spikes.all = data.frame(dat2[c("left","right","well.y","well.x","cell.count","cells",spikes)])    # Subset the relevant columns of the data table
melted.spikes.all = melt(spikes.all,id.vars=c("left","right","well.y","well.x","cell.count","cells"))    # Melt the four spke-ins in a single variable



# Per RNA, determine rows/columns with extremely low values
# First take out the outliers (use data between 0.05 and 0.95 quantiles) and then apply t.test between the chosen barcode and rest of the barcodes
# This is done without scaling because it is checking very low read counts only (not the bad combination effects)
  
AllBcs = c(unique(as.character(melted.spikes.all$left)), unique(as.character(melted.spikes.all$right)))

resPvs = do.call('rbind',lapply(AllBcs, function(Bc){
  res = do.call('rbind',lapply(spikes, function(spike){
    RNA1Vals = as.matrix(melted.spikes.all[(melted.spikes.all$left == Bc | melted.spikes.all$right == Bc) & melted.spikes.all$variable == spike,"value"])
    RNA1Vals = as.matrix(RNA1Vals[(quantile(RNA1Vals,0.05) < RNA1Vals) & (RNA1Vals < quantile(RNA1Vals,0.95))])
    RNA2Vals = as.matrix(melted.spikes.all[(melted.spikes.all$left != Bc & melted.spikes.all$right != Bc) & melted.spikes.all$variable == spike,"value"])
    RNA2Vals = as.matrix(RNA2Vals[(quantile(RNA2Vals,0.05) < RNA2Vals) & (RNA2Vals < quantile(RNA2Vals,0.95))])
    return(data.frame(SP = spike, Bc = Bc, pv = t.test(RNA1Vals, RNA2Vals)[[3]]))
  }))
  return(res)
}))



thresh = 1E-200   # Set a threshold after checking the resPvs and the raw data
flagPerSpike = resPvs[resPvs$pv <= thresh, ]   # Determine row/column and spike-in combination with very low read counts



# Filter out the spike-in-row/column combination with very low reads

datclean = dat    # To create a cleaned file, first take the raw data

for (i in 1:nrow(flagPerSpike)) {    # Set the very bad (flagged - low read) combinations as NA 
  datclean[[as.character(flagPerSpike[i,"SP"])]][which(datclean$left == as.character(flagPerSpike[i,"Bc"]) | datclean$right == as.character(flagPerSpike[i,"Bc"]))] = NA
}

datclean$nas = apply(is.na(datclean[spikes]), 1, sum)    # Count how many NAs are there per well
datclean = subset(datclean, nas < 2)[,which(names(datclean) != "nas")]    # Remove the row, where 2+ spike-ins are NAs



##############################################################################################
#--------------------------------   Clean up the Outliers    --------------------------------#
##############################################################################################


# Take the spike-in part of the cleaned data

spikes.clean = data.frame(datclean[c("left","right","well.y","well.x","cell.count","cells",spikes)])
melted.spikes.clean = melt(spikes.clean,id.vars=c("left","right","well.y","well.x","cell.count","cells"), na.rm=T)    # Do not take the NAs into the melted data

melt.plate.0 = melted.spikes.clean



# Determine the outliers that do not fit into the full model

mp2 = glm.nb(value ~ variable + left + right + cells + left:variable + right:variable, data = melt.plate.0)    # Full model


melt.plate.0$res = abs(residuals(mp2))   # Absolute value of residuals
melt.plate.0$rem = "no"   # Should this row be removed? First set everything to "no"
melt.plate.0$rem[melt.plate.0$res > quantile(melt.plate.0$res, c(0.95))] = "yes"   # Label the outliers (95th quantile of residuals) that should be removed: "yes"



# Check how does filtering out outliers improve the model 

melt.plate.f = subset(melt.plate.0, rem == "no")    # Without outliers
mp3 = glm.nb(value ~ variable + left + right + cells + left:variable + right:variable, data = melt.plate.f)


pdf(paste0(NGS,"_",lib," Residuals mp2 (Full).pdf"), width=10,height=10)    # Compare the model before removing outliers...
plot(mp2)
dev.off()

pdf(paste0(NGS,"_",lib," Residuals mp3 (Filtered).pdf"), width=10,height=10)    # ...with the model after removing outliers
plot(mp3)
dev.off()



# Filter the data once more, this time for outliers, for subsequent steps

datclean.2 = datclean

additionalFlags = subset(melt.plate.0, rem == "yes")[c("left","right","variable")]    # Flags for the outliers

for (i in 1:nrow(additionalFlags)) {    # Set the outliers as NA
datclean.2[[as.character(additionalFlags$variable[i])]][which((datclean.2$left == additionalFlags$left[i]) & (datclean.2$right == additionalFlags$right[i]))] = NA
}


# Filter the wells that has more than one NAs

datclean.2$nas = apply(is.na(datclean.2[spikes]), 1, sum)    # Count how many NAs are there per well
datclean.2 = subset(datclean.2, nas < 2)[,which(names(datclean.2) != "nas")]    # Remove the row, where 2+ spike-ins are NAs




#################################################################################################################################
###################################   MODEL & CORRECT FOR BARCODE-PRIMER COMBINATION EFFECT   ###################################
#################################################################################################################################



# Take the spike-in part of the cleaned.2 data

spikes.clean.2 = data.frame(datclean.2[c("left","right","well.y","well.x","cell.count","cells",spikes)])


# Scale each spike-in to the same range, so that, when taking the geometric mean, the number/identity of spike-ins in the well doesn't matter
# when one of them is NA

minC = ifelse(min(spikes.clean.2[spikes], na.rm=T) == 0, 1, min(spikes.clean.2[spikes], na.rm=T))   # If min is 0, replace it with 1
maxC = max(spikes.clean.2[spikes], na.rm=T)

spikes.clean.2[spikes] = sapply(spikes.clean.2[spikes], function(x) rescale(x, to = c(minC,maxC)))

melted.spikes.clean.2 = melt(spikes.clean.2,id.vars=c("left","right","well.y","well.x","cell.count","cells"), na.rm=T)
melted.spikes.clean.2$value = round(melted.spikes.clean.2$value, 0)    # Since glm.nb requires integers, round the values



# Calculate the correction factor

mp1 = glm.nb(value ~ variable + left + right + left:variable + right:variable, data = melted.spikes.clean.2)   # Model with interactions
mp0 = glm.nb(value ~ variable + left + right, data = melted.spikes.clean.2)   # Model without interactions


melt.plate = melted.spikes.clean.2

melt.plate$prd3 = predict.glm(mp3, melt.plate, type="response")   # Full predictions (filtered model)
melt.plate$prd1 = predict.glm(mp1, melt.plate, type="response")   # With interactions
melt.plate$prd0 = predict.glm(mp0, melt.plate, type="response")   # Without interactions


melt.plate$fold = melt.plate$prd1/melt.plate$prd0   # Correction factor: fold-change of the interactions relative to the average (division)



# Correct the spike-in read counts

melt.plate$corF = melt.plate$value/melt.plate$fold



# Save and re-load the workspace image

# save.image(file=paste0(NGS,"_",lib," AfterFiltering.Rdata"))
# load(paste0(NGS,"_",lib,"AfterFiltering.Rdata"))


##############################################################################################
#--------------------------   Heatmaps of the correction steps   ----------------------------# 
##############################################################################################



# Heatmaps of the correction steps

steps = c("value", "prd3", "prd1", "prd0", "corF")
steplabels = c("Raw", "(filtered) Model Predictions", "with Interactions", "without Interactions (Average)",  "Corrected Read Counts")

pdf(paste0(NGS,"_",lib," Spike-in Correction Steps.pdf"), width=14,height=10)
for (gene in spikes) {
  
  for (step in 1:5) {
    
    sp = subset(melt.plate, variable==gene)
    frame = dcast(sp, left + right + well.y + well.x + cell.count ~ variable, value.var=steps[[step]])
    try(heatmap.2(log2(acast(frame, as.character(well.y)~well.x, value.var=gene)+1), Rowv=F, Colv=F, dendrogram='none',
                  trace="none", density.info="density", key.title=NA, cexRow=1.4, cexCol=1.4, adjCol=c(NA,0.5),
                  main=paste0("\n\n\n\n",NGS," - ",lib,"\n\n",gene, "\n\n", steplabels[[step]]),
                  labRow = acast(frame, as.character(well.y)~well.x, value.var="left")[,12],
                  labCol = acast(frame, as.character(well.y)~well.x, value.var="right")[8,],
                  col=colorRampPalette(c("white","Azure","PowderBlue","SteelBlue","MidnightBlue"))(300), lwid=c(3,10), lhei=c(2,7),
                  cellnote=acast(frame, as.character(well.y)~well.x, value.var="cell.count"), notecol="gray30", notecex=1.2))
  }
}
dev.off()



# Data frame of corrected spike-ins

corrected = dcast(melt.plate, left + right + well.y + well.x + cell.count + cells ~ variable, value.var="corF")



##############################################################################################
#----------------   Correlograms of spike-ins before & after correction   -------------------# 
##############################################################################################


pdf(paste0(NGS,"_",lib," Correlation Plots (Raw).pdf"), width=6, height=6)
corrplot(cor(spikes.all[spikes],spikes.all[spikes]), type="upper", order="hclust", tl.col="black", tl.srt=45, mar=c(0,4,0,4), title="\n\n\n\n\n\nRaw Reads")
ggpairs(spikes.all, columns = which(colnames(spikes.all) %in% spikes), title = "Raw Reads", axisLabels = "show",
        upper=list(continuous = wrap("cor",size=4,alignPercent=1,color="steelblue")), lower = list(continuous = wrap("points",size=0.5))) +
          theme(plot.title=element_text(hjust=0.5), axis.text = element_text(size=9), axis.text.x = element_text(angle=45, hjust=1))
dev.off()

pdf(paste0(NGS,"_",lib," Correlation Plots (Simple Division).pdf"), width=6, height=6)
corrplot(cor(corrected[spikes],corrected[spikes], use = "pairwise.complete.obs"), type="upper", order="hclust", tl.col="black", tl.srt=45, mar=c(0,4,0,4), title="\n\n\n\n\n\nCorrected Reads (Simple Division)")
ggpairs(corrected, columns = which(colnames(corrected) %in% spikes), title = "Corrected Reads", axisLabels = "show",
        upper=list(continuous = wrap("cor",size=4,alignPercent=1,color="steelblue")), lower = list(continuous = wrap("points",size=0.5))) +
  theme(plot.title=element_text(hjust=0.5), axis.text = element_text(size=9), axis.text.x = element_text(angle=45, hjust=1))
dev.off()



#################################################################################################################################
###################################################   NORMALIZE   THE   DATA  ###################################################
#################################################################################################################################



# Calculate the normalization factor per well using geometric mean (all four spike-ins per well, or three if one is NA)

corrected2 = corrected    # Copy the corrected spike-ins
corrected2$factor = 0    # Initiate normalization factor column

corspikes = paste0(spikes, "cor")    # Labels for corrected spike-ins

for(i in 1:nrow(corrected2)){    # Calculate the normalization factor
  corrected2$factor[i] = geoMean(as.numeric(corrected2[c("RNA1","RNA2","RNA6")][i,]), na.rm = T)    # In this specific example RNA8 is excluded due to low read counts
}


corrected2$factor = corrected2$factor / mean(corrected2$factor)   # Center


# Normalize the data (cleaned.2)

normalized = datclean.2
normalized = merge(normalized, corrected, by = c("left","right","well.y","well.x","cell.count","cells"), suffixes = c("","cor"))   # Add the corrected spike-ins to the data
normalized$factor = corrected2$factor    # Add the factors

normalized2 = normalized
normalized2[c(genes, corspikes)] = normalized2[c(genes, corspikes)] / normalized2$factor   # Normalize



##############################################################################################
#------------------------------   Heatmaps of Normalized Data  ------------------------------#
##############################################################################################


# Heatmaps of the normalized data

pdf(paste0(NGS,"_",lib," Heatmaps Normalized (Simple Division).pdf"), width=14,height=10)
for (gene in c(genes, corspikes)) {
  
  frame = normalized2[c("well.y", "well.x", "left", "right", "cell.count", gene)]
  try(heatmap.2(log2(acast(frame, as.character(well.y)~well.x, value.var=gene)+1), Rowv=F, Colv=F, dendrogram='none',
                trace="none", density.info="density", key.title=NA, cexRow=1.4, cexCol=1.4, adjCol=c(NA,0.5),
                main=paste0("\n\n\n\n",NGS," - ",lib,"\n\n",gene, "\n\nNormalized (Simple Division)"),
                labRow = acast(frame, as.character(well.y)~well.x, value.var="left")[,8],
                labCol = acast(frame, as.character(well.y)~well.x, value.var="right")[12,],
                col=colorRampPalette(c("white","Azure","PowderBlue","SteelBlue","MidnightBlue"))(300), lwid=c(3,10), lhei=c(2,7),
                cellnote=acast(frame, as.character(well.y)~well.x, value.var="cell.count"), notecol="gray30", notecex=1.2))
}
dev.off()



#################################################################################################################################
#######################################################   WRITE TO EXCEL  #######################################################
#################################################################################################################################


write.xlsx2(as.data.frame(normalized2), paste0(NGS," Normalized Simple Method.xlsx"), sheetName = lib, col.names=TRUE, row.names=FALSE, append=TRUE)

