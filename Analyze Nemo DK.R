##### DATA NOTES #####

#' Data from Dom 1
#' 
#' Channels:
#' 1 = BF
#' 2 = p65
#' 5 = DRAQ5
#' 6 = eFLUOR780
#' 
#' Prefiltering:
#' Gradient RMS BF? (maybe DRAQ5)
#' eFLUOR+


rm(list=ls())
gc()
source('~/GitHub/R-General/.Rprofile')
# source('~/GitHub/R-Informatics-private/HuMoments.R')
source('~/GitHub/R-Informatics-private/Zernike.R')
source('~/GitHub/R-Cytoprofiling/PreProcessingHelpers.R')
# source('~/Documents/GitHub/R-General/.Rprofile')
# source('~/Documents/GitHub/R-Informatics/HuMoments.R')
# source('~/Documents/GitHub/R-Informatics-private/Zernike.R')
# source('~/Documents/GitHub/R-Cytoprofiling/PreProcessingHelpers.R')
# source('~/Documents/GitHub/R-Informatics-private/machlearn.R')
library(bit64)
library(rgl)
library(data.table)



##### READ FEATURE DATA #####
jData <- list()
jData[['0']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=0, type='File', name='Feature CSV Table', list(Mouse='WT',Rep='1',Valid='true'))
jData[['1']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=1, type='File', name='Feature CSV Table', list(Mouse='DK',Rep='1',Valid='true'))
jData[['2']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=2, type='File', name='Feature CSV Table', list(Mouse='NES',Rep='1',Valid='true'))
jData[['3']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=3, type='File', name='Feature CSV Table', list(Mouse='Naive',Rep='1',Valid='true'))
jData[['4']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=1, y=0, type='File', name='Feature CSV Table', list(Mouse='WT',Rep='2',Valid='true'))
jData[['5']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=1, y=1, type='File', name='Feature CSV Table', list(Mouse='DK',Rep='2',Valid='true'))
jData[['6']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=1, y=2, type='File', name='Feature CSV Table', list(Mouse='NES',Rep='2',Valid='true'))
jData[['8']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=2, y=0, type='File', name='Feature CSV Table', list(Mouse='WT',Rep='3',Valid='true'))
jData[['9']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=2, y=1, type='File', name='Feature CSV Table', list(Mouse='DK',Rep='3',Valid='true'))
jData[['10']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=2, y=2, type='File', name='Feature CSV Table', list(Mouse='NES',Rep='3',Valid='true'))
jData[['12']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=3, y=0, type='File', name='Feature CSV Table', list(Mouse='WT',Rep='4',Valid='true'))
jData[['13']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=3, y=1, type='File', name='Feature CSV Table', list(Mouse='DK',Rep='4',Valid='true'))
jData <- rbindlist(jData)
x <- readJEXDataTables(jData, idCols=c('Id','Page'))

# Remove unnecessary measures
x <- removeMeasurementNamesContaining(x, c('Symmetry','ImageMoments','Circle'))

##### READ COLOCALIZATION DATA #####
jData <- list()
jData[['0']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=0, type='File', name='Coloc CSV Table', list(Mouse='WT',Rep='1',Valid='true'))
jData[['1']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=1, type='File', name='Coloc CSV Table', list(Mouse='DK',Rep='1',Valid='true'))
jData[['2']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=2, type='File', name='Coloc CSV Table', list(Mouse='NES',Rep='1',Valid='true'))
jData[['3']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=0, y=3, type='File', name='Coloc CSV Table', list(Mouse='Naive',Rep='1',Valid='true'))
jData[['4']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=1, y=0, type='File', name='Coloc CSV Table', list(Mouse='WT',Rep='2',Valid='true'))
jData[['5']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=1, y=1, type='File', name='Coloc CSV Table', list(Mouse='DK',Rep='2',Valid='true'))
jData[['6']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=1, y=2, type='File', name='Coloc CSV Table', list(Mouse='NES',Rep='2',Valid='true'))
jData[['8']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=2, y=0, type='File', name='Coloc CSV Table', list(Mouse='WT',Rep='3',Valid='true'))
jData[['9']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=2, y=1, type='File', name='Coloc CSV Table', list(Mouse='DK',Rep='3',Valid='true'))
jData[['10']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=2, y=2, type='File', name='Coloc CSV Table', list(Mouse='NES',Rep='3',Valid='true'))
jData[['12']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=3, y=0, type='File', name='Coloc CSV Table', list(Mouse='WT',Rep='4',Valid='true'))
jData[['13']] <- readJEXData(dbPath='Y:/Jay/A JEX Databases/20180810 - Nemo DK', ds='Plate', x=3, y=1, type='File', name='Coloc CSV Table', list(Mouse='DK',Rep='4',Valid='true'))
jData <- rbindlist(jData)
y <- readJEXDataTables(jData, idCols=c('Page','Id'), samples.to.match.and.append=x)

##### GATHER/PREPARE SOME IMPORTANT COLUMN NAMES #####

extras <- getAllColNamesExcept(y, c('cId','ImageChannel','MaskChannel','Measurement','Value'))
idCols <- c('x','y','Page','Id')

makeComplexId(y, cols=idCols)

x1 <- copy(y)#x[cId %in% sample(unique(cId), size=3000)]

# Avoid standardizing the id cols
lapply.data.table(x1, as.character, cols=c('x','y','Page','Id','ImageChannel','MaskChannel'), in.place=T)

# Get rid of extraneous text in the measurement names that are carried over from JEX/Java
fixLongTableStringsInCol(x1, 'Measurement')

##### CONVERT TO WIDE TABLE WITH IMAGECHANNEL AND MASKCHANNEL #####

# First reorganize, keeping ImageChannel and MaskChannel as columns for various calculations
x2 <- reorganize(x1, idCols=c('cId','ImageChannel','MaskChannel',extras), measurementCols=c('Measurement'), valueCols=c('Value'), sep='.')

##### PRE-STANDARDIZATION CALCULATIONS #####

# Transform the correlation coefficients if necessary (not SymmetryCorrelation yet. Still need to use for other calcs.)
correlationNames <- getColNamesContaining(x2,'CorrelationCoefficient')
correlationNames <- correlationNames[!grepl('Similarity',correlationNames)]
if(length(correlationNames) > 0)
{
	lapply.data.table(x2, FUN=sim.transform, cols=correlationNames, in.place=T)
	setnames(x2, correlationNames, paste0(correlationNames, '.Similarity'))
}

# Summarize the geometric data that currently has data for each subregion of each mask.
x2 <- summarizeGeometry(x2, cellIdCols=c('cId', extras), removeXY=T) # Expect a warning from duplication of values (these duplicates are then removed)

# Normalize the Zernike mangnitudes to make them intensity invariant
x2 <- meanNormalizeWideZernikeMoments(x2, meanName='Stats.Mean')

# Copy data at this point
x3 <- copy(x2)

##### REORGANIZE #####

# For any cols that are factors, refactor to remove non-existent possible values
x3 <- refactor(x3)

# Reorganize the table
x3 <- reorganize(x3, idCols=c('cId',extras), measurementCols=c('ImageChannel','MaskChannel'), valueCols=getAllColNamesExcept(x3, c('cId','ImageChannel','MaskChannel',extras)), sep='.')

# Remove all columns that are ALL NON-FINITE
removeColsMatching(x3, col.test=function(a){all(!is.finite(a))})

##### RELABEL STUFF #####

# Change values of conditions for readability.
x3[, Mouse:=factor(Mouse, levels=c('WT','DK','NES','Naive'))]
x3 <- x3[order(Mouse, Rep)]

##### CALCULATE SOME METRICS TO PLOT #####

x4 <- copy(x3)
fileToSave <- 'Y:/Jay/R Projects/20180810 - Nemo DK/x3_All.csv'
# fileToSave <- '/Volumes/Miyamoto/Jay/R Projects/Dom1/x4.csv'
fwrite(x3, file=fileToSave)
# x4 <- fread(fileToSave, header=T)

# Plot things for gating
# x4[, Focus:=log((Stats.Sum.NSR.NMem/Stats.Sum.NSR.Nuc)/(Geometric.SizeIterable.None.NMem/Geometric.SizeIterable.None.WholeCell))]
# data.table.plot.all(x4[sample(.N, 10000)], type='p', alpha=0.1, xcol='Focus', ycol='Geometric.Roundness.None.WholeCell')
data.table.plot.all(x4[sample(.N, 10000)], type='p', alpha=0.1, xcol='ZernikeMag20_THISwSEC.NSR.Nuc', ycol='Geometric.Roundness.None.WholeCell')

# Take the top 95% in mean Nuc staining intensity (i.e., needs to stain robustly for nucleus)
# Take the top 85% in size
# Take the top 66% in roundness
# Take CD8+ and CD44+
# Take the top 300 in focus
threshStain <- getPercentileValues(x4$Stats.Mean.2.WholeCell, 0.1)
threshNuc <- getPercentileValues(x4$Stats.Mean.5.Nuc, 0.05)
threshSize <- getPercentileValues(x4$Geometric.SizeIterable.None.WholeCell, 0.15)
threshRound <- getPercentileValues(x4$Geometric.Roundness.None.WholeCell, 0.33)
threshNSR <- getPercentileValues(x4$ZernikeMag20_THISwSEC.NSR.Nuc, 0.66)
temp <- x4[Stats.Sum.3.Mem > 0 & Stats.Sum.4.Mem > 0]
clust <- assignToClustersXY(xData=logicle(temp$Stats.Sum.3.Mem, base=exp(1)), yData=logicle(temp$Stats.Sum.4.Mem, base=exp(1)))
plotClustersN(xdata=logicle(temp$Stats.Sum.3.Mem, base=exp(1)), ydata=logicle(temp$Stats.Sum.4.Mem, base=exp(1)), xlab='CD44', ylab='CD8', clusterNums=clust$clusters)
thresh3_CD44 <- exp(clust$threshX)
thresh4_CD8 <- exp(clust$threshY)
abline(h=log(thresh4_CD8))
abline(v=log(thresh3_CD44))
x4[, CD8:=Stats.Sum.4.Mem > thresh4_CD8]
x4[, CD44:=Stats.Sum.3.Mem > thresh3_CD44]
x5 <- x4[Stats.Mean.5.Nuc > threshNuc & CD8==T & Geometric.SizeIterable.None.WholeCell > threshSize & Geometric.Roundness.None.WholeCell > threshRound]
x5[, list(n=.N), by=.(Mouse,Rep)]
setorder(x5, -ZernikeMag20_THISwSEC.NSR.Nuc)
x5 <- x5[, .SD[1:200], by=.(Mouse,Rep)]
hist(x4$ZernikeMag20_THISwSEC.NSR.Nuc, col=adjustColor('red',0.5), xlab='yo', ylab='yo')
hist(x5$ZernikeMag20_THISwSEC.NSR.Nuc, col=adjustColor('black', 0.2), add=T)
x5[, list(n=.N), by=.(Mouse,Rep)]

data.table.plot.all(x4, type='d', xcol='Stats.Mean.2.WholeCell')
# 
# data.table.plot.all(x5, type='p', alpha=0.1, xcol='Stats.Mean.5.Nuc', ycol='ZernikeMag20_THISwSEC.NSR.Nuc')
# 
# data.table.plot.all(x4, type='p', alpha=0.1, xcol='Stats.Mean.NSR.NMem', ycol='Stats.Mean.5.NMem')
# data.table.plot.all(x4, type='p', alpha=0.1, xcol='ZernikeMag20_THISwSEC.NSR.Nuc', ycol='Stats.Mean.5.Nuc')
x4 <- copy(x5)

# Calculate metrics
#/ZernikeMag20_THISwFIXED.NSRp6Nuc.Nuc)]
# Check how many are not finite for RL
sum(!is.finite(x4$RL))
# x4[, Log.Inv.RL:=log(1/RL)]#/ZernikeMag20_THISwFIXED.NSRp6Nuc.Nuc)]
# x4[, PCC.Nuc:=Stats.PearsonsCorrelationCoefficient.Similarity.2_5.Nuc]
x4[, PCC.WC:=Stats.PearsonsCorrelationCoefficient.Similarity.2_5.WholeCell]
x4[, NormSum:=Stats.Sum.2.WholeCell/median(Stats.Sum.2.WholeCell[Mouse=='WT'])]
x4[, CellVolume:=(4*pi/3)*(Geometric.SizeIterable.None.WholeCell/pi)^(3/2)/10000]
x4[, NucVolume:=(4*pi/3)*(Geometric.SizeIterable.None.Nuc/pi)^(3/2)/10000]
x4[, NucVolFrac:=NucVolume/CellVolume]
x4[, NucRadFrac:=(NucVolume/CellVolume)^(1/3)]
x4[, sSampled:=3]
x4[, sActual:=NucVolume/CellVolume]
# x4[, CRatio:=calcCRatio(RL, 1, getPercentileValues(RL, levels=0.99), sSampled), by=c('Stain','Stim')]
# x4[, AFrac:=calcARatio(RL, 1, getPercentileValues(RL, levels=0.99), sSampled, sActual), by=c('Stain','Stim')]
# x4[, lCRatio:=log2(CRatio)]
# x4[, lAFrac:=log2(AFrac)]

# # Calculate normalize radial localization, standardizing per ds then rescaling by means of the medians and mads of all the ds's
# x4[, SRL:=(Log.Inv.RL-median(Log.Inv.RL[Mouse=='WT']))/mad(Log.Inv.RL[Mouse=='WT']), by=c('Stim','ds','Stain')]
# duh <- x4[, list(med=median(Log.Inv.RL[Mouse=='WT']), mad=mad(Log.Inv.RL[Mouse=='WT'])), by=c('Stim','ds','Stain')]
# duh <- duh[, list(med=mean(med), mad=mean(mad)), by=c('Stim','Stain')]
# x4[, NRL:=SRL*duh$mad[duh$Stim==.BY[[1]] & duh$Stain==.BY[[2]]]+duh$med[duh$Stim==.BY[[1]] & duh$Stain==.BY[[2]]], by=c('Stim','Stain')]
# x4[, SRL2:=(Log.Inv.RL-median(Log.Inv.RL[Mouse=='WT']))/mad(Log.Inv.RL[Mouse=='WT']), by=c('Stim','ds','Stain')]
# x4 <- x4[order(Stim2, ds, Stain, -Mouse, LMB)]


# # Calculate normalize RADIAL LOCALIZATION, standardizing per ds then rescaling by means of the medians and mads of all the ds's
x4[, RL:=RadialLocalization.2_5.Nuc]
x4[, NRL:=(RL)/median(RL[Mouse=='WT'])]
# duh <- x4[, list(med=median(RL[Mouse=='WT'])), by=c('Rep')]
# duh <- duh[, list(med=mean(med, na.rm=T))]
# x4[, NRL:=SRL*duh$med]
# x4 <- x4[order(Rep, -Mouse)]
# 
# # Calculate normalize SUM INTENSITY, standardizing per ds then rescaling by means of the medians and mads of all the ds's
x4[, NucAmount:=Stats.Sum.2.Nuc]
x4[, NNucAmount:=NucAmount/median(NucAmount[Mouse=='WT'])]
# duh <- x4[, list(med=median(NucAmount[Mouse=='WT'])), by=c('Rep')]
# duh <- duh[, list(med=mean(med, na.rm=T))]
# x4[, NNucAmount:=SNucAmount*duh$med]
# x4 <- x4[order(Rep, -Mouse)]
# 
# # Calculate normalize SUM INTENSITY, standardizing per ds then rescaling by means of the medians and mads of all the ds's
x4[, WCAmount:=Stats.Sum.2.WholeCell]
x4[, NWCAmount:=WCAmount/median(WCAmount[Mouse=='WT'])]
# duh <- x4[, list(med=median(WCAmount[Mouse=='WT'])), by=c('Rep')]
# duh <- duh[, list(med=mean(med, na.rm=T))]
# x4[, NWCAmount:=SWCAmount*duh$med]
# x4 <- x4[order(Rep, -Mouse)]
# 
# # Calculate normalize MEAN INTENSITY, standardizing per ds then rescaling by means of the medians and mads of all the ds's
x4[, NucMean:=Stats.Mean.2.Nuc]
x4[, NNucMean:=NucMean/median(NucMean[Mouse=='WT'])]
# duh <- x4[, list(med=median(NucMean[Mouse=='WT'])), by=c('Rep')]
# duh <- duh[, list(med=mean(med, na.rm=T))]
# x4[, NNucMean:=SNucMean*duh$med]
# x4 <- x4[order(Rep, -Mouse)]

fileToSave <- 'Y:/Jay/R Projects/Dom1/x4_300.csv'
# fileToSave <- '/Volumes/Miyamoto/Jay/R Projects/Dom1/x4.csv'
fwrite(x4, file=fileToSave)
# x4 <- fread(fileToSave, header=T)


##### PLOT STUFF #####

# Additional Calcs
z <- copy(x4)
z[, CRatio:=calcCRatio2(RL=NRL, RLMin=getPercentileValues(z$NRL, levels=c(0.01)), RLMax=getPercentileValues(z$NRL, levels=c(0.99)), s=sSampled)]
z[, AFrac:=calcARatio2(NRL, 1, getPercentileValues(NRL[Mouse=='WT'], levels=0.99), sSampled, sActual)]
z[, lCRatio:=logit.transform(CRatio)]
z[, lAFrac:=logit.transform(AFrac)]
z[, Temp:=ZernikeMag20_THISwSEC.2.Nuc]
# z[, NNucAmount:=exp(NlNucAmount)/median(exp(NlNucAmount[Mouse=='WT']), na.rm=T)]
# z[, NWCAmount:=exp(NlWCAmount)/median(exp(NlWCAmount[Mouse=='WT']), na.rm=T)]

# Labels for plots
z[, Mouse:=factor(Mouse, levels=c('Naive','WT','DK','NES'))]
setorder(z, Mouse, -CD44)
z[, Mouse:=as.character(Mouse)]
# setnames(z, c('Stats.Sum.2.Nuc','Stats.Mean.2.Nuc','Stats.Sum.2.WholeCell','Stats.Mean.2.WholeCell',''))

par(mfrow=c(1,4))
ycol='Stats.Sum.2.Nuc'#'ZernikeMag20_THISwSEC.2.Nuc'#
xcol='Geometric.SizeIterable.None.WholeCell'#'Stats.Mean.2.Nuc'#
xlims <- getPercentileValues(z[[xcol]], c(0.01,0.99))
ylims <- getPercentileValues(z[[ycol]], levels=c(0.02,0.99))
logXY='y'
if(ylims[1] < 0 && (logXY=='y' | logXY=='xy')){ylims[1] <- 1}
percentiles <- c(0.01,0.99,0.01,0.99)
data.table.plot.all(z[Mouse %in% c('Naive','WT','DK','NES') & NRL > 0 & NRL < 2], log=logXY, percentile.limits = percentiles, xcol=xcol, trans.logit=c(F,F), contour.adj=c(1.5,1.5), ylim=ylims, xlim=xlims, ycol=ycol, type='c', contour.ngrid = 40, contour.quantiles = T, by='CD44', plot.by='Mouse')

par(mfrow=c(1,4))
data.table.plot.all(z, xcol='Geometric.SizeIterable.None.WholeCell', type='d', percentile.limits = c(0.01,0.99,0.01,0.99), by='CD44', plot.by='Mouse')


z[, lNNucAmount:=log(NNucAmount)]
# data.table.plot.all(z[Mouse %in% c('WT','DK')], xcol='Geometric.SizeIterable.None.WholeCell', ycol='NNucAmount', contour.quantiles=T, contour.ngrid=40, type='c', percentile.limits = c(0,1,0,1), by='Mouse')
data.table.plot.all(z[Mouse %in% c('WT','DK')], xcol='NNucAmount', type='d', percentile.limits = c(0.01,0.99,0.01,0.99), by='Mouse')
data.table.plot.all(z[Mouse %in% c('WT','DK')], xcol='NNucAmount', type='c', percentile.limits = c(0.01,0.99,0.01,0.99), by='Mouse', plot.by='Rep')
data.table.plot.all(z[Mouse %in% c('WT','DK')], xcol='Stats.Sum.2.WholeCell', type='c', percentile.limits = c(0.01,0.99,0.01,0.99), by='Mouse', plot.by='Rep')


# data.table.plot.all(z[Mouse %in% c('WT','DK') & Stats.Sum.2.Nuc > 5], xcol='Geometric.SizeIterable.None.WholeCell', ycol='NNucAmount', contour.quantiles=T, contour.ngrid=40, type='c', percentile.limits = c(0,1,0,1), by='Mouse')
# data.table.plot.all(z[Mouse %in% c('WT','DK')], xcol='NRL', ycol='Geometric.SizeIterable.None.WholeCell', contour.quantiles=T, contour.ngrid=40, type='c', percentile.limits = c(0,1,0,1), by='Mouse')
# data.table.plot.all(z[Mouse %in% c('WT','DK') & NRL > 0], xcol='lCRatio', type='d', percentile.limits = c(0,1,0,1), by='Mouse')
data.table.plot.all(z[Mouse %in% c('WT','DK') & NRL > 0 & NRL < 2], xcol='Stats.PearsonsCorrelationCoefficient.Similarity.2_5.WholeCell', type='d', percentile.limits = c(0.01,0.99,0.01,0.99), by='Mouse', plot.by='Rep')
data.table.plot.all(z[Mouse %in% c('WT','DK') & NRL > 0 & NRL < 2], xcol='NRL', type='d', percentile.limits = c(0.01,0.95,0.01,0.95), by='Mouse')


# Plot Cell Area vs NNucAmount
fileToSave <- 'Y:/Jay/R Projects/Dom1/Nuc Amount vs Localization - ZoomB - HiRes'
# fileToSave <- '/Volumes/Miyamoto/Jay/R Projects/Dom1/Nuc Amount vs Localization - Zoom.pdf'
scale <- 1.5
nlevels <- 4 #c(0.1,0.20,0.40,0.60,0.80,1.00)
ngrid <- 40
dev.off()
png(file=paste0(fileToSave, ' - ', max(c(nlevels,length(nlevels))), '.png'), width=6*scale, height=scale*4, res=300, units='in', family='TT Arial')
# pdf(file=fileToSave, width=6*scale, height=scale*4, family='sans')
# dev.off()
par(mfcol=c(2,3))
stains <- c('IkBa','RelA','cRel')
stims <- c(T,F)
ymins <- c(0.4, 0.4, 0.2, 0.4, 0.5, 0.5)
ymaxs <- c(4.2, 1.9, 2.1, 1.9, 10.33, 1.9)
xmins <- c(300, 300, 300, 300, 300, 300)
xmaxs <- c(875, 495, 810, 495, 920, 495)
# ylims <- data.table(lo=)
n <- 1
i <- 1
for(stain in stains[1])
{
	j <- 1
	for(stim in stims[1])
	{
		par(mfg=c(j,i))
		data.table.plot.all(z[ds %in% DS[1:3] & Stain==stain & Stim==stim & LMB==F], contour.quantiles=T, add=F, cex.lab=1.5, main.show=F, mar=c(3.6, 3.6, 1.5, 1.5), oma=c(0,0,0,0), ycol='NNucAmount', yaxs='i', xaxs='i', type='c', las=1, legend=F, h=1, h.lty=2, h.lwd=1, h.col='black', xlim=c(300,920), ylim=c(0.2, 10.3), contour.levels=nlevels, contour.ngrid=ngrid, mgp=c(2.4,1,0), xcol='Geometric.SizeIterable.None.WholeCell', colors=c('firebrick1','gray31'), alpha=0.2, xlab='Cell Area [pixels]', by=c('Mouse'), plot.by=c('Stain','Stim'), cross.lwd=3, cross.fun=mean, cross.args=list(na.rm=T), cross.plot=F, ylab='Norm. Nuc. Expression', log='', logicle.params=list(transX=0.0000001, transY=0.01, base=2), legend.cex=1.0)
		axis(2, at=1, las=1)
		if(!(stain == 'cRel' && stim==T))
		{
			require(TeachingDemos)
			w <- grconvertX(par('usr')[1:2], from='user', to='in')
			w <- w[2]-w[1]
			h <- grconvertY(par('usr')[3:4], from='user', to='in')
			h <- h[2]-h[1]
			subplot({
				par(mfg=c(j,i))
				data.table.plot.all(z[ds %in% DS[1:3] & Stain==stain & Stim==stim & LMB==F], contour.quantiles=T, cex.lab=1.5, add=F, mar=NULL, main.show=F, ycol='NNucAmount', yaxs='i', xaxs='i', type='c', las=1, legend=F, h=1, h.lty=2, h.lwd=1, h.col='black', xlim=c(xmins[n], xmaxs[n]), ylim=c(ymins[n], ymaxs[n]), contour.levels=nlevels, contour.ngrid=ngrid, xcol='Geometric.SizeIterable.None.WholeCell', colors=c('firebrick1','gray31'), alpha=0.2, xlab='', by=c('Mouse'), plot.by=c('Stain','Stim'), cross.lwd=3, cross.fun=mean, cross.args=list(na.rm=T), cross.plot=F, ylab='', log='', logicle.params=list(transX=0.0000001, transY=0.01, base=2), legend.cex=1.0)
			}, pars=list(mar=c(2,2,0,0), mgp=c(3.1,1,0), oma=c(0,0,0,0)), size=c(w*0.5, h*0.5), x='topright', hadj=1, vadj=1)
		}
		j <- j + 1
		n <- n + 1
	}
	i <- i + 1
}
dev.off()

fileToSave <- 'Y:/Jay/R Projects/Dom1/Amount Ratio.png'
scale <- 1.5
# png(file=fileToSave, width=6*scale, height=scale*4, res=600, units='in', family='TT Arial')
par(mfcol=c(2,3))
lims <- getPercentileValues(z$AFrac[is.finite(z$AFrac)], levels=c(0.001, 0.999))
lims <- c(0.01, 0.99)
data.table.plot.all(z[ds %in% DS[1:3] & is.finite(AFrac) & LMB==F], cex.lab=1.5, legend.cex=1.2, density.args=list(adjust=1.2), legend=T, trans.logit=c(T,T), xcol='AFrac', colors=c('coral4','coral1','gray11','gray61'), lwd=2, type='d', by=c('Mouse','LMB'), xlim=lims, xlab='Nuc. Amount / Tot. Amount', plot.by=c('Stain','Stim'))
# dev.off()

fileToSave <- 'Y:/Jay/R Projects/Dom1/z.csv'
# fileToSave <- '/Volumes/Miyamoto/Jay/R Projects/Dom1/x4.csv'
fwrite(z, file=fileToSave)
# x4 <- fread(fileToSave, header=T)


z.sum <- z[LMB==F, list(med=median(NNucAmount)), by=c('Stain','Stim','Mouse','LMB','ds')]
z.sum <- z.sum[, t.test(med[Mouse=='WT'], med[Mouse=='NES'])$p.value, by=c('Stain','Stim')]
plot.func <- function(dt, ...)
{
	par(mgp=c(3,1,0), mar=c(4, 4.7, 2, 1.5))
	boxplot(med ~ Mouse, data = dt, mgp=c(2.9,1,0), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, ...)
	stripchart(med ~ Mouse, vertical = TRUE, data = dt, 
			 method = "jitter", jitter=0, add = TRUE, pch = 21, cex=1.8, bg = 'gray')
	return(1)
}
# daPath <- '/Volumes/Miyamoto/Jay/R Projects/Dom1/PopulationDiffs.pdf'
# pdf(file=daPath, family='Arial', width=6*1.8, height=1.8*2*6/3)
par(mfcol=c(2,3))
z.sum[, list(plot.func(copy(.SD), lwd = 1, ylab = 'Norm. Mann-Whitney U Stat.', las=1, main=paste0(.BY[[1]], ' : ', .BY[[2]]))), by=c('Stain','Stim')]

z[, grp:=factor(paste(Stain, Stim, Mouse, LMB, sep='.'))]

kruskal.test(NNucAmount ~ grp, data=z[Stain=='cRel' & !(Mouse=='NES' & LMB==T)])
kruskal.test(NNucAmount ~ grp, data=z[Stain=='RelA' & !(Mouse=='NES' & LMB==T)])
kruskal.test(NNucAmount ~ grp, data=z[Stain=='IkBa' & !(Mouse=='NES' & LMB==T)])

kruskal.test(Geometric.SizeIterable.None.WholeCell ~ grp, data=z[Stain=='cRel' & !(Mouse=='NES' & LMB==T)])
kruskal.test(Geometric.SizeIterable.None.WholeCell ~ grp, data=z[Stain=='RelA' & !(Mouse=='NES' & LMB==T)])
kruskal.test(Geometric.SizeIterable.None.WholeCell ~ grp, data=z[Stain=='IkBa' & !(Mouse=='NES' & LMB==T)])

z.sum <- z[LMB==F, list(med.WT=median(NNucAmount[Mouse=='WT']), med.NES=median(NNucAmount[Mouse=='NES']), p.value=wilcox.test(NNucAmount[Mouse=='WT'], NNucAmount[Mouse=='NES'])$p.value), by=c('Stain','Stim')]

z.sum <- z[LMB==F, list(med.Stim=median(NNucAmount[Stim==T]), med.Unstim=median(NNucAmount[Stim==F]), p.value=wilcox.test(NNucAmount[Stim==T], NNucAmount[Stim==F])$p.value), by=c('Stain','Mouse')]

z.sum <- z[Mouse=='WT', list(med.LMBplus=median(NNucAmount[LMB==T]), med.LMBminus=median(NNucAmount[LMB==F]), p.value=wilcox.test(NNucAmount[LMB==T], NNucAmount[LMB==F])$p.value), by=c('Stain','Stim')]

z.sum <- z[LMB==F, list(med.WT=median(Geometric.SizeIterable.None.WholeCell[Mouse=='WT']), med.NES=median(Geometric.SizeIterable.None.WholeCell[Mouse=='NES']), p.value=wilcox.test(Geometric.SizeIterable.None.WholeCell[Mouse=='WT'], Geometric.SizeIterable.None.WholeCell[Mouse=='NES'])$p.value), by=c('Stain','Stim')]
z.sum
z.sum <- z[LMB==F, list(med.Stim=median(Geometric.SizeIterable.None.WholeCell[Stim==T]), med.Unstim=median(Geometric.SizeIterable.None.WholeCell[Stim==F]), p.value=wilcox.test(Geometric.SizeIterable.None.WholeCell[Stim==T], Geometric.SizeIterable.None.WholeCell[Stim==F])$p.value), by=c('Stain','Mouse')]
z.sum
z.sum <- z[Mouse=='WT', list(med.LMBplus=median(Geometric.SizeIterable.None.WholeCell[LMB==T]), med.LMBminus=median(Geometric.SizeIterable.None.WholeCell[LMB==F]), p.value=wilcox.test(Geometric.SizeIterable.None.WholeCell[LMB==T], Geometric.SizeIterable.None.WholeCell[LMB==F])$p.value), by=c('Stain','Stim')]
z.sum



# toTest <- x4[(LMB==FALSE) & is.finite(Test) & Valid==T, list(Test=median(Test)), by=c('Stain','Stim','ds','Mouse','LMB2')]
toTest <- z[(LMB==FALSE), getWilcoxStatForEachGroup(.SD, valCol='NRL', by=c('Mouse','LMB2','ds')), by=c('Stain','Stim')]
toTest[, Wmin:=0]
toTest[, Wmax:=n.x*n.y]
toTest <- toTest[order(Stain,Stim,W)]
toTest
toTest[, W2:=(W-Wmin)/(Wmax-Wmin)]
par(mfrow=c(1,1))
plot.func <- function(dt, ...)
{
	par(mgp=c(3,1,0), mar=c(4, 4.7, 2, 1.5))
	boxplot(W2 ~ Mouse, data = dt, mgp=c(2.9,1,0), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, ...)
	stripchart(W2 ~ Mouse, vertical = TRUE, data = dt, 
			 method = "jitter", jitter=0, add = TRUE, pch = 21, cex=1.8, bg = 'gray')
	return(1)
}
daPath <- '/Volumes/Miyamoto/Jay/R Projects/Dom1/PopulationDiffs.pdf'
pdf(file=daPath, family='Arial', width=6*1.8, height=1.8*2*6/3)
par(mfcol=c(2,3))
toTest[, list(plot.func(copy(.SD), lwd = 1, ylab = 'Norm. Mann-Whitney U Stat.', las=1, ylim=range(toTest$W2), main=paste0(.BY[[1]], ' : ', .BY[[2]]))), by=c('Stain','Stim')]
dev.off2(daPath)

wilcox.test.combined(data=x4[(LMB==FALSE) & Stim==T & Stain=='IkBa'], replCols = c('ds'), condCol='Mouse', valCol='NRL', exact=F, two.tailed=F)

