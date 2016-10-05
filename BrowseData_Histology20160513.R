rm(list=ls())

source('D:/GitHub/R-General/.Rprofile')
source('D:/GitHub/R-Informatics-private/HuMoments.R')
source('D:/GitHub/R-Informatics-private/Zernike.R')
source('D:/GitHub/R-Cytoprofiling/PreProcessingHelpers.R')

# source('~/.Rprofile')
# source('~/Public/DropBox/GitHub/R-Informatics-private/HuMoments.R')
# source('/Users/jaywarrick/Public/DropBox/GitHub/R-Informatics-private/Zernike.R')
# source('~/Public/DropBox/GitHub/R-Cytoprofiling/PreProcessingHelpers.R')

#sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')
library(data.table)
library(foreign) # read.arff
library(gtools) # mixedorder

# A1,2,3 - WT
# A4,5,6 - NES
# A7,9,10 - WT SecOnly
# A8,11,12 - NES SecOnly
#
# B10,11,12 - WT
# B7,8,9 - NES
# B3,4,6 - WT SecOnly
# B1,2,5 - NES SecOnly

# Read in the clump data to just analyze single cells.
clump <- NULL
temp <- read.arff('F:/BoneMarrowSmears/Pt407/File - Clump Data/x0_y1.arff')
temp$Pt <- '407'
clump <- rbind(clump, temp)
temp <- read.arff('F:/BoneMarrowSmears/Pt608/File - Clump Data/x0_y1.arff')
temp$Pt <- '608'
clump <- rbind(clump, temp[,names(temp) %in% names(clump)])
temp <- read.arff('F:/BoneMarrowSmears/Pt609All/File - Clump Data/x0_y1.arff')
temp$Pt <- '609All'
clump <- rbind(clump, temp[,names(temp) %in% names(clump)])
temp <- read.arff('F:/BoneMarrowSmears/Pt610/File - Clump Data/x0_y1.arff')
temp$Pt <- '610'
clump <- rbind(clump, temp[,names(temp) %in% names(clump)])
clump$cId <- paste(clump$Pt, clump$Loc, clump$Id, sep='.')
clump$Id <- as.numeric(as.character(clump$Id))
clump$Loc <- as.numeric(as.character(clump$Loc))
singles <- clump[clump$Value == 1,]
singles$Value <- NULL
singles <- singles[with(singles, order(Pt, Loc, Id)), ]
chosen <- singles[singles$cId %in% resample(singles$cId, 100),]

x <- NULL
temp <- fread('F:/BoneMarrowSmears/Pt407/File - Output CSV Table/x0_y1.csv')
temp$Pt <- '407'
x <- rbind(x, temp)
temp <- fread('F:/BoneMarrowSmears/Pt608/File - Output CSV Table/x0_y1.csv')
temp$Pt <- '608'
x <- rbind(x, temp, use.names=T, fill=F)
temp <- fread('F:/BoneMarrowSmears/Pt609All/File - Output CSV Table/x0_y1.csv')
temp$Pt <- '609All'
temp$C <- NULL
x <- rbind(x, temp, use.names=T, fill=F)
temp <- fread('F:/BoneMarrowSmears/Pt610/File - Output CSV Table/x0_y1.csv')
temp$Pt <- NULL
temp$Pt <- '610'
temp$C <- NULL
x <- rbind(x, temp, use.names=T, fill=F)
x$cId <- paste(x$Pt, x$Loc, x$Id, sep='.')
x$Id <- as.numeric(as.character(x$Id))
x$Loc <- as.numeric(as.character(x$Loc))

x <- x[x$cId %in% chosen$cId]

x1 <- data.table(x)

##### DATA PREPROCESSING #####

fixColNames(x1)

# Make things easier to peruse
setorder(x1, Pt, Loc, Id, Label, Measurement, MaskChannel, ImageChannel)
x1$Id <- as.character(x1$Id) # Avoid standardizing the Id

save(x1, file='F:/BoneMarrowSmears/R Analysis/x1.Rdata')
save(clump, file='F:/BoneMarrowSmears/R Analysis/clump.Rdata')

load(file='F:/BoneMarrowSmears/R Analysis/x1.Rdata')
set.seed(1234)
x1b <- fixLongTableStringsInCol(x1, 'Measurement')
#x1b <- removeMeasurementNamesContaining(x1b,'DNZ')
#x1b <- removeMeasurementNamesContaining(x1b,'NUCw')
x1b <- removeMeasurementNamesContaining(x1b,'ZernikeCircleX')
x1b <- removeMeasurementNamesContaining(x1b,'ZernikeCircleY')
x1b <- removeMeasurementNamesContaining(x1b,'DNZernikeInnerCircleX')
x1b <- removeMeasurementNamesContaining(x1b,'DNZernikeOuterCircleY')
#x1b <- replaceSubStringInAllRowsOfCol(x1b, '395 X 455M', '390 X 440', 'ImageChannel')
# x1b <- replaceSubStringInAllRowsOfCol(x1b, '485 X 525 M', '485 X 525', 'ImageChannel')
# x1b <- replaceSubStringInAllRowsOfCol(x1b, '560 X 607 M', '560 X 607', 'ImageChannel')
# x1b <- replaceSubStringInAllRowsOfCol(x1b, '650 X 705 M', '648 X 684', 'ImageChannel')
# x1b <- x1b[ImageChannel != '560 X 607']
# x1b <- x1b[MaskChannel == 'WholeCell']

# Get rid of cells that have some NA data
# naData <- unique(x1b[!is.finite(Value)]$cId)
# x1b <- x1b[!(cId %in% naData)]
# goodData <- unique(x1b[Measurement == 'ZernikeMag11_NUCwFIXED']$cId)
# x1b <- x1b[(cId %in% goodData)]

# Do some calculations
x2 <- intIntensityNormalizeCentralMoments(x1b)
x2 <- meanNormalizeZernikeMoments(x2)
x2 <- calculateHuMoments(x2)
x2 <- calculateZernikeDotProduct(x2)

# Tempororarily make the table wide to calculate averages of Haralick over the different directions
x2 <- getWideTable(x2)
x2 <- calculateRMSofHaralick(x2)
x2 <- removeExtraneousColumns(x2)

# Get our long table back and reorder
x3 <- getLongTableFromTemplate(x2, x1b)
setorder(x3, Pt, Loc, Id, Label, Measurement, MaskChannel, ImageChannel)

# Perform robust standardization (x-median)/mad (entertain idea of not applying to histogram bins)
x3 <- standardizeLongData(x3, by=c('MaskChannel','ImageChannel','Measurement','Pt'))

# Generate a table of differences between measures for each MaskChannel/ImageChannel/Measurement combination
x3 <- refactor(x3)
diffs <- calculateChannelDifferences(x3)

# Standardize the difference data
diffs <- standardizeLongData(diffs, by=c('MaskChannel','ImageChannel','Measurement','Pt'))

# Merge it with the original dataset, merging MaskChannel and ImageChannel into MeasurementName
x3 <- rbindlist(list(x3,diffs), use.names = TRUE)
x3$Measurement <- paste(x3$Measurement, x3$MaskChannel, x3$ImageChannel, sep='_')
x3[,MaskChannel:=NULL]
x3[,ImageChannel:=NULL]

# Get a wide table for machine learning and plotting
x3 <- getWideTable(x3)

# Fix naming issues introduced by merging MaskChannel and ImageChannel names with Measurement name
x3 <- fixColNames(x3)

# Perform final sorting of columns of data for easier perusing
x3 <- sortColsByName(x3)
# shinyData[,lapply(.SD, function(x){if(is.factor(x)){return(as.factor(x))}else{return(x)}})]
shinyData <- copy(x3)

# Detect an columns with non-finite data and print their names
getNumericCols(shinyData)[as.logical(as.vector(shinyData[,lapply(.SD, function(x){length(which(!is.finite(x)))>0}),.SDcols=getNumericCols(shinyData)]))]

# Write the data for potential analysis outside R
#write.csv(shinyData, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/shinyData_Dom3.csv', row.names=FALSE)
write.csv(shinyData, file='F:/BoneMarrowSmears/R Analysis/shinyData.Rdata', row.names=FALSE)
#shinyData <- read.csv(file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/test.csv')
# Look at the data
#browseShinyData()

###### RANDOM FOREST MACHINE LEARNING #####
library(randomForest)

# Need to remove names like Id, Label, ImRow, ImCol, Z,
dataToTest <- shinyData[Expt==1]
dataToTest[, c('cId','Id','Label','ImCol','ImRow','Z','Loc','file','Expt'):=NULL]
removeColsWithInfiniteVals(dataToTest)
dataToTest$Class <- as.factor(dataToTest$Class)

dataToTrain <- shinyData[Expt==3]
dataToTrain[, c('cId','Id','Label','ImCol','ImRow','Z','Loc','file','Expt'):=NULL]
removeColsWithInfiniteVals(dataToTrain)
dataToTrain$Class <- as.factor(dataToTrain$Class)

#write.csv(dataToTest, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/dataToTest_Dom2.csv', row.names = FALSE)
write.csv(dataToTest, file='G:/Data Tables/dataToTest_Dom1to3_r1234.csv', row.names = FALSE)

dataToTest <- fread('G:/Data Tables/dataToTest_Dom1to3_r1234.csv')

# Set the random seed to reproduce results
set.seed(1234)

# Learn the trees

rf <- randomForest(formula= Class ~ ., data=dataToTrain, ntree=100, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)
rf2 <- randomForest(formula= Class ~ ., data=dataToTrain, ntree=25, maxnodes=5, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)
rf3 <- randomForest(formula= Class ~ ., data=dataToTrain, ntree=25, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)
# Creat interactive plot to browse importance results

length(which(predict(rf2, dataToTest)==dataToTest$Class))/nrow(dataToTest)
library(plotly)
rfImp <- data.frame(rf2$importance)
rfImp$name <- row.names(rfImp)
rfImp <- rfImp[order(rfImp$MeanDecreaseAccuracy, decreasing=TRUE),]
plot(rfImp$MeanDecreaseAccuracy, pch=20, col='deepskyblue3', ylab='Importance', xlab='Feature Rank')
plot_ly(rfImp, mode='markers', x=row.names(rfImp), y=rfImp$MeanDecreaseAccuracy, text=row.names(rfImp))
layout(hovermode="closest")
rfImp$RelMeanDecreaseAccuracy <- rfImp$MeanDecreaseAccuracy/(max(rfImp$MeanDecreaseAccuracy))
rfImp[1:50,c('name','RelMeanDecreaseAccuracy')]
which(grepl('Dot',rfImp$name))
rfImp$name[which(grepl('Dot',rfImp$name))]

i <- 1
#dir.create('/Users/jaywarrick/Documents/MMB/Projects/Dominique/Dom1to3_r1234')
dir.create('G:/Data Tables/Dom1to3_r1234')
for(feature in rfImp$name[1:100])
{
	temp <- gsub(".","_",feature, fixed=T)
	#filePath <- file.path('/Users/jaywarrick/Documents/MMB/Projects/Dominique/Dom2',paste0(i, "_", temp,'.pdf'))
	filePath <- file.path('G:/Data Tables/Dom1to3_r1234',paste0(i, "_", temp,'.pdf'))
	pdf(file=filePath, width=6, height=5)
	plotHist(shinyData,feature)
	dev.off()
	i <- i + 1
}
plotHist(x3,'ZernikeDot40_THISwSEC_WholeCell_390X440_dot_648X684')
plotHist(x3,'ZernikeDot20_THISwSEC_WholeCell_390X440_dot_648X684')
plotHist(x3,'Geometric.SizeIterable_WholeCell_None')
plotHist(x3,'Geometric.SizeIterable_Nuc_None')
plotHist(x3,'Stats.Sum_Nuc_390X440')
plotHist(x3,'Stats.Sum_WholeCell_390X440')
plotHist(x3,'NCR_Sum')
plotHist(x3,'NCR_Mean')

# Plot N/C
x1b[MaskChannel=='Nuc' & ImageChannel=='648 X 684' & Measurement=='Stats.Sum']/x1b[MaskChannel=='Cyt' & ImageChannel=='648 X 684' & Measurement=='Stats.Sum']
x3[,NCR_Sum:=Stats.Sum_Nuc_648X684/Stats.Sum_Cyt_648X684,by='Class']
x3[,NCR_Mean:=Stats.Mean_Nuc_648X684/Stats.Mean_Cyt_648X684,by='Class']

x2[MaskChannel=='Nuc',mean(Geometric.SizeIterable, na.rm=T),by='Class']
x2[MaskChannel=='WholeCell',mean(Geometric.SizeIterable, na.rm=T),by='Class']
x3[,mean(Stats.Sum_Nuc_390X440, na.rm=T),by='Class']


##### DEBUGGING & TESTS #####

names(dataToTest)[as.logical(as.vector(dataToTest[,lapply(.SD, function(x){length(which(is.na(x)))>1})]))]
names(dataToTest)[as.logical(as.vector(dataToTest[,lapply(.SD, function(x){length(which(is.null(x)))>1})]))]

#### Test for index issue #####
order1 <- data.table(ImRow=0:5, ImCol=rep(0:5, each=6), index1=0:35)
order2 <- data.table(ImRow=rep(0:5,each=6), ImCol=0:5, index2=0:35)
order3 <- data.table(ImRow=0:5, ImCol=rep(0:5, each=6), index3=c(0:5,rev(6:11),12:17,rev(18:23),24:29,rev(30:35)))
duh <- data.table(ImRow=c(1,3,5,0),ImCol=c(0,5,4,1), other=c('a','b','c','d'))
shinyData <- merge(order1, shinyData, by=c('ImRow','ImCol'), all=FALSE)
shinyData <- merge(order2, shinyData, by=c('ImRow','ImCol'), all=FALSE)
shinyData <- merge(order3, shinyData, by=c('ImRow','ImCol'), all=FALSE)
library(Hmisc)
mt <- shinyData[Class=='MT']
wt <- shinyData[Class=='WT']
X <- mt[,c('index1','index2','index3'),with=F]
yNames <- getNumericCols(mt)[!(getNumericCols(mt) %in% c('index1','index2','index3'))]
Y <- mt[,yNames, with=F]
temp <- as.data.frame(cor(data.frame(X),data.frame(Y)))
cors <- data.frame(var=names(temp), cor1=as.numeric(as.vector(temp[1,])), cor2=as.numeric(as.vector(temp[2,])), cor3=as.numeric(as.vector(temp[3,])))
cors1 <- cors[order(abs(cors$cor1), decreasing=TRUE),][,c('var','cor1')]
cors2 <- cors[order(abs(cors$cor2), decreasing=TRUE),][,c('var','cor2')]
cors3 <- cors[order(abs(cors$cor3), decreasing=TRUE),][,c('var','cor3')]
nameToGet1 <- as.character(cors1$var[2])
nameToGet2 <- as.character(cors2$var[2])
nameToGet3 <- as.character(cors3$var[2])
plot(mt[,get('index1')], mt[,get(nameToGet1)])
plot(mt[,get('index2')], mt[,get(nameToGet2)])
plot(mt[,get('index3')], mt[,get(nameToGet3)])
# Yay, no trend!

##### Debug stuff #####
x2[Measurement=='HistogramBin_13' & MaskChannel=='WholeCell']

MADs <- shinyData[,lapply(.SD, mad, na.rm=TRUE), .SDcols=getNumericCols(shinyData)]
SDs <- shinyData[,lapply(.SD, sd, na.rm=TRUE), .SDcols=getNumericCols(shinyData)]
SDs$net.imagej.ops.Ops.Stats.StdDev_WholeCell_650X705M

# Probability of flipping 57 heads out of 100 flips given a fair (0.5 fraction) coin is
1-pbinom(4,6,0.30)
qbinom(c(0.025,0.975),100,0.5)

duh <- copy(x1a)
fixColNames(duh)
fixNames(duh, col=c('MaskChannel','ImageChannel','Measurement'))
replaceStringInAllRowsOfCol(duh,'Measurement',"net.imagej.ops.Ops.","")
duh$Measurement <- paste(duh$Measurement, duh$MaskChannel, duh$ImageChannel, sep='_')
duh[,c('MaskChannel','ImageChannel'):=NULL]
duh <- getWideTable(duh)
hist(duh$Stats.Mean_WholeCell_485X525M, breaks=100)
plot_ly(mode='markers',x=duh$Stats.Mean_WholeCell_485X525M-500, y=duh$Stats.Mean_WholeCell_395X455M-500)
layout(p, xaxis = list(type = "log"),
	  yaxis = list(type = "log"), hovermode="closest")
shinyData <- duh
shinyData$Class <- 'MT'
shinyData[800:1600]$Class <- 'WT'
