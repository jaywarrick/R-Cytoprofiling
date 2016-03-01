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
library(foreign)

# A1,2,3 - WT
# A4,5,6 - NES
# A7,9,10 - WT SecOnly
# A8,11,12 - NES SecOnly
#
# B10,11,12 - WT
# B7,8,9 - NES
# B3,4,6 - WT SecOnly
# B1,2,5 - NES SecOnly

#dir1 <- '/Volumes/Seagate Backup Plus Drive/Data Tables/Dominique 1'
dir1 <- 'G:/Data Tables/Dominique 1'
MTFileList1 <- c('x0_y0.csv')
WTFileList1 <- c('x1_y0.csv')

#dir2 <- '/Volumes/Seagate Backup Plus Drive/Data Tables/Dominique 2'
dir2 <- 'G:/Data Tables/Dominique 2'
MTFileList2 <- c('x0_y0.csv','x0_y1.csv','x0_y2.csv','x0_y3.csv')
WTFileList2 <- c('x0_y4.csv','x0_y5.csv','x0_y6.csv')

#dir3 <- '/Volumes/Seagate Backup Plus Drive/Data Tables/Dominique 3'
dir3 <- 'G:/Data Tables/Dominique 3'
MTFileList3 <- c('x0_y3.csv','x0_y4.csv','x0_y5.csv')#,'x1_y6.csv','x1_y7.csv','x1_y8.csv')
WTFileList3 <- c('x0_y0.csv','x0_y1.csv','x0_y2.csv')#,'x1_y9.csv','x1_y10.csv','x1_y11.csv')
##### DATA PREPROCESSING #####

set.seed(1234)
n <- 3000
tableList <- list()
tableList <- append(tableList, getTableList(dir1, MTFileList1, 'MT', expt=1, sampleSize = n))
tableList <- append(tableList, getTableList(dir1, WTFileList1, 'WT', expt=1, sampleSize = n))
tableList <- append(tableList, getTableList(dir2, MTFileList2, 'MT', expt=2, sampleSize = n))
tableList <- append(tableList, getTableList(dir2, WTFileList2, 'WT', expt=2, sampleSize = n))
tableList <- append(tableList, getTableList(dir3, MTFileList3, 'MT', expt=3, sampleSize = n))
tableList <- append(tableList, getTableList(dir3, WTFileList3, 'WT', expt=3, sampleSize = n))
x1 <- rbindlist(tableList, use.names=TRUE)
# Check number sampled
#length(unique(x1$cId[grepl('2.x',x1$cId,fixed=T) & x1$Class=='MT']))
fixColNames(x1)

# Make things easier to peruse
setorder(x1, Expt, file, Id, Label, MaskChannel, Measurement, ImageChannel)
x1$Id <- as.character(x1$Id) # Avoid standardizing the Id

save(x1, file='/Volumes/Seagate Backup Plus Drive/Data Tables/x1_sampled.Rdata')
save(x1, file='G:/Data Tables/x1_sampled.Rdata')

load(file='G:/Data Tables/x1_sampled.Rdata')
set.seed(1234)
x1b <- fixLongTableStringsInCol(x1, 'Measurement')
#x1b <- removeMeasurementNamesContaining(x1b,'DNZ')
#x1b <- removeMeasurementNamesContaining(x1b,'NUCw')
x1b <- removeMeasurementNamesContaining(x1b,'ZernikeCircleX')
x1b <- removeMeasurementNamesContaining(x1b,'ZernikeCircleY')
x1b <- removeMeasurementNamesContaining(x1b,'DNZernikeInnerCircleX')
x1b <- removeMeasurementNamesContaining(x1b,'DNZernikeOuterCircleY')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '395 X 455M', '390 X 440', 'ImageChannel')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '485 X 525 M', '485 X 525', 'ImageChannel')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '560 X 607 M', '560 X 607', 'ImageChannel')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '650 X 705 M', '648 X 684', 'ImageChannel')
x1b <- x1b[ImageChannel != '560 X 607']
#x1b <- x1b[MaskChannel == 'WholeCell']

# Get rid of cells that have some NA data
naData <- unique(x1b[!is.finite(Value)]$cId)
x1b <- x1b[!(cId %in% naData)]
goodData <- unique(x1b[Measurement == 'ZernikeMag11_NUCwFIXED']$cId)
x1b <- x1b[(cId %in% goodData)]

# Do some calculations
x2 <- intIntensityNormalizeCentralMoments(x1b)
x2 <- meanNormalizeZernikeMoments(x2)
x2 <- calculateHuMoments(x2)
x2 <- calculateZernikeDotProduct(x2)

# Tempororarily make the table wide to calculate averages of Haralick over the different directions
x2b <- getWideTable(x2)
x2b <- calculateRMSofHaralick(x2b) 
x2b <- removeExtraneousColumns(x2b)

# Get our long table back and reorder
x3 <- getLongTableFromTemplate(x2b, x2)
setorder(x3, Expt, file, Id, Label, MaskChannel, Measurement, ImageChannel)

# Perform robust standardization (x-median)/mad (entertain idea of not applying to histogram bins)
x3 <- standardizeLongData(x3)

# Generate a table of differences between measures for each MaskChannel/ImageChannel/Measurement combination
x3 <- refactor(x3)
diffs <- calculateChannelDifferences(x3)

# Standardize the difference data
diffs <- standardizeLongData(diffs)

# Merge it with the original dataset, merging MaskChannel and ImageChannel into MeasurementName
x4 <- rbindlist(list(x3,diffs), use.names = TRUE)
x4$Measurement <- paste(x4$Measurement, x4$MaskChannel, x4$ImageChannel, sep='_')
x4[,MaskChannel:=NULL]
x4[,ImageChannel:=NULL]

# Get a wide table for machine learning and plotting
x5 <- getWideTable(x4)

# Fix naming issues introduced by merging MaskChannel and ImageChannel names with Measurement name
x5 <- fixColNames(x5)

# Perform final sorting of columns of data for easier perusing
x5 <- sortColsByName(x5)
# shinyData[,lapply(.SD, function(x){if(is.factor(x)){return(as.factor(x))}else{return(x)}})]
shinyData <- copy(x5)

getNumericCols(shinyData)[as.logical(as.vector(shinyData[,lapply(.SD, function(x){length(which(!is.finite(x)))>0}),.SDcols=getNumericCols(shinyData)]))]

# Write the data for potential analysis outside R
#write.csv(shinyData, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/shinyData_Dom3.csv', row.names=FALSE)
write.csv(shinyData, file='G:/Data Tables/shinyData_Dom1to3_r1234.csv', row.names=FALSE)
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
plotHist(x5,'ZernikeDot40_THISwSEC_WholeCell_390X440_dot_648X684')
plotHist(x5,'ZernikeDot20_THISwSEC_WholeCell_390X440_dot_648X684')
plotHist(x5,'Geometric.SizeIterable_')


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
