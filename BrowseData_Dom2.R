rm(list=ls())
source('~/.Rprofile')
source('~/Public/DropBox/GitHub/R-Informatics-private/HuMoments.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-Informatics-private/Zernike.R')
# source('D:/GitHub/R-General/.Rprofile')
source('~/Public/DropBox/GitHub/R-Cytoprofiling/PreProcessingHelpers.R')
#sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')
library(data.table)
library(foreign)

fileList <- c('x0_y0.csv','x0_y1.csv','x0_y2.csv','x0_y3.csv','x0_y4.csv','x0_y5.csv','x0_y6.csv')
# dataMT <- getData(db='/Volumes/JEX Cruncher/JEX Databases/Dominique', ds='Mutant vs WT', x=0, y=0, type='File', name='Output CSV Table')
# dataWT <- getData(db='/Volumes/JEX Cruncher/JEX Databases/Dominique', ds='Mutant vs WT', x=1, y=0, type='File', name='Output CSV Table')
# dataFE <- getData(db='/Users/jaywarrick/Documents/JEX/Feature Extraction', ds='Dataset Name', x=0, y=0, type='File', name='Output CSV Table')
##### DATA PREPROCESSING #####

# Read in the data into a single table
# x1a <- fread('/Users/jaywarrick/Desktop/A Sandbox/JEXData0000000003.csv')

tableList <- list()
for(i in 1:7)
{
     temp <- fread(file.path('/Users/jaywarrick/Documents/MMB/Projects/Dominique/x0_y6', fileList[i]))
     if(i <= 4)
     {
          temp$Class <- 'WT'
     }
     else
     {
          temp$Class <- 'MT'
     }
     fixColNames(temp)
     fixNames(temp, c('Measurement','ImageChannel','MaskChannel'))
     tableList[[i]] <- copy(temp)
}
x1 <- rbindlist(tableList, use.names = TRUE)
replaceSubStringInAllRowsOfCol(x1,'net.imagej.ops.Ops.','','Measurement')
replaceSubStringInAllRowsOfCol(x1,'_Order_','','Measurement')
replaceSubStringInAllRowsOfCol(x1,'_Rep_','','Measurement')

# Make things easier to peruse
setorder(x1, Id, Label, MaskChannel, Measurement, ImageChannel)
x1$Id <- as.character(x1$Id) # Avoid standardizing the Id

# Do some calculations
x2 <- copy(x1)
x2 <- intIntensityNormalizeCentralMoments(x2)
x2 <- meanNormalizeZernikeMoments(x2)
x2 <- calculateHuMoments(x2)
x2 <- calculateZernikeDotProduct(x2)

# Tempororarily make the table wide to calculate averages of Haralick over the different directions
x2b <- getWideTable(x2)
x2b <- calculateRMSofHaralick(x2b)
x2b <- removeExtraneousColumns(x2b)

# Get our long table back and reorder
x3 <- getLongTableFromTemplate(x2b, x2)
setorder(x3, Id, Label, MaskChannel, Measurement, ImageChannel)

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

# Write the data for potential analysis outside R
write.csv(shinyData, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/shinyData_Dom2.csv', row.names=FALSE)
#shinyData <- read.csv(file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/test.csv')
# Look at the data
#browseShinyData()

###### RANDOM FOREST MACHINE LEARNING #####
library(randomForest)

# Need to remove names like Id, Label, ImRow, ImCol, Z,
dataToTest <- shinyData
dataToTest[, c('cId','Id','Label','ImCol','ImRow','Z'):=NULL]
removeColsWithInfiniteVals(dataToTest)

# The 560 channel is potentially suspect due to image acquistion issues. Remove to avoid potential bias.
dataToTest <- removeColNamesContaining(dataToTest, "560")
dataToTest$Class <- as.factor(dataToTest$Class)

write.csv(dataToTest, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/dataToTest_Dom2.csv', row.names = FALSE)

# Set the random seed to reproduce results
set.seed(416)

# Learn the trees

rf <- randomForest(formula= Class ~ ., data=dataToTest, ntree=100, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)
rf2 <- randomForest(formula= Class ~ ., data=dataToTest, ntree=25, maxnodes=10, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)
rf3 <- randomForest(formula= Class ~ ., data=dataToTest, ntree=25, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)
# Creat interactive plot to browse importance results
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
dir.create('/Users/jaywarrick/Documents/MMB/Projects/Dominique/Dom2')
for(feature in rfImp$name[1:100])
{
     temp <- gsub(".","_",feature, fixed=T)
     filePath <- file.path('/Users/jaywarrick/Documents/MMB/Projects/Dominique/Dom2',paste0(i, "_", temp,'.pdf'))
     pdf(file=filePath, width=6, height=5)
     plotHist(shinyData,feature)
     dev.off()
     i <- i + 1
}


##### DEBUGGING & TESTS #####

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
