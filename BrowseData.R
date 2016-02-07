rm(list=ls())
source('~/.Rprofile')
# source('D:/GitHub/R-General/.Rprofile')
sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')

library(foreign)
fileTable <- read.arff('/Volumes/JEX Cruncher/JEX Databases/Dominique/temp/JEXData0000000000.arff')

##### DATA PREPROCESSING #####

# Read in the data into a single table
x1a <- fread(fileTable$Value[1])
x1a$Class <- 'MT'
x1b <- fread(fileTable$Value[2])
x1b$Class <- 'WT'
x1 <- rbindlist(list(x1a,x1b), use.names = TRUE)

# Fix some naming issues from ImageJ
x1 <- fixColNames(x1)
x1 <- fixMeasurementNames(x1)

# Make things easier to peruse
setorder(x1, Id, Label, MaskChannel, Measurement, ImageChannel)
x1$Id <- as.character(x1$Id) # Void standardizing the Id

# Remove measures that vascilate between 180 deg and 0 deg (i.e., no diff but numberwise they have a non-zero difference)
x2 <- removeMeasurementNamesContaining(x1, "Phase_Order_2_Rep_0")
x2 <- removeMeasurementNamesContaining(x2, "Phase_Order_4_Rep_0")

# Tempororarily make the table wide to calculate averages of Haralick over the different directions
x2b <- getWideTable(x2)
x2b <- calculateRMSofHaralick(x2b)

# Get our long table back and reorder
x3 <- getLongTableFromTemplate(x2b, x2)
setorder(x3, Id, Label, MaskChannel, Measurement, ImageChannel)

# Perform robust standardization (x-median)/mad (entertain idea of not applying to histogram bins)
x3 <- standardizeLongData(x3)

# Generate a table of differences between measures for each MaskChannel/ImageChannel/Measurement combination
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

# Convert ImRow and ImCol to numeric for looking at locational correlations
x5$ImRow <- as.numeric(as.character(x5$ImRow))
x5$ImCol <- as.numeric(as.character(x5$ImCol))

# Generate a more informative text string for pt labels in plotly
x5$cId <- paste0(x5$Id, ' RCClass[', x5$ImRow, ',', x5$ImCol, ',', x5$Class, ']')
x5$cId <- paste0(x5$Id, ' RCClass[', x5$ImRow, ',', x5$ImCol, ',', x5$Class, ']')

# Perform final sorting of columns of data for easier perusing
x5 <- sortColsByName(x5)
shinyData <- x5

# Write the data for potential analysis outside R
write.csv(shinyData, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/test.csv')

# Look at the data
browseShinyData()

###### RANDOM FOREST MACHINE LEARNING #####
library(randomForest)

# Need to remove names like Id, Label, ImRow, ImCol, Z,
dataToTest <- shinyData
dataToTest[, c('cId','Id','Label','ImCol','ImRow','Z'):=NULL]

# The 560 channel is potentiall suspect due to image acquistion issues. Remove to avoid potential bias.
dataToTest <- removeColNamesContaining(dataToTest, "560")
dataToTest$Class <- as.factor(dataToTest$Class)

# Set the random seed to reproduce results
set.seed(416)

# Learn the trees
rf <- randomForest(formula= Class ~ ., data=dataToTest, ntree=50, maxnodes=3, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)

# Get an ordered table of importances (accuracy and gini)
rfImp <- data.frame(rf$importance)
rfImp <- rfImp[order(rfImp$MeanDecreaseAccuracy, decreasing=TRUE),]

# Creat interactive plot to browse importance results
library(plotly)
plot_ly(rfImp, mode='markers', x=row.names(rfImp), y=rfImp$MeanDecreaseAccuracy, text=row.names(rfImp))

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
