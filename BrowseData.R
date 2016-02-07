rm(list=ls())
source('~/.Rprofile')
# source('D:/GitHub/R-General/.Rprofile')
sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')

library(foreign)
fileTable <- read.arff('/Volumes/JEX Cruncher/JEX Databases/Dominique/temp/JEXData0000000000.arff')

# Preprocess the data
x1a <- fread(fileTable$Value[1])
x1a$Class <- 'MT'
x1b <- fread(fileTable$Value[2])
x1b$Class <- 'WT'
x1 <- rbindlist(list(x1a,x1b), use.names = TRUE)
x1 <- fixColNames(x1)
x1 <- fixMeasurementNames(x1)
setorder(x1, Id, Label, MaskChannel, Measurement, ImageChannel)
x1$Id <- as.character(x1$Id)
x2 <- removeMeasurementNamesContaining(x1, "Phase_Order_2_Rep_0")
x2 <- removeMeasurementNamesContaining(x2, "Phase_Order_4_Rep_0")
x3 <- standardizeLongData(x2)
diffs <- getDifferences(x3)
diffs <- standardizeLongData(diffs)
x4 <- rbindlist(list(x3,diffs), use.names = TRUE)
x4$Measurement <- paste(x4$Measurement, x4$MaskChannel, x4$ImageChannel, sep='_')
x4[,MaskChannel:=NULL]
x4[,ImageChannel:=NULL]
x5 <- getWideTable(x4)
x5 <- fixColNames(x5)


# Fix a few things for plotting etc
x5$ImRow <- as.numeric(as.character(x5$ImRow))
x5$ImCol <- as.numeric(as.character(x5$ImCol))
x5$cId <- paste0(x5$Id, ' RCClass[', x5$ImRow, ',', x5$ImCol, ',', x5$Class, ']')
x5$cId <- paste0(x5$Id, ' RCClass[', x5$ImRow, ',', x5$ImCol, ',', x5$Class, ']')
x5 <- sortColsByName(x5)
shinyData <- x5

write.csv(shinyData, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/test.csv')

# Look at the data
browseShinyData()

library(randomForest)
# Need to remove names like Id, Label, ImRow, ImCol, Z,
dataToTest <- shinyData
dataToTest[, c('cId','Id','Label','ImCol','ImRow','Z'):=NULL]

##### Try some subsets #####
dataToTest <- removeColNamesContaining(dataToTest, "560")
dataToTest$Class <- as.factor(dataToTest$Class)
set.seed(416)
rf <- randomForest(formula= Class ~ ., data=dataToTest, ntree=5, maxnodes=10, importance=TRUE, proximity=TRUE, do.trace=TRUE, keep.forest=TRUE)
rfImp <- data.frame(rf$importance)
rfImp <- rfImp[order(rfImp$MeanDecreaseGini, decreasing=TRUE),]
library(plotly)
plot_ly(rfImp, mode='markers', x=row.names(rfImp), y=rfImp$MeanDecreaseAccuracy, text=row.names(rfImp))

##### Test for index issue ######
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



