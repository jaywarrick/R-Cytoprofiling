rm(list=ls())
source('~/.Rprofile')
sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')

library(foreign)
fileTable <- read.arff('/Volumes/JEX Cruncher/JEX Databases/Dominique/temp/JEXData0000000000.arff')

# Preprocess the data
x1a <- fread(fileTable$Value[1])
x1a$Class <- 'MT'
x1b <- fread(fileTable$Value[2])
x1b$Class <- 'WT'
x1 <- rbindlist(list(x1a,x1b), use.names = TRUE)
setorder(x1, Id, Label, MaskChannel, Measurement, ImageChannel)
x1$Id <- as.character(x1$Id)
x2 <- removeMeasurementNamesContaining(x1, "Phase_Order_2_Rep_0")
x2 <- removeMeasurementNamesContaining(x2, "Phase_Order_4_Rep_0")
x3 <- standardizeLongData(x2)
diffs1 <- getDifferences(x1)
x3 <- getWideTable(diffs)
standardizeWideData(x3)
x4 <- merge(x2, x3, by=getAllColNamesExcept(x1, 'Value'))
x3$Measurement <- paste(x3$Measurement, x3$MaskChannel, x3$ImageChannel, sep='_')
x3[,MaskChannel:=NULL]
x3[,ImageChannel:=NULL]
x4 <- getWideTable(x3)


# Fix a few things for plotting etc
x4$ImRow <- as.numeric(as.character(x4$ImRow))
x4$ImCol <- as.numeric(as.character(x4$ImCol))
replaceCharacterInColNames(x4, '\\$', '.')
x4$cId <- paste0(x4$Id, ' RC[', x4$ImRow, ',', x4$ImCol, ']')
x4$cId <- paste0(x4$Id, ' RC[', x4$ImRow, ',', x4$ImCol, ']')
x4 <- sortColsByName(x4)
shinyData <- x4

colsToRemove <- x4[,lapply(.SD, function(x){length(which(is.na(x)))>0})]
colsToRemove <- names(colsToRemove)[which(as.logical(as.vector(colsToRemove)))]
x4[,(colsToRemove):=NULL]
x4$Class <- as.factor(x4$Class)

# Get rid of spaces in names and colons as well.
names(shinyData) <- gsub(' ', '', names(shinyData))
names(shinyData) <- gsub(':', '_', names(shinyData))

# Look at the data
browseShinyData()

library(randomForest)
# Need to remove names like Id, Label, ImRow, ImCol, Z,
dataToTest <- shinyData
dataToTest[, c('cId','Id','Label','ImCol','ImRow','Z'):=NULL]
rf <- randomForest(Class ~ ., dataToTest)

library(Hmisc)
# Make the image indicies
shinyData$index1 <- shinyData$ImRow + shinyData$ImCol*max(shinyData$ImRow)
mt <- dataToTest[Class=='MT']
wt <- dataToTest[Class=='WT']

X <- dataToTest[,get(c('ImRow'))]

# Stuff

alpha <- 7/806
ls <- c(1035.29,1012.13, 1014.69)
V<- (pi*(.666/2)^2)*ls*alpha - 2*pi*(0.075/2)^2*ls*alpha

