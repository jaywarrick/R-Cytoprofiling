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
diffs <- getDifferences(x3)
standardizeLongData(diffs)
x4 <- rbindlist(list(x3,diffs), use.names = TRUE)
x4$Measurement <- paste(x4$Measurement, x4$MaskChannel, x4$ImageChannel, sep='_')
x4[,MaskChannel:=NULL]
x4[,ImageChannel:=NULL]
x5 <- getWideTable(x4)


# Fix a few things for plotting etc
x5$ImRow <- as.numeric(as.character(x5$ImRow))
x5$ImCol <- as.numeric(as.character(x5$ImCol))
replaceCharacterInColNames(x5, '\\$', '.')
x5$cId <- paste0(x5$Id, ' RC[', x5$ImRow, ',', x5$ImCol, ']')
x5$cId <- paste0(x5$Id, ' RC[', x5$ImRow, ',', x5$ImCol, ']')
names(x5) <- gsub(' ', '', names(x5))
names(x5) <- gsub(':', '_', names(x5))
x5 <- sortColsByName(x5)
shinyData <- x5

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

