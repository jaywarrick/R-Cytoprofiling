rm(list=ls())
source('D:/GitHub/R-General/.Rprofile')
sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')

x1 <- getXYArffAsTable(dir='D:/Jay/ValidationDataset/Validation/File - Output ARFF Table', file='x19_y0.arff')

# Preprocess the data
x1 <- rbindlist(list(x1a,x1b), use.names = TRUE)
setorder(x1, Id, Label, MaskChannel, Measurement, ImageChannel)
x1$Id <- as.character(x1$Id)
x1$'NUCLEAR LOCALIZATION' <- as.character(x1$'NUCLEAR LOCALIZATION')
x1$'NUCLEAR RADIUS' <- as.character(x1$'NUCLEAR RADIUS')
x1[,c('BLUR','SNR'):=NULL]
x2 <- removeMeasurementNamesContaining(x1, "Phase_Order_2_Rep_0")
x2 <- removeMeasurementNamesContaining(x2, "Phase_Order_4_Rep_0")
x3 <- standardizeLongData(x2)
diffs <- getDifferences(x3)
standardizeLongData(diffs)
x4 <- rbindlist(list(x3,diffs), use.names = TRUE)
x4$Measurement <- paste(x4$Measurement, x4$MaskChannel, x4$ImageChannel, sep='_')
x4[,MaskChannel:=NULL]
x4[,ImageChannel:=NULL]
x4$'NUCLEAR LOCALIZATION' <- as.numeric(x4$'NUCLEAR LOCALIZATION')
x4$'NUCLEAR RADIUS' <- as.numeric(x4$'NUCLEAR RADIUS')
x4 <- fixColNames(x4)
x5 <- getWideTable(x4)


# Fix a few things for plotting etc
x5 <- fixColNames(x5)
x5 <- sortColsByName(x5)
shinyData <- x5

write.csv(shinyData, file='D:/Jay/ValidationDataset/test.csv')

# Look at the data
browseShinyData()

library(randomForest)
# Need to remove names like Id, Label, ImRow, ImCol, Z,
dataToTest <- shinyData
dataToTest[, c('cId','Id','Label','ImCol','ImRow','Z'):=NULL]
rf <- randomForest(Class ~ ., dataToTest)
