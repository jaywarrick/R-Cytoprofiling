rm(list=ls())
source('~/.Rprofile')
sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')

library(foreign)
fileTable <- read.arff('/Volumes/JEX Cruncher/JEX Databases/Dominique/temp/JEXData0000000022.arff')

# Need to set groupCols as "everything but" Value and ImageChannel
# Should be simple, just make vector of names... !(names(duh) %in% c('Value','ImageChannel'))
# temp <- c('a','b')
# duh[,mean(c),by=temp]
shinyData <- getFileTable(fileTable)

shinyData$cId <- paste0(shinyData$Id, ' RC[', shinyData$ImRow, ',', shinyData$ImCol, ']')

browseShinyData()

