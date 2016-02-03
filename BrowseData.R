rm(list=ls())
source('~/.Rprofile')
sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')

library(foreign)
fileTable <- read.arff('/Users/jaywarrick/Documents/JEX/Feature Extraction/temp/JEXData0000000053.arff')

shinyData <- getFileTable(fileTable)

browseShinyData()
