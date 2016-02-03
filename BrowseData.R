rm(list=ls())
source('~/.Rprofile')
sourceGitHubFile(user)

library(foreign)
fileTable <- read.arff('/Users/jaywarrick/Documents/JEX/Feature Extraction/temp/JEXData0000000053.arff')

x <- getFileTable(fileTable)

browseData(x)
