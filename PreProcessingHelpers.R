library(data.table)
library(foreign)

##### Visualization #####

browseShinyData <- function()
{
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/ui.R')
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/server.R')
     shinyApp(ui=myUI, server=myServer)
}

##### Table IO #####

getXYCSVsAsTableFromDir <- function(dir, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
{
     ret <- list()
     fList <- list.files(path = dir, recursive = TRUE)
     for(f in fList)
     {
          if((grepl('x', f) || grepl('y', f)) & grepl('.csv', f))
          {
               fileName <- strsplit(f, "\\.")[[1]][1]
               ret[[fileName]] <- getXYCSVAsTable(dir, f, xName, xExpression, yName, yExpression)
          }
     }
     retTable <- rbindlist(ret)
     return(retTable)
}

getXYCSVAsTable <- function(dir, file, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
{
     fileName <- strsplit(file, "\\.")[[1]][1]
     xy <- strsplit(fileName, "_")[[1]]
     y <- as.numeric(substr(xy[1],2,nchar(xy[1])))
     x <- as.numeric(substr(xy[2],2,nchar(xy[2])))
     xVal <- eval(parse(text=xExpression))
     yVal <- eval(parse(text=yExpression))
     print(paste0('Reading ', file.path(dir,file), ' as ', xName, '=', xVal, ', ', yName, '=', yVal, '.'))
     theTable <- fread(file.path(dir,file))
     theTable[,(xName),with=FALSE] <- xVal
     theTable[,(yName),with=FALSE] <- yVal
     return(theTable)
}

getXYArffsAsTableFromDir <- function(dir, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
{
     ret <- list()
     fList <- list.files(path = dir, recursive = TRUE)
     for(f in fList)
     {
          if((grepl('x', f) || grepl('y', f)) & grepl('.arff', f))
          {
               fileName <- strsplit(f, "\\.")[[1]][1]
               ret[[fileName]] <- getXYArffAsTable(dir, f, xName, xExpression, yName, yExpression)
          }
     }
     retTable <- rbindlist(ret)
     return(retTable)
}

getXYArffAsTable <- function(dir, file, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
{
     fileName <- strsplit(file, "\\.")[[1]][1]
     xy <- strsplit(fileName, "_")[[1]]
     y <- as.numeric(substr(xy[1],2,nchar(xy[1])))
     x <- as.numeric(substr(xy[2],2,nchar(xy[2])))
     xVal <- eval(parse(text=xExpression))
     yVal <- eval(parse(text=yExpression))
     print(paste0('Reading ', file.path(dir,file), ' as ', xName, '=', xVal, ', ', yName, '=', yVal, '.'))
     theTable <- read.arff(file.path(dir,file))
     theTable[,xName] <- xVal
     theTable[,yName] <- yVal
     return(data.table(theTable))
}

##### Wide Table Operations #####

removeColsWithInfiniteVals <- function(x)
{
     duh <- x[,lapply(.SD, function(y){length(which(!is.finite(y))) > 0}), .SDcols=getNumericCols(x)]
     duh2 <- getNumericCols(x)[as.logical(as.vector(duh))]
     if(length(duh2 > 0))
     {
          print("Removing cols with infinite values...")
     }
     for(col in duh2)
     {
          print(col)
          x[,(col):=NULL]
     }
}

getColNamesContaining <- function(x, name)
{
     return(names(x)[grepl(name,names(x))])
}

removeColNamesContaining <- function(x, name)
{
     colsToRemove <- getColNamesContaining(x,name)
     print(paste0("Removing colums with names containing '", name, "'"))
     for(colToRemove in colsToRemove)
     {
          print(colToRemove)
          x[,(colToRemove):=NULL]
     }
     return(x)
}

fixColNames <- function(x)
{
     replaceStringInColNames(x, ' ', '')
     replaceStringInColNames(x, '\\$', '.')
     replaceStringInColNames(x, ':', '_')
}

getAllColNamesExcept <- function(x, names)
{
     return(names(x)[!(names(x) %in% names)])
}

getNumericCols <- function(x)
{
     return(names(x)[unlist(x[,lapply(.SD, is.numeric)])])
}

getNonNumericCols <- function(x)
{
     return(names(x)[!unlist(x[,lapply(.SD, is.numeric)])])
}

replaceStringInColNames <- function(x, old, new)
{
     oldNames <- names(x)
     newNames <- gsub(old, new, names(x))
     setnames(x, oldNames, newNames)
}

getWideTable <- function(x)
{
     idCols <- getAllColNamesExcept(x, c('Value','Measurement'))
     x <- reorganize(x, idCols)
     x <- sortColsByName(x);
     return(x)
}

sortColsByName <- function(x)
{
     setcolorder(x, sort(names(x)))
}

standardizeWideData <- function(x)
{
     removeNoVarianceCols(x)
     robustScale <- function(x)
     {
          m <- median(x, na.rm=TRUE)
          return((x-m)/mad(x, center=m, na.rm=TRUE))
     }
     x[,lapply(.SD, function(x){if(is.numeric(x)){return(robustScale(x))}else{return(x)}})]
}

removeNoVarianceCols <- function(x)
{
     namesToRemove <- getNoVarianceCols(x)
     if(length(namesToRemove) > 0)
     {
          print("Removing cols with a variance of zero...")
          for(name in namesToRemove)
          {
               print(name)
               x[,(name):=NULL]
          }
     }
}

getNoVarianceCols <- function(x)
{
     tempSD <- function(y){sd(y, na.rm = TRUE)}
     tempNames <- x[,lapply(.SD, tempSD), .SDcols=getNumericCols(x)]
     return(names(tempNames)[as.numeric(as.vector(tempNames))==0])
}

##### Long Table Operations #####

getLongTable <- function(x, idCols, measurementName='Measurement', valueName='Value')
{
     return(melt(x, getAllColNamesExcept(x, idCols), variable.name=measurementName, value.name=valueName, na.rm=TRUE))
}

getLongTableFromTemplate <- function(x, longTemplate)
{
     return(getLongTable(x, idCols=getAllColNamesExcept(x, getAllColNamesExcept(longTemplate, c('Measurement','Value')))))
}

getMeasurementNamesContaining <- function(x, name)
{
     ms <- unique(x$Measurement)
     return(ms[grepl(name,ms)])
}

removeMeasurementNamesContaining <- function(x, name)
{
     namesToRemove <- getMeasurementNamesContaining(x, name)
     print("Removing the following Measurements...")
     for(name in namesToRemove)
     {
          print(name)
     }
     x <- x[!(Measurement %in% namesToRemove)]
     return(x)
}

standardizeLongData <- function(x, by=c('MaskChannel','ImageChannel','Measurement'))
{
     robustScale <- function(x)
     {
          m <- median(x, na.rm=TRUE)
          return((x-m)/mad(x, center=m, na.rm=TRUE))
     }
     x <- removeNoMADMeasurements(x)
     x[,Value:=robustScale(Value),by=by]
     return(x)
}

removeNoVarianceMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement'))
{
     # See if we have any columns to remove and record the info for reporting
     temp <- x[,list(stdev=sd(get(val))), by=by]
     temp <- data.frame(temp[stdev == 0])
     print("Removing measurements with 0 variance...")
     print(temp)
     # Tempororarily add a column in the table with stdev in it
     x[,stdev:=sd(get(val)), by=by]
     y <- x[stdev != 0]
     x[, stdev:=NULL]
     y[, stdev:=NULL]
     return(y)
}

removeNoMADMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement'))
{
     # See if we have any columns to remove and record the info for reporting
     temp <- x[,list(MAD=mad(get(val), na.rm=TRUE)), by=by]
     temp <- data.frame(temp[MAD == 0])
     print("Removing measurements with 0 MAD...")
     print(temp)
     # Tempororarily add a column in the table with stdev in it
     x[,MAD:=mad(get(val), na.rm=TRUE), by=by]
     y <- x[MAD != 0]
     x[, MAD:=NULL]
     y[, MAD:=NULL]
     return(y)
}

removeNoVarianceMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement'))
{
     # See if we have any columns to remove and record the info for reporting
     temp <- x[,list(stdev=sd(get(val))), by=by]
     temp <- data.frame(temp[stdev == 0])
     print("Removing measurements with 0 variance...")
     print(temp)
     # Tempororarily add a column in the table with stdev in it
     x[,stdev:=sd(get(val)), by=by]
     y <- x[stdev != 0]
     x[, stdev:=NULL]
     y[, stdev:=NULL]
     return(y)
}

replaceSubStringInAllRowsOfCol <- function(x, old, new, col)
{
     x[,c(col):=gsub(old,new,get(col))]
}

fixNames <- function(x, col)
{
     for(colName in col)
     {
          replaceSubStringInAllRowsOfCol(x,' ','',colName)
          replaceSubStringInAllRowsOfCol(x,'\\$','.',colName)
          replaceSubStringInAllRowsOfCol(x,':','_',colName)
     }
}

##### Feature Calculations #####

calculateChannelDifferences <- function(x)
{
     if(length(unique(x$ImageChannel)) > 1)
     {
          # Calculate differences between channels
          idCols <- getAllColNamesExcept(x, c('Value','ImageChannel'))
          x2 <- x[,list(ImageChannel=getComboNames(ImageChannel), Value=getComboDifferences(Value)), by=idCols]
     }else
     {
          # return an empty table with the same columns as provided
          return(x[FALSE])
     }
}

getComboNames <- function(x)
{
     if(length(x) < 2)
     {
          return(NULL)
     }
     temp <- combn(x, 2)
     #print(temp)
     temp <- paste0(temp[1,],"_minus_",temp[2,])
     return(temp)
}

getComboDifferences <- function(x)
{
     if(length(x) < 2)
     {
          return(NULL)
     }
     temp <- combn(x, 2)
     temp <- temp[1,]-temp[2,]
     return(temp)
}

calculateRMSofHaralick <- function(x, removeOriginalHaralickMeasures=FALSE)
{
     # If keeping Haralick features, combine measures for each direction by averaging to make "rotationally invariant".
     # Find all names with Horizontal in them
     hNames <- getColNamesContaining(x, 'Horizontal')
     vNames <- gsub("Horizontal", "Vertical", hNames)
     dNames <- gsub("Horizontal", "Diagonal", hNames)
     adNames <- gsub("Horizontal", "AntiDiagonal", hNames)
     avgNames <- gsub("Horizontal", "Avg", hNames)

     haralickNames <- data.frame(H=hNames, V=vNames, D=dNames, AD=adNames, avg=avgNames, stringsAsFactors=FALSE)
     myfunc <- function(row, theNames)
     {
          return(mean(row[,theNames$H] + row[,theNames$V] + row[,theNames$D] + row[,theNames$AD]))
     }

     x <- data.frame(x)
     for(i in 1:nrow(haralickNames))
     {
          x[,haralickNames[i,5]] <- (x[,haralickNames[i,1]] + x[,haralickNames[i,2]] + x[,haralickNames[i,3]] + x[,haralickNames[i,4]])/4
          if(removeOriginalHaralickMeasures)
          {
               x <- x[,!(names(x) %in% as.character(haralickNames[i,1:4]))]
          }
     }

     return(data.table(x))
}

getColors <- function(pointClasses)
{
     ret <- rep('rgb(0,0,1,0.2)', length(pointClasses))
     ret[pointClasses == 'MT'] <- 'rgb(1,0,0,0.2)'
     return(ret)
}

##### Testing #####

# testFunc2 <- function(x, measurement)
# {
#      sdx <- sd(x, na.rm=TRUE)
#      if(is.na(sdx) || sdx == 0 || is.nan(sdx))
#      {
#           print(paste0("Removing zero variance measure: ", measurement, '.'))
#           return(NULL)
#      }else
#      {
#           return(x)
#      }
# }
# duh2 <- data.table(a=rep(1:3,each=3), b=c(1:3,c(1,1,1),1:3), c=c('a','b','c','d','e','f','g','h','i'))
# duh2[,list(Value=testFunc2(b, a)), by=c('a')]

