library(data.table)
library(foreign)

browseShinyData <- function()
{
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/ui.R')
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/server.R')
     shinyApp(ui=myUI, server=myServer)
}

getComboNames <- function(x, prefix='_')
{
     if(length(x) < 2)
     {
          return(NULL)
     }
     temp <- combn(x, 2)
     #print(temp)
     temp <- paste0(temp[1,],"_minus_",temp[2,])
     temp <- paste0(prefix, temp)
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

getFileTable <- function(fileTable)
{
     temp <- data.table(read.arff(fileTable$Value[1]))
     setorder(temp, Id, Label, MaskChannel, Measurement, ImageChannel)
     return(temp)
}

getColNamesContaining <- function(x, name)
{
     return(names(x)[grepl(name,names(x))])
}

getMeasurementNamesContaining <- function(x, name)
{
     ms <- unique(x$Measurement)
     return(ms[grepl(name,ms)])
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

testFunc2 <- function(x, measurement)
{
     sdx <- sd(x, na.rm=TRUE)
     if(is.na(sdx) || sdx == 0 || is.nan(sdx))
     {
          print(paste0("Removing zero variance measure: ", measurement, '.'))
          return(NULL)
     }else
     {
          return(x)
     }
}
duh2 <- data.table(a=rep(1:3,each=3), b=c(1:3,c(1,1,1),1:3), c=c('a','b','c','d','e','f','g','h','i'))
duh2[,list(Value=testFunc2(b, a)), by=c('a')]

removeNoVarianceMeasurements <- function(x)
{
     # Get a vector with
     temp <- x[,list(stdev=sd(Value)), by=c('MaskChannel','ImageChannel','Measurement')]
     temp <- data.frame(temp[stdev == 0])
     print("Removing measurements with 0 variance...")
     print(temp)
     x[,stdev:=sd(Value), by=c('MaskChannel','ImageChannel','Measurement')]
     x <- x[stdev != 0]
     x[, stdev:=NULL]
     return(x)
}

getNoVarianceCols <- function(x)
{
     tempSD <- function(y){sd(y, na.rm = TRUE)}
     tempNames <- x[,lapply(.SD, tempSD), .SDcols=getNumericCols(x)]
     return(names(tempNames)[as.numeric(as.vector(tempNames))==0])
}

getDifferences <- function(x)
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
     x[,lapply(.SD, function(x){if(is.numeric(x)){return(scale(x))}else{return(x)}})]
}

standardizeLongData <- function(x)
{
     x <- removeNoVarianceMeasurements(x)
     x[,Value:=scale(Value),by=c('MaskChannel','ImageChannel','Measurement')]
     return(x)
}

replaceStringInMeasurementNames <- function(x, old, new)
{
     x[,Measurement:=gsub(old,new,Measurement)]
     return(x)
}

fixMeasurementNames <- function(x)
{
     #      x <- replaceStringInMeasurementNames(x, ' ', '')
     #      x <- replaceStringInMeasurementNames(x, '\\$', '.')
     #      x <- replaceStringInMeasurementNames(x, ':', '_')
     #      return(x)
     x[,Measurement:=gsub(' ','',Measurement)]
     x[,Measurement:=gsub('\\$','.',Measurement)]
     x[,Measurement:=gsub(':','_',Measurement)]
     return(x)
}

fixColNames <- function(x)
{
	replaceStringInColNames(x, ' ', '')
	replaceStringInColNames(x, '\\$', '.')
	replaceStringInColNames(x, ':', '_')
	return(x)
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
     names(x) <- gsub(old, new, names(x))
}

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

