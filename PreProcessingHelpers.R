library(data.table)

browseShinyData <- function()
{
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/ui.R')
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/server.R')
     shinyApp(ui=myUI, server=myServer)
}

getComboNames <- function(x, prefix='')
{
     if(length(x) < 2)
     {
          return(NULL)
     }
     temp <- combn(x, 2)
     #print(temp)
     temp <- paste0(temp[1,],"_minus_",temp[2,])
     temp <- paste0(prefix, '_', temp)
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

getNoVarianceMeasurements <- function(x)
{
     tempSD <- function(x){sd(x, na.rm=TRUE)}
     temp <- x[,lapply(.SD, tempSD),.SDcols=getNumericCols(x),by=c('MaskChannel','ImageChannel','Measurement')]
     return(unique(temp[Value==0]$Measurement))
}

removeNoVarianceMeasurements <- function(x)
{
     msToRemove <- getNoVarianceMeasurements(x)
     print('Removing the following Measurements with zero variance...')
     for(m in msToRemove)
     {
          print(m)
     }
     return(x[!(Measurement %in% msToRemove)])
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

replaceCharacterInColNames <- function(x, old, new)
{
     names(x) <- gsub(old, new, names(x))
}



