library(data.table)

browseShinyData <- function()
{
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataBrowser/ui.R')
     sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataBrowser/server.R')
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
     setcolorder(x, sort(names(x)))
     return(x)
}

standardizeWideData <- function(x)
{
     x[,lapply(.SD, function(x){if(is.numeric(x)){return(scale(x))}else{return(x)}})]
}

standardizeLongData <- function(x)
{
     x[,Value:=scale(Value),by=c('MaskChannel','ImageChannel','Measurement')]
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

