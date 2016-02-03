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
     print(temp)
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

     if(length(unique(temp$ImageChannel)) > 1)
     {
          # Calculate differences between channels
          temp2 <- temp[,list(Combo=getComboNames(ImageChannel), Value=getComboDifferences(Value)), by=.(Id,MaskChannel,Measurement)]
          temp2$Measurement <- paste0(temp2$Measurement, '_', temp2$Combo)
          temp2[,Combo:=NULL]
          temp$Measurement <- paste0(temp$Measurement, '_', temp$MaskChannel, '_', temp$ImageChannel)
          temp[,MaskChannel:=NULL]
          temp[,ImageChannel:=NULL]
          temp <- rbindlist(list(temp, temp2))
     }else
     {
          # We only have 1 ImageChannel (or 'none'), so just remove the column and only use the Mask Channel
          temp$Measurement <- paste0(temp$Measurement, '_', temp$MaskChannel)
          temp[,MaskChannel:=NULL]
          temp[,ImageChannel:=NULL]
     }

     temp <- reorganize(temp, c('Id','Label'))
     setcolorder(temp, sort(names(temp)))
     return(temp)
}