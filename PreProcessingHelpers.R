library(data.table)
library(foreign)
library(spatstat)
library(shotGroups)

##### GitHub Tools #####

sourceGitHubFile <- function(user, repo, branch, file)
{
     require(curl)
     destfile <- tempfile()
     fileToGet <- paste0("https://raw.githubusercontent.com/", user, "/", repo, "/", branch, "/", file)
     curl_download(url=fileToGet, destfile)
     source(destfile)
}

##### Visualization #####

browseShinyData <- function()
{
	sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/ui.R')
	sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/server.R')
	shinyApp(ui=myUI, server=myServer)
}

##### General #####

resample <- function(x, ...)
{
     x[sample.int(length(x), ...)]
}

getLocsFromRCs <- function(r, c, numRows)
{
	r + max(numRows) * c
}

sind <- function(x)
{
	return(sin(x*pi/180))
}

cosd <- function(x)
{
	return(cos(x*pi/180))
}

tand <- function(x)
{
	return(tan(x*pi/180))
}

refactor <- function(x)
{
	return(x[,lapply(.SD, function(x){if(is.factor(x)){factor(x)}else{x}})])
}

##### Table IO #####

getTableListFromDB <- function(db, ds, x, y, objectName, isArff=F, storeFilePath=F, class=NULL, assignClass=T, expt=NULL, repl=NULL, sampleSize=NULL, colsToRemove = c(), cIdCols = c(), fsep='\\\\')
{
     tableList <- list()
     
     jData <- getData(db=db, ds=ds, x=x, y=y, type='File', name=objectName)
     fileList <- jData$fileList
     
     if(!is.null(sampleSize))
     {
          subSampleSize <- sampleSize / length(fileList)
     }
     
     # For each file in the fileList
     for(f in fileList)
     {
          # Read the file in
          print(paste0('Reading file: ', f))
          
          if(isArff)
          {
               library(foreign)
               temp <- data.table(read.arff(f))
          }
          else
          {
               temp <- fread(f)
          }
          
          # Store the filepath that was imported if desired
          if(storeFilePath)
          {
               temp$File <- f	
          }
          
          # Store the name/number of the experiment/replicate associated with this file
          if(!is.null(expt))
          {
               temp$Expt <- expt
          }
          if(!is.null(replicate))
          {
               temp$Repl <- repl
          }
          
          # Create/Assign a 'Class' column
          if(!is.null(class) && assignClass)
          {
               temp$Class <- class
          }
          else if(!is.null(class) && !assignClass)
          {
               setnames(temp,class,'Class')
               temp$Class <- as.character(temp$Class)
          }
          
          # Create a column with a complex Id that will be completely unique for each sample
          idColsFound <- cIdCols[cIdCols %in% names(temp)]
          cat(names(temp))
          if(length(idColsFound) == 0)
          {
               cat(names(temp))
               stop('Must specify cIds for this function to enable sampling.')
          }
          if(length(idColsFound) != length(cIdCols))
          {
               warning(cat('The specified cIdCols (', cIdCols[!(cIdCols %in% names(temp))], 'is/are not column names of the table being retrieved... (', names(temp), ')'))
          }
          temp[,c('cId'):=paste(mapply(function(x){unique(as.character(x))}, mget(idColsFound)), collapse='.'), by=idColsFound]
          
          # put the complex Id first and the class column last
          setcolorder(temp, c('cId', names(temp)[names(temp) != 'cId']))
          
          # Put the 'Class' column as the last column of the table
          if('Class' %in% names(temp))
          {
               setcolorder(temp, c(names(temp)[names(temp) != 'Class'], 'Class'))
          }
          
          # Remove specified columns from the data
          for(tempCol in colsToRemove)
          {
               if(tempCol %in% names(temp))
               {
                    temp[,c(tempCol) := NULL]
               }
               else
               {
                    warning(paste(tempCol, 'is not a column of the data table so it cannot be removed'))
               }
          }
          
          # Grab the randomly sampled rows of the file
          if(!is.null(sampleSize))
          {
               rIds <- trySample(unique(temp$cId), subSampleSize)
               temp <- temp[cId %in% rIds]
          }
          
          # Print the column names for a little feedback
          print(names(temp))
          
          # Append this table to the list of tables provided.
          tableList <- append(tableList, list(temp))
     }
     return(tableList)
}

getTableList <- function(dir, fileList, class, expt, sampleSize=NULL, cellIds=NULL)
{
	if(!is.null(sampleSize))
	{
		subSampleSize <- sampleSize / length(fileList)
	}
	tableList <- list()
	for(f in fileList)
	{
		print(paste0('Reading file: ', file.path(dir, f)))
		temp <- fread(file.path(dir, f))
		temp$Class <- class
		temp$Expt <- expt
		if('Z' %in% names(temp))
		{
			temp[,Z:=NULL]
		}
		if('A' %in% names(temp))
		{
			temp[,A:=NULL]
		}
		if('B' %in% names(temp))
		{
			temp[,B:=NULL]
		}
		if(!('ImRow' %in% names(temp)))
		{
			temp$ImRow <- 1
		}
		if(!('ImCol' %in% names(temp)))
		{
			temp$ImCol <- 1
		}
		if(!('Loc' %in% names(temp)))
		{
			temp[,Loc:=getLocsFromRCs(ImRow, ImCol, max(ImRow) + 1)]
		}
		temp$file <- f
		temp$cId <- paste(temp$Expt, temp$file, temp$Loc, temp$Id, sep='.')
		temp[,ImRow:=NULL]
		temp[,ImCol:=NULL]
		if(!is.null(cellIds))
		{
		     rIds <- cellIds
		     temp <- temp[cId %in% rIds]
		}
		if(!is.null(sampleSize))
		{
			rIds <- trySample(unique(temp$cId), subSampleSize)
			temp <- temp[cId %in% rIds]
		}
		tableList <- append(tableList, list(temp))
	}
	return(tableList)
}

getTableList <- function(dir, featureTable, colocTable, class, expt, sampleSize=NULL, cellIds=NULL)
{
     if(!is.null(sampleSize))
     {
          subSampleSize <- sampleSize / length(fileList)
     }
     tableList <- list()
     for(f in fileList)
     {
          print(paste0('Reading file: ', file.path(dir, f)))
          temp <- fread(file.path(dir, f))
          temp$Class <- class
          temp$Expt <- expt
          if('Z' %in% names(temp))
          {
               temp[,Z:=NULL]
          }
          if('A' %in% names(temp))
          {
               temp[,A:=NULL]
          }
          if('B' %in% names(temp))
          {
               temp[,B:=NULL]
          }
          if(!('ImRow' %in% names(temp)))
          {
               temp$ImRow <- 1
          }
          if(!('ImCol' %in% names(temp)))
          {
               temp$ImCol <- 1
          }
          if(!('Loc' %in% names(temp)))
          {
               temp[,Loc:=getLocsFromRCs(ImRow, ImCol, max(ImRow) + 1)]
          }
          temp$file <- f
          temp$cId <- paste(temp$Expt, temp$file, temp$Loc, temp$Id, sep='.')
          temp[,ImRow:=NULL]
          temp[,ImCol:=NULL]
          if(!is.null(cellIds))
          {
               rIds <- cellIds
               temp <- temp[cId %in% rIds]
          }
          if(!is.null(sampleSize))
          {
               rIds <- trySample(unique(temp$cId), subSampleSize)
               temp <- temp[cId %in% rIds]
          }
          tableList <- append(tableList, list(temp))
     }
     return(tableList)
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

##### Wide Table Operations #####

divideColAByColB <- function(x, colA, colB)
{
	x[get(colB)==0,(colA):=NA]
	x[get(colB)!=0,(colA):=get(colA)/get(colB)]
	return(x)
}

#' The function grabs each column one by one. The function
#' is applied to the first two columns to create a result
#' (e.g., fun(c1, c2, ret) = ret). This result is then used for combining with the
#' next column (e.g., fun(c2, c3, ret) = ret). The final vector result is then
#' returned. This is usually used in conjunction with lapply.data.table
#' to test values in each column, creating logical columns.
#' 
#' So, if you want to check for each row if ...
#' 
#'  All cols are T --> summarizeLogicalsForEachRow(x, FUN=function(a,b,ret){ret <- ret & a & b})
#'  Any cols are T --> summarizeLogicalsForEachRow(x, FUN=function(a,b,ret){!a & !b)})
#'  Num cols that are T --> summarizeLogicalsForEachRow(x, FUN=function(a,b,ret){a | b})
#'  
#' @param col.fun a function that is first applied to all the specified columns
#' @param row.fun a function that is used after applying the col.fun, it is sequentially applied to each column and takes two vector parameters the first representing a column from the table, the second being the return vector of row.fun that is recycled to accumulate results across columns
#' @param ret.init the initial value provided for the ret parameter to row.fun
#' @param mCols the columns to analyze
#' @param mColsContaining if a col name contains this string, that column will be analyzed, otherwise it is left out
#' @param mColFilter a function that when applied using lapply, returns a vector of logicals indicating which cols to analyze
#' @export
#' 
summarizeRows <- function(x, col.fun=NULL, row.fun, ret.init=NULL, mCols=NULL, mColsContaining=NULL, mColFilter=NULL)
{
     if(is.null(mCols))
     {
          if(!is.null(mColsContaining))
          {
               mCols <- getColNamesContaining(x, mColsContaining)
          }
          else if(!is.null(mColFilter))
          {
               mCols <- as.logical(lapply(x, mColFilter))
               mCols <- names(x)[mCols]
          }
          else
          {
               mCols <- names(x)
          }
     }
     
     # Apply the col.fun
     y <- lapply.data.table(x, FUN=col.fun, cols=mCols, in.place=F)
     
     n <- seq_along(mCols)
     for(i in n)
     {
          ret.init <- row.fun(y[[mCols[i]]], ret.init)
     }
     return(ret.init)
}

allColsTrue <- function(x, test, mCols=NULL, mColsContaining=NULL, mColFilter=NULL)
{
     ret <- summarizeRows(x, col.fun=test, row.fun=function(a, ret){a & ret}, ret.init=rep(T, nrow(x)), mCols=mCols, mColsContaining=mColsContaining, mColFilter=mColFilter)
     return(ret)
}

anyColsTrue <- function(x, test, mCols=NULL, mColsContaining=NULL, mColFilter=NULL)
{
     ret <- summarizeRows(x, col.fun=test, row.fun=function(a, ret){a | ret}, ret.init=rep(F, nrow(x)), mCols=mCols, mColsContaining=mColsContaining, mColFilter=mColFilter)
     return(ret)
}

numColsTrue <- function(x, test, mCols=NULL, mColsContaining=NULL, mColFilter=NULL)
{
     ret <- summarizeRows(x, col.fun=test, row.fun=function(a, ret){ret + a}, ret.init=rep(0, nrow(x)), mCols=mCols, mColsContaining=mColsContaining, mColFilter=mColFilter)
     return(ret)
}

removeColsWithNonFiniteVals <- function(x, cols=NULL)
{
  removeColsMatching(x, cols=cols, col.test=function(n){any(!is.finite(n))})

	# duh <- x[,lapply(.SD, function(y){length(which(!is.finite(y))) > 0}), .SDcols=getNumericCols(x)]
	# duh2 <- getNumericCols(x)[as.logical(as.vector(duh))]
	# if(length(duh2 > 0))
	# {
	# 	print("Removing cols with infinite values...")
	# }
	# for(col in duh2)
	# {
	# 	print(col)
	# 	x[,(col):=NULL]
	# }
}

removeColsMatching <- function(x, cols=NULL, col.test=function(n){all(!is.finite(n))}, ...)
{
  # Remove rows and columns of data contining non-finite data (typically inverses etc.)
  temp.names <- copy(names(x))
  lapply.data.table(x, FUN=function(a){if(col.test(a, ...)){return(NULL)}else{return(a)}}, cols=cols, in.place=T)
  print('Removed the following columns.')
  temp.names <- temp.names[!(temp.names %in% names(x))]
  print(temp.names)
}

getColNamesContaining <- function(x, name)
{
     return(names(x)[grepl(name,names(x),fixed=TRUE)])
}

removeCols <- function(x, colsToRemove)
{
	
     colsToRemove <- colsToRemove[colsToRemove %in% names(x)]
     if(length(colsToRemove) == 0)
     {
          print(paste0("Didn't find any columns to remove."))
     }
     else
     {
          print(paste0("Removing columns"))
          for(colToRemove in colsToRemove)
          {
               print(colToRemove)
               x[,(colToRemove):=NULL]
          }
     }
	return(x)
}

removeColsContaining <- function(x, stringsToMatch)
{
	print(paste0("Removing colums with names containing..."))
     colsToRemove <- c()
	for(stringToMatch in stringsToMatch)
	{
		print(stringToMatch)
		colsToRemove <- c(colsToRemove, getColNamesContaining(x, stringToMatch))
	}
     colsToRemove <- unique(colsToRemove[colsToRemove %in% names(x)])
     if(length(colsToRemove) == 0)
     {
          print(paste0("Didn't find any columns to remove."))
     }
     else
     {
          print(paste0(""))
          print(paste0("Removing colums..."))
          for(colToRemove in colsToRemove)
          {
               print(colToRemove)
               x[,(colToRemove):=NULL]
          }
     }
	return(x)
}

fixColNames <- function(x)
{
	replaceStringInColNames(x, ' ', '')
	replaceStringInColNames(x, '$', '.')
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

getNumericColsOfInterest <- function(x, data.cols=NULL, data.cols.contains=NULL)
{
     ret <- getNumericCols(x)
     if(!is.null(data.cols))
     {
          ret <- intersect(data.cols, ret)
     }
     else if(is.null(data.cols) & !is.null(data.cols.contains))
     {
          ret <- intersect(getColNamesContaining(x, data.cols.contains), ret)
     }
     return(ret)
}

getNonNumericCols <- function(x)
{
	return(names(x)[!unlist(x[,lapply(.SD, is.numeric)])])
}

replaceStringInColNames <- function(x, old, new)
{
	oldNames <- names(x)
	newNames <- gsub(old, new, names(x), fixed=T)
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


# Function for calculating point stats in data.table call.
getPointStats <- function(x, y, weights)
{
     require(spatstat)
     if(length(x) > 1)
     {
          pts <- matrix(c(x, y), ncol=2)
          circ <- getMinCircle(pts)
     }
     else if(length(x) == 1)
     {
          circ <- list(rad=1L, ctr=c(x,y))
     }
     else
     {
          return(lapply(list(Intensity=NA, WeightedIntensity=NA, ConvexArea=NA, ConvexPerimeter=NA), as.double))
     }
     pp <- ppp(x, y, window=disc(radius=circ$rad*1.01, centre=circ$ctr))
     hull <- convexhull(pp)
     ret <- lapply(list(Intensity=intensity(pp), WeightedIntensity=intensity(pp, weights=weights), ConvexArea=area.owin(hull), ConvexPerimeter=perimeter(hull), Diameter=circ$rad), as.double)
     if(ret$ConvexArea == 0 || ret$ConvexPerimeter == 0)
     {
          ret$Circularity <- 0
     }
     else
     {
          ret$Circularity <- as.double(4*pi*ret$ConvexArea/(ret$ConvexPerimeter*ret$ConvexPerimeter))
     }
     return(ret)
}

#' Takes a table with cId, ImageChannel, and MaskChannels (+ feature columns)
#' Geometry.<X> features have values for sub-regions of the masks so each
#' MaskChannel value has encoded the subregion number as well.
#' 
#' This function collapses the subregion data to just single values for each cell.
#' 
summarizeGeometry <- function(x, cellIdCols='cId', removeXY=T)
{
     idCols <- cellIdCols
     
     # If we haven't run this function before
     if(!('Geometric.MaximumFeretsDiameter' %in% names(x)))
     {
          stop("It looks like you already ran this function on the table given it is missing columns that are deleted by this function. Aborting.")
     }
     
     # Gather important columns for segregating the data for calculations
     colsToSummarize <- getColNamesContaining(x, 'Geometric.')
     colsToKeep <- getAllColNamesExcept(x, colsToSummarize)
     
     # Table to keep
     tableToKeep <- x[, mget(colsToKeep)]
     rowsToDiscard <- allColsTrue(tableToKeep, test=is.na, mCols=getAllColNamesExcept(tableToKeep, c(idCols,'ImageChannel','MaskChannel')))
     tableToKeep <- tableToKeep[!rowsToDiscard]
     
     # Table to summarize
     x <- x[, mget(colsToSummarize), by=c(idCols,'MaskChannel','ImageChannel')] # Use by statement to keep idcols
     rowsToDiscard <- allColsTrue(x, test=is.na, mCols=getAllColNamesExcept(x, c(idCols,'ImageChannel','MaskChannel')))
     x <- x[!rowsToDiscard]
     
     # If there isn't one already, create a column that just has the MaskChannel information split from the subregion ids (i.e., p1, p2, ..., pn)
     if(!('MaskChannel2' %in% names(x)))
     {
          splitColumnAtString(x, colToSplit='MaskChannel', sep='.p', newColNames=c('MaskChannel2'), keep=c(1L))
     }
     
     # # Check to make sure that eventual remerging of the tables will be legitimate.
     # uniqueKeep <- tableToKeep[, list(v=T), by=c('ImageChannel','MaskChannel2')]
     # uniqueSummary <- x[, list(v=T), by=c('ImageChannel','MaskChannel2')]
     # testTable <- rbindlist(list(uniqueKeep,uniqueSummary), use.names=T)
     # testTable <- testTable[, list(v=T), by=c('ImageChannel','MaskChannel')]
     # if((nrow(uniqueKeep) + nrow(uniqueSummary)) != nrow(testTable))
     # {
     #      stop("We use rbindlist to merge the tables later which assumes unique cId, 
     #           ImageChannel, MaskChannel combinations for geometric data compared to 
     #           other data. This doesn't appear to be true so throwing error and aborting. 
     #           Need to use merge function for portions of the table that share cId, 
     #           ImageChannel, MaskChannel combinations and rbindlist for unique ones.")
     # }
     
     # Make columns with weights for each subregion
     x[, ':='(weights=Geometric.SizeIterable/(sum(Geometric.SizeIterable, na.rm=T)), countWeights=Geometric.SizeIterable/(max(Geometric.SizeIterable, na.rm=T))), by=c(idCols, 'ImageChannel', 'MaskChannel2')]
     
     # Calculate ratios and remove the corresponding parent metrics
     x[, ':='(Geometric.FeretsAspectRatio = Geometric.MaximumFeretsDiameter/Geometric.MinimumFeretsDiameter, Geometric.EllipseAspectRatio = Geometric.MajorAxis/Geometric.MinorAxis)]
     removeCols(x, c('Geometric.MaximumFeretsDiameter', 'Geometric.MinimumFeretsDiameter', 'Geometric.MajorAxis', 'Geometric.MinorAxis', 'Geometric.MaximumFeretsAngle', 'Geometric.MinimumFeretsAngle'))
     
     # Remove subregion data where any Geometric feature measure is NA (typically small regions with no area etc.)
     rowsToDiscard <- anyColsTrue(x, test=is.na, mColsContaining='Geometric.')
     x <- x[!rowsToDiscard]
     # Now we have only subregion data where the geometric information is fully defined
     
     # Decide how to combine different geometric features
     #geomFeatures <- c('Convexity', 'Solidity', 'SizeIterable', 'BoundarySize', 'MainElongation', 'Circularity',
     #'Boxivity', 'Eccentricity', 'MajorAxis', 'MaximumFeretsDiameter', 'MinimumFeretsDiameter',
     #'MinorAxis', 'Roundness', 'X', 'Y')
     geomFeatures_Total <- c('Geometric.SizeIterable', 'Geometric.BoundarySize')
     geomFeatures_SizeWeightedMean <- c('Geometric.SizeIterable', 'Geometric.Convexity', 'Geometric.Solidity', 'Geometric.MainElongation', 'Geometric.Circularity', 'Geometric.Boxivity', 'Geometric.Eccentricity', 'Geometric.BoundarySize','Geometric.FeretsAspectRatio','Geometric.EllipseAspectRatio','Geometric.Roundness')
     
     # Aggregate geometry data from the different subregions (this created duplicate row information)
     # all na.rm's are removed to cause errors if any NA's are found.
     for(feature in geomFeatures_Total)
     {
          x[, (paste0(feature, '.Total')):=sum(get(feature)), by=c(idCols, 'ImageChannel', 'MaskChannel2')]
     }
     for(feature in geomFeatures_SizeWeightedMean)
     {
          x[, (feature):=sum(get(feature)*weights), by=c(idCols, 'ImageChannel', 'MaskChannel2')]
     }
     
     # Create new count columns (this also results in duplicate row information)
     x[, ':='(N=as.double(.N), weightedN=sum(countWeights)), by=c(idCols, 'ImageChannel', 'MaskChannel2')]
     
     # Calculate point stats
     x[, c('Geometric.SubRegionIntensity', 'Geometric.SubRegionWeightedIntensity', 'Geometric.SubRegionConvexArea', 'Geometric.SubRegionConvexPerimeter', 'Geometric.SubRegionRadius', 'Geometric.SubRegionCircularity'):=getPointStats(Geometric.X, Geometric.Y, weights), by=c(idCols, 'ImageChannel', 'MaskChannel2')]
     
     # Overwrite Maskchannel that currently encodes subregion as well and replace with just MaskChannel information.
     # Also remove helper columns
     if(removeXY)
     {
          x[, ':='(MaskChannel=MaskChannel2, MaskChannel2=NULL, weights=NULL, countWeights=NULL, Geometric.X=NULL, Geometric.Y=NULL)]
     }
     else
     {
          x[, ':='(MaskChannel=MaskChannel2, MaskChannel2=NULL, weights=NULL, countWeights=NULL)]
     }
     
     
     # Remove all the duplicate information that was created during calculations
     x <- unique(x)
     
     # Now merge tableToKeep and x
     x <- merge(tableToKeep, x, all=T)
     
     # Return the result
     return(x)
}

standardizeWideData <- function(x, row.normalize=F, row.use.median=F, col.use.median=F, col.use.mad=F, data.cols=NULL, data.cols.contains=NULL, by=NULL, trySDIfNeeded=T)
{
     # x <- data.table(a=1.1:3.1, b=4.1:6.1, c=c(100.1,110.1,120.1)); duh <- copy(x); duh2 <- as.data.table(lapply(x, log))
     # data.cols <- c('b','c')
     # row.normalize=F
     # row.use.median=F
     # col.use.median=F
     # col.use.mad=F
     # data.cols.contains = NULL
     # by=NULL

     # Get the names of the columns of interest
     cols <- getNumericColsOfInterest(x, data.cols=data.cols, data.cols.contains=data.cols.contains)

     # Remove cols with now variance
     removeNoVarianceCols(x, use.mad=col.use.mad, cols=cols, by=by, trySDIfNeeded=trySDIfNeeded)

     # Redefine cols just in case
     cols <- getNumericColsOfInterest(x, data.cols=data.cols, data.cols.contains=data.cols.contains)

     if(row.normalize)
     {
          # Take each row and divide it by the row mean/median (i.e., row normalize)
          x[, rowNum:=as.double(1:nrow(x))] # create a rowNum column to do operations by row
          if(row.use.median)
          {
               x[, c(cols) := lapply(.SD, "-", median(unlist(.SD))), .SDcols=cols, by=.(rowNum)]
          }else
          {
               x[, c(cols) := lapply(.SD, "-", mean(unlist(.SD))), .SDcols=cols, by=.(rowNum)]
          }
          x[, rowNum:=NULL] # Remove the rowNum column
     }

     x[,c(cols):=lapply(.SD, robustScale, use.median=col.use.median, use.mad=col.use.mad, trySDIfNeeded=trySDIfNeeded), .SDcols=cols, by=by]
     return(x)
}

# Now perform regular column standardization (grouping as appropriate using 'by')
robustScale <- function(x, use.median, use.mad, trySDIfNeeded=T)
{
     if(use.median)
     {
          m <- median(x[is.finite(x)], na.rm=TRUE)
     }
     else
     {
          m <- mean(x[is.finite(x)], na.rm=TRUE)
     }
     
     if(use.mad)
     {
          sig <- mad(x[is.finite(x)], center=m, na.rm=TRUE)
          if(!is.na(sig) & sig == 0 & trySDIfNeeded)
          {
               sig <- sd(x[is.finite(x)], na.rm=TRUE)
          }
     }
     else
     {
          sig <- sd(x[is.finite(x)], na.rm=TRUE)
     }
     
     return((x-m)/sig)
}

removeNoVarianceCols <- function(x, use.mad=F, cols=NULL, by=NULL, trySDIfNeeded=T)
{
     namesToRemove <- getNoVarianceCols(x, use.mad=use.mad, cols=cols, by=by, trySDIfNeeded=trySDIfNeeded)
     if(length(namesToRemove) > 0)
     {
          if(use.mad)
          {
               print("Removing cols with a MAD of zero...")
          }
          else
          {
               print("Removing cols with a variance of zero...")
          }
          for(name in namesToRemove)
          {
               print(name)
               x[,(name):=NULL]
          }
     }
}

tempSD1 <- function(y, trySD)
{
     temp <- mad(y[is.finite(y)], na.rm = TRUE)
     if(!is.na(temp) & temp == 0 & trySD)
     {
          temp <- sd(y[is.finite(y)], na.rm = TRUE)
     }
     return(temp)
}

tempSD2 <- function(y, trySD)
{
     return(sd(y, na.rm = TRUE))
}

getNoVarianceCols <- function(x, use.mad, cols=NULL, by=NULL, trySDIfNeeded=T)
{
     cols <- getNumericColsOfInterest(x, data.cols=cols)
     if(use.mad)
     {
          tempSD <- tempSD1
          if(is.null(by))
          {
               tempNames <- x[,lapply(.SD, tempSD, trySD=trySDIfNeeded), .SDcols=cols]
          }
          else
          {
               tempNames <- x[,lapply(.SD, tempSD, trySD=trySDIfNeeded), .SDcols=cols, by=by]
               tempNames <- tempNames[, lapply(.SD, min), .SDcols=cols]
          }
     }
     else
     {
          tempSD <- tempSD2
          if(is.null(by))
          {
               tempNames <- x[,lapply(.SD, tempSD, trySD=trySDIfNeeded), .SDcols=cols]
          }
          else
          {
               tempNames <- x[,lapply(.SD, tempSD, trySD=trySDIfNeeded), .SDcols=cols, by=by]
               tempNames <- tempNames[, lapply(.SD, min), .SDcols=cols]
          }
     }
     return(names(tempNames)[as.numeric(as.vector(tempNames))==0])
}

##### Long Table Operations #####

divideMAbyMBbyRef <- function(x, mA, mB)
{
	mATable <- x[Measurement==mA]
	mBTable <- x[Measurement==mB]
	if(nrow(mATable) != nrow(mBTable))
	{
		# Try to perform the operation on the subset of the mB column (can't do reverse because we are editing the mA column)
		mBTable <- mBTable[MaskChannel %in% unique(mATable$MaskChannel)]
		if(nrow(mATable) != nrow(mBTable))
		{
			stop('Number of rows for these measurements do not match! Aborting operation.')
		}
	}
	ret <- mATable$Value / mBTable$Value
	x[Measurement==mA]$Value <- ret
	return(x)
}

# intIntensityNormalizeCentralMoments <- function(x)
# {
# 	mNames <- getMeasurementNamesContaining(x, 'ImageMoments.CentralMoment')
# 	for(mName in mNames)
# 	{
# 		x <- divideMAbyMBbyRef(x, mName, 'Stats.Sum')
# 	}
# 	return(x)
# }

# meanNormalizeZernikeMoments <- function(x)
# {
# 	mNames <- getMeasurementNamesContaining(x, 'ZernikeMag')
# 	for(mName in mNames)
# 	{
# 		x <- divideMAbyMBbyRef(x, mName, 'Stats.Mean')
# 	}
# 	return(x)
# }

getRowsMatching <- function(x, col, baseName)
{
	return(x[grepl(baseName, x[[col]])])
}

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

getMomentTable <- function(x, baseName='ImageMoments.CentralMoment')
{
	theNames <- unique(x[['Measurement']])
	theNames <- theNames[grepl(baseName,theNames, fixed=TRUE)]
	start <- nchar(baseName)
	orders <- substr(theNames, start+1, start+2)
	ret <- data.frame(Measurement=theNames, orderx=as.numeric(substr(orders,1,1)), ordery=as.numeric(substr(orders,2,2)))
	return(ret)
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

removeMeasurementNames <- function(x, names)
{
     print("Removing the following Measurements...")
     for(name in names)
     {
          print(name)
     }
     x <- x[!(Measurement %in% names)]
     return(x)
}

standardizeLongData <- function(x, by=c('MaskChannel','ImageChannel','Measurement','Expt'))
{
	robustScale <- function(x, measurement)
	{
		if(substr(measurement,1,12) == 'ZernikePhase')
		{
			return(x)
		}
		else
		{
			m <- median(x, na.rm=TRUE)
			return((x-m)/mad(x, center=m, na.rm=TRUE))
		}
	}
	x <- removeNoMADMeasurements(x, by=by)
	x[,Value:=robustScale(Value,Measurement),by=by]
	return(x)
}

standardizeLongData2 <- function(x, val.col='Value', by)
{
	robustScale <- function(x)
	{
		m <- median(x, na.rm=TRUE)
		return((x-m)/mad(x, center=m, na.rm=TRUE))
	}
	x[,c(val.col):=robustScale(get(val.col)),by=by]
	return(x)
}

removeNoVarianceMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'))
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

removeNoMADMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'))
{
	# Tempororarily add a column in the table with stdev in it
	x[,MAD:=mad(get(val), na.rm=TRUE), by=by]
	toRemove <- unique(x[MAD == 0]$Measurement)
	if(length(toRemove)>0)
	{
		print("Removing measurements with 0 MAD...")
		for(m in toRemove)
		{
			print(m)
		}
		y <- x[!(Measurement %in% toRemove)]
		x[, MAD:=NULL]
		y[, MAD:=NULL]
		return(y)
	}else
	{
		x[, MAD:=NULL]
		return(x)
	}
}

# removeNoVarianceMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'))
# {
# 	# See if we have any columns to remove and record the info for reporting
# 	temp <- x[,list(stdev=sd(get(val))), by=by]
# 	temp <- data.frame(temp[stdev == 0])
# 	print("Removing measurements with 0 variance...")
# 	print(temp)
# 	# Tempororarily add a column in the table with stdev in it
# 	x[,stdev:=sd(get(val)), by=by]
# 	y <- x[stdev != 0]
# 	x[, stdev:=NULL]
# 	y[, stdev:=NULL]
# 	return(y)
# }

replaceSubStringInAllRowsOfCol <- function(x, old, new, col)
{
	x[,c(col):=gsub(old,new,get(col),fixed=TRUE)]
}

trySample <- function(x, n, replace=F, prob=NULL)
{
	if(n > length(x))
	{
		return(x)
	}
	else
	{
		return(sample(x, n, replace, prob))
	}
}

fixLongTableStringsInCol <- function(x, col)
{
	replaceSubStringInAllRowsOfCol(x,'_Order_','',col)
	replaceSubStringInAllRowsOfCol(x,'_Rep_','',col)
	replaceSubStringInAllRowsOfCol(x,'$','.',col)
	replaceSubStringInAllRowsOfCol(x,'net.imagej.ops.Ops.','',col)
	replaceSubStringInAllRowsOfCol(x,'function.ops.JEXOps.','',col)
	replaceSubStringInAllRowsOfCol(x,' ','',col)
	replaceSubStringInAllRowsOfCol(x,':','_',col)
}

##### Feature Calculations #####

filterLBPCodes <- function(x, nSigma=-2)
{
     LBPCols <- getColNamesContaining(x, 'LBP')
     LBPTots <- as.numeric(as.vector(lapply.data.table(x, FUN=function(x){log(sum(x, na.rm=T)+1)}, cols=LBPCols)))
     LBPColsToDelete <- LBPCols[LBPTots < (median(LBPTots) + nSigma*mad(LBPTots))]
     removeCols(x, LBPColsToDelete)
}

splitColumnAtString <- function(x, colToSplit, newColNames, sep='.', keep=NULL)
{
     x[, (newColNames) := tstrsplit(get(colToSplit), sep, fixed=T, keep=keep)]
}

calculateChannelProducts <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN=getComboProducts, sep='_')
{
     return(performChannelComboCalcs(x, comboCol=comboCol, valsToPermute=valsToPermute, idCols=idCols, mCols=mCols, mColContains=mColContains, FUN=FUN, sep=sep))
}

calculateChannelDifferences <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN=getComboDifferences, sep='_')
{
     return(performChannelComboCalcs(x, comboCol=comboCol, valsToPermute=valsToPermute, idCols=idCols, mCols=mCols, mColContains=mColContains, FUN=FUN, sep=sep))
}

calculateChannelRatios <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN=getComboRatios, sep='_')
{
     return(performChannelComboCalcs(x, comboCol=comboCol, valsToPermute=valsToPermute, idCols=idCols, mCols=mCols, mColContains=mColContains, FUN=FUN, sep=sep))
}

calculateChannelLogRatios <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN=getComboLogRatios, sep='_')
{
     return(performChannelComboCalcs(x, comboCol=comboCol, valsToPermute=valsToPermute, idCols=idCols, mCols=mCols, mColContains=mColContains, FUN=FUN, sep=sep))
}

#' Intended to be used with a wide table that still has ImageChannel and MaskChannel information.
performChannelComboCalcs <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN, sep='_', ...)
{
     uniqueChannels <- valsToPermute
     
     if(is.null(mCols))
     {
          mCols <- getColNamesContaining(x, mColContains)
     }
     
	if(length(uniqueChannels) > 1)
	{
	     validRows <- x[[comboCol]] %in% uniqueChannels
	     for(mCol in mCols)
	     {
	          validRows <- validRows & !is.na(x[[mCol]])
	     }
	     x2 <- x[validRows, lapply(.SD, FUN=FUN, ...), .SDcols=mCols, by=idCols]
	     x2[[comboCol]] <- x[validRows, lapply(.SD, FUN=getComboNames, sep=sep), .SDcols=comboCol, by=idCols][[comboCol]]
	    return(x2)
	     
	}else
	{
		# return an empty table with the same columns as provided
		return(x[FALSE])
	}
}

mergeComboData <- function(x, comboData, mCol.old, mCol.new)
{
     setnames(comboData, old=mCol.old, new=mCol.new)
     x <- merge(x, comboData, by=c('cId','ImageChannel','MaskChannel'), all=T)
     return(x)
}

# Meant to be called on a subset of the main table
getComboNames <- function(x, sep='_')
{
	if(length(x) < 2)
	{
		return(NULL)
	}
	temp <- combn(x, 2)
	#print(temp)
	temp <- paste0(temp[1,],sep,temp[2,])
	return(temp)
}

# Meant to be called on a subset of the main table
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

# Meant to be called on a subset of the main table
getComboProducts <- function(x)
{
	if(length(x) < 2)
	{
		return(NULL)
	}
	temp <- combn(x, 2)
	temp <- temp[1,]*temp[2,]
	return(temp)
}

# Meant to be called on a subset of the main table
getComboRatios <- function(x)
{
     if(length(x) < 2)
     {
          return(NULL)
     }
     temp <- combn(x, 2)
     temp <- temp[1,]/temp[2,]
     return(temp)
}

# Meant to be called on a subset of the main table
getComboLogRatios <- function(x)
{
     if(length(x) < 2)
     {
          return(NULL)
     }
     temp <- combn(x, 2)
     temp <- log(temp[1,]/temp[2,])
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

calculateGeometricFeaturesFromDNZernikeCircles <- function(x)
{
     if(sum(grepl('DNZernike',names(x),fixed=T)) == 0)
     {
          warning("Didn't find any DNZernike circles to perform calculations. Aborting.")
     }
     else
     {
          # the outercircle is 1.5 times the equivalent radius of the whole cell, so adjust appropriately to get actual cell radius
          x[, ':='(Geometric.Area.Nuc=pi*DNZernikeInnerCircleR^2, Geometric.Area.WholeCell=pi*((2.0/3.0)*DNZernikeOuterCircleR)^2)]
          x[, ':='(Geometric.NucAreaFraction=Geometric.Area.Nuc/Geometric.Area.WholeCell)]
          x[, ':='(Geometric.NucNormalizedOffset=sqrt((DNZernikeInnerCircleX-DNZernikeOuterCircleX)^2 + (DNZernikeInnerCircleY-DNZernikeOuterCircleY)^2)/DNZernikeInnerCircleR)]
          x[, ':='(Geometric.WholeCellNormalizedOffset=sqrt((DNZernikeInnerCircleX-DNZernikeOuterCircleX)^2 + (DNZernikeInnerCircleY-DNZernikeOuterCircleY)^2)/DNZernikeOuterCircleR)]
          x[, ':='(Geometric.NormalizedOffset=sqrt((DNZernikeInnerCircleX-DNZernikeOuterCircleX)^2 + (DNZernikeInnerCircleY-DNZernikeOuterCircleY)^2)/(DNZernikeOuterCircleR-DNZernikeInnerCircleR))]
     }
}

normalizeColsToOtherCol <- function(x, numeratorCols, denominatorCol='Stats.Sum')
{
     FUN <- function(a,b)
     {
          ret <- a/b
          ret[!is.finite(ret)] <- NA
          return(ret)
     }
     x[, c(numeratorCols):=lapply(.SD, FUN=FUN, get(denominatorCol)), .SDcols=numeratorCols][]
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

