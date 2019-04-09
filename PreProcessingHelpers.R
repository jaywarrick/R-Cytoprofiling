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

##### Wide Table Operations #####

divideManyColsByOneCol <- function(x, numeratorCols, denominatorCol='Stats.Sum')
{
	FUN <- function(a,b)
	{
		ret <- a/b
		ret[!is.finite(ret)] <- NA
		return(ret)
	}
	x[, c(numeratorCols):=lapply(.SD, FUN=FUN, get(denominatorCol)), .SDcols=numeratorCols][]
}

divideColAByColB <- function(x, colA, colB)
{
	return(divideManyColsByOneCol(x, numeratorCols = colA, denominatorCol=colB))
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

removeColsWithAnyNonFiniteVals <- function(x, cols=NULL)
{
	removeColsMatching(x, cols=cols, col.test=function(n){any(!is.finite(n))})
}

removeColsWithAllNonFiniteVals <- function(x, cols=NULL)
{
	removeColsMatching(x, cols=cols, col.test=function(n){all(!is.finite(n))})
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

getColNamesContaining <- function(x, names, and=T)
{
	matchingCols <- rep(T, length(names(x)))
	for(stringToMatch in names)
	{
		if(and)
		{
			matchingCols <- matchingCols & grepl(stringToMatch,names(x),fixed=TRUE)
		}
		else
		{
			matchingCols <- matchingCols | grepl(stringToMatch,names(x),fixed=TRUE)
		}
	}
	matchingNames <- names(x)[matchingCols]
	if(length(matchingCols) == 0)
	{
		print(paste0("Didn't find any matching columns!!!"))
	}
	else
	{
		return(matchingNames)
	}
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
	replaceSubstringInColNames(x,'_Order_','')
	replaceSubstringInColNames(x,'_Rep_','')
	replaceSubstringInColNames(x,'$','.')
	replaceSubstringInColNames(x,'net.imagej.ops.Ops.','')
	replaceSubstringInColNames(x,'function.ops.JEXOps.','')
	replaceSubstringInColNames(x,' ','')
	replaceSubstringInColNames(x,':','_')
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

replaceSubstringInColNames <- function(x, old, new)
{
	oldNames <- names(x)
	newNames <- gsub(old, new, names(x), fixed=T)
	setnames(x, oldNames, newNames)
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
#' Note: The names must be 'fixed' as well using either 'fixColNames' or 'fixLongTableStringsInCol'
#' 
#' This function collapses the subregion data to just single values for each cell.
#' 
summarizeGeometry <- function(x, cellIdCols='cId', removeXY=T)
{
     idCols <- cellIdCols
     
     # If we haven't run this function before
     if(!('Geometric.SizeIterable' %in% names(x)))
     {
          stop("It looks like you already ran this function on the table given it is missing columns that are deleted by this function. Aborting.")
     }
     
     # Gather important columns for segregating the data for calculations
     colsToSummarize <- getColNamesContaining(x, 'Geometric.')
     colsToSummarize <- colsToSummarize[!(colsToSummarize %in% c('Geometric.COMX','Geometric.COMY','Geometric.MaximaX','Geometric.MaximaY'))]
     colsToKeep <- getAllColNamesExcept(x, colsToSummarize)
     
     # Table to keep
     tableToKeep <- x[, mget(colsToKeep)]
     rowsToDiscard <- allColsTrue(tableToKeep, test=is.na, mCols=getAllColNamesExcept(tableToKeep, c(idCols,'ImageChannel','MaskChannel')))
     tableToKeep <- tableToKeep[!rowsToDiscard]
     if(any(grepl('.p', tableToKeep$MaskChannel, fixed=T)))
     {
     	stop("There shouldn't be any .p nomenclature left over in the MaskChannel column at this step within the function. Check to see if all of the id columns were provided. Aborting.")
     }
     
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
     if(all(c('Geometric.MaximumFeretsDiameter', 'Geometric.MinimumFeretsDiameter') %in% names(x)))
     {
     	x[, ':='(Geometric.FeretsAspectRatio = Geometric.MaximumFeretsDiameter/Geometric.MinimumFeretsDiameter)]
     }
     if(all(c('Geometric.MajorAxis', 'Geometric.MinorAxis') %in% names(x)))
     {
     	x[, ':='(Geometric.EllipseAspectRatio = Geometric.MajorAxis/Geometric.MinorAxis)]
     }
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
     geomFeatures_Total <- geomFeatures_Total[geomFeatures_Total %in% names(x)]
     geomFeatures_SizeWeightedMean <- c('Geometric.SizeIterable', 'Geometric.Convexity', 'Geometric.Solidity', 'Geometric.MainElongation', 'Geometric.Circularity', 'Geometric.Boxivity', 'Geometric.Eccentricity', 'Geometric.BoundarySize','Geometric.FeretsAspectRatio','Geometric.EllipseAspectRatio','Geometric.Roundness')
     geomFeatures_SizeWeightedMean <- geomFeatures_SizeWeightedMean[geomFeatures_SizeWeightedMean %in% names(x)]
     
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
     
     if(all(c('Geometric.X','Geometric.Y') %in% names(x)))
     {
       # Calculate point stats
       x[, c('Geometric.SubRegionIntensity', 'Geometric.SubRegionWeightedIntensity', 'Geometric.SubRegionConvexArea', 'Geometric.SubRegionConvexPerimeter', 'Geometric.SubRegionRadius', 'Geometric.SubRegionCircularity'):=getPointStats(Geometric.X, Geometric.Y, weights), by=c(idCols, 'ImageChannel', 'MaskChannel2')]
     }
     
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
     x <- merge(x=tableToKeep, y=x, by=c(idCols, 'MaskChannel', 'ImageChannel'), all=T)
     
     # Return the result
     return(x)
}

summarizeSymmetryData <- function(x, sim.trans=T, logit.trans=T)
{
	# Summarize Symmetry Data
	toDelete <- getColNamesContaining(x,'SymmetryCorrelation')
	setnames(x, 'SymmetryCorrelation.1.1', 'SymmetryCorrelation.1.Avg')
	sym2 <- getColNamesContaining(x,'SymmetryCorrelation.2')
	sym3 <- getColNamesContaining(x,'SymmetryCorrelation.3')
	sym4 <- getColNamesContaining(x,'SymmetryCorrelation.4')
	x[, SymmetryCorrelation.2.Avg:=apply(.SD,1,mean,na.rm=T), .SDcols=sym2]
	x[, SymmetryCorrelation.3.Avg:=apply(.SD,1,mean,na.rm=T), .SDcols=sym3]
	x[, SymmetryCorrelation.4.Avg:=apply(.SD,1,mean,na.rm=T), .SDcols=sym4]
	x[, c(sym2,sym3,sym4):=NULL]
	x[, SymmetryLobe.1:=SymmetryAmplitude.1*SymmetryCorrelation.1.Avg]
	x[, SymmetryLobe.2:=SymmetryAmplitude.2*SymmetryCorrelation.2.Avg]
	x[, SymmetryLobe.3:=SymmetryAmplitude.3*SymmetryCorrelation.3.Avg]
	x[, SymmetryLobe.4:=SymmetryAmplitude.4*SymmetryCorrelation.4.Avg]
	correlationNames <- c(getColNamesContaining(x,'SymmetryCorrelation'), getColNamesContaining(x,'SymmetryLobe'))
	if(length(correlationNames) > 0 & sim.trans)
	{
		lapply.data.table(x, FUN=sim.transform, cols=correlationNames, in.place=T)
		setnames(x, correlationNames, paste0(correlationNames, '.Similarity'))
	}
	ampNames <- getColNamesContaining(x, 'SymmetryAmplitude')
	if(length(ampNames) > 0 & logit.trans)
	{
		# Use logit transform for amplitudes because from 0 to 1
		lapply.data.table(x, FUN=logit.transform, cols=ampNames, in.place=T)
		setnames(x, ampNames, paste0(ampNames, '.Logit'))
	}
	return(x)
}

standardizeWideData <- function(x, row.normalize=F, row.use.median=F, col.use.median=F, col.use.mad=F, col.use.percentiles=F, percentiles=c(0.05,0.95), data.cols=NULL, data.cols.contains=NULL, by=NULL, trySDIfNeeded=T, na.rm.no.variance.cols=F, suffix=NULL)
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
	removeNoVarianceCols(x, use.mad=col.use.mad, cols=cols, by=by, trySDIfNeeded=trySDIfNeeded, na.rm=na.rm.no.variance.cols)
	
	# Redefine cols just in case
	cols <- getNumericColsOfInterest(x, data.cols=data.cols, data.cols.contains=data.cols.contains)
	
	# New cols
	if(!is.null(suffix))
	{
		newCols <- paste(cols, suffix, sep='.')
	}
	
	if(row.normalize)
	{
		# Take each row and divide it by the row mean/median (i.e., row normalize)
		x[, rowNum:=as.double(1:nrow(x))] # create a rowNum column to do operations by row
		if(row.use.median)
		{
			x[, c(newCols) := lapply(.SD, "-", median(unlist(.SD))), .SDcols=cols, by=.(rowNum)]
		}else
		{
			x[, c(newCols) := lapply(.SD, "-", mean(unlist(.SD))), .SDcols=cols, by=.(rowNum)]
		}
		x[, rowNum:=NULL] # Remove the rowNum column
	}
	
	if(row.normalize)
	{
		x[,c(newCols):=lapply(.SD, FUN=robustScale, use.median=col.use.median, use.mad=col.use.mad, use.percentiles=col.use.percentiles, percentiles=percentiles, trySDIfNeeded=trySDIfNeeded), .SDcols=newCols, by=by]
	}
	else
	{
		x[,c(newCols):=lapply(.SD, FUN=robustScale, use.median=col.use.median, use.mad=col.use.mad, use.percentiles=col.use.percentiles, percentiles=percentiles, trySDIfNeeded=trySDIfNeeded), .SDcols=cols, by=by]
	}
	return(x)
}

# Now perform regular column standardization (grouping as appropriate using 'by')
robustScale <- function(x, use.median, use.mad, use.percentiles=F, percentiles=c(0.05,0.95), trySDIfNeeded=T)
{
	if(use.median)
	{
		m <- median(x[is.finite(x)], na.rm=TRUE)
	}
	else if(use.percentiles)
	{
		l(lo, hi) %=% getPercentileValues(x[is.finite(x)], c(percentiles[1],percentiles[2]))
		l(daMin, daMax) %=% qnorm(percentiles)
		return(adjustIntensity(x, lo, hi, daMin, daMax))
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

removeNoVarianceCols <- function(x, use.mad=F, cols=NULL, by=NULL, trySDIfNeeded=T, na.rm=F)
{
	namesToRemove <- getNoVarianceCols(x, use.mad=use.mad, cols=cols, by=by, trySDIfNeeded=trySDIfNeeded, na.rm=na.rm)
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
	return(sd(y[is.finite(y)], na.rm = TRUE))
}

#'@param na.rm is whether or not to keep cols with groups for which an SD or MAD cannot be calculated (e.g., only 1 value or no values).
#'This is different than if the SD or MAD calculates to 0. If a 0 is encountered for a group, the col will be removed.
#'
getNoVarianceCols <- function(x, use.mad, cols=NULL, by=NULL, trySDIfNeeded=T, na.rm=F)
{
	cols <- getNumericColsOfInterest(x, data.cols=cols)
	tempMin <- function(x, na.rm)
	{
		ret <- min(x, na.rm=na.rm)
		if(is.na(ret))
		{
			return(0)
		}
		return(ret)
	}
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
			tempNames <- tempNames[, lapply(.SD, tempMin, na.rm=na.rm), .SDcols=cols]
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
			tempNames <- tempNames[, lapply(.SD, tempMin, na.rm=na.rm), .SDcols=cols]
		}
	}
	return(names(tempNames)[as.numeric(as.vector(tempNames))==0])
}

##### Long Table Operations #####

#' This is verified to work with randomly ordered tables.
divideMAbyMBbyRef <- function(x, mA, mB, idCols='cId')
{
	tempDivide <- function(values, Measurement, mA, mB)
	{
		a <- values[Measurement==mA]
		b <- values[Measurement==mB]
		if(any(!is.finite(c(a,b))))
		{
			return(NA)
		}
		if(length(a) != 1 || length(b) != 1)
		{
			stop("Id cols were specific enough to produce one value for each id")
		}
		else
		{
			if(Measurement[1]==mA)
			{
				return(c(a/b, b))	
			}
			else
			{
				return(c(b, a/b))
			}
		}
	}
	x[Measurement %in% c(mA,mB), Value:=tempDivide(Value, Measurement, mA, mB), by=idCols]
	# mATable <- x[Measurement==mA]
	# mBTable <- x[Measurement==mB]
	# if(nrow(mATable) != nrow(mBTable))
	# {
	# 	# Try to perform the operation on the subset of the mB column (can't do reverse because we are editing the mA column)
	# 	mBTable <- mBTable[MaskChannel %in% unique(mATable$MaskChannel)]
	# 	if(nrow(mATable) != nrow(mBTable))
	# 	{
	# 		stop('Number of rows for these measurements do not match! Aborting operation.')
	# 	}
	# }
	# ret <- mATable$Value / mBTable$Value
	# x[Measurement==mA]$Value <- ret
	# return(x)
}

#' JEX outputs the normalized (size scale invariant) central moments already
#' To make them intensity invariant, we just need to divide by the mean
#' Not sure why we were normalizing by the total.
#' Plus, the zernike's already give us a rotationally invariant measure and can be
#' made to be contrast invariant by dividing by the mean
#' so, tending to avoid use of these basic moments now since they are directional.
#' If we want rotationally invariant, then we have Hu's.
# intIntensityNormalizeCentralMoments <- function(x, idCols='cId')
# {
# 	mNames <- getMeasurementNamesContaining(x, 'ImageMoments.CentralMoment')
# 	for(mName in mNames)
# 	{
# 		x <- divideMAbyMBbyRef(x, mName, 'Stats.Sum', idCols=idCols)
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

getMeasurementNamesContaining <- function(x, names)
{
	ms <- unique(x$Measurement)
	ret <- character(0)
	for(name in names)
	{
		ret <- c(ret, ms[grepl(name,ms,fixed=T)])
	}
	ret <- unique(ret)
	return(ret)
}

removeMeasurementNamesContaining <- function(x, names)
{
	namesToRemove <- getMeasurementNamesContaining(x, names)
	return(removeMeasurementNames(x, namesToRemove))
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

replaceSubStringInAllRowsOfCol <- function(x, old, new, col)
{
	x[,c(col):=gsub(old,new,get(col),fixed=TRUE)]
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

#' Intended to be used with a wide table that still has ImageChannel and MaskChannel information.
calculateChannelProducts <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN=getComboProducts, sep='_')
{
	return(performChannelComboCalcs(x, comboCol=comboCol, valsToPermute=valsToPermute, idCols=idCols, mCols=mCols, mColContains=mColContains, FUN=FUN, sep=sep))
}

#' Intended to be used with a wide table that still has ImageChannel and MaskChannel information.
calculateChannelDifferences <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN=getComboDifferences, sep='_')
{
	return(performChannelComboCalcs(x, comboCol=comboCol, valsToPermute=valsToPermute, idCols=idCols, mCols=mCols, mColContains=mColContains, FUN=FUN, sep=sep))
}

#' Intended to be used with a wide table that still has ImageChannel and MaskChannel information.
calculateChannelRatios <- function(x, comboCol, valsToPermute, idCols, mCols=NULL, mColContains=NULL, FUN=getComboRatios, sep='_')
{
	return(performChannelComboCalcs(x, comboCol=comboCol, valsToPermute=valsToPermute, idCols=idCols, mCols=mCols, mColContains=mColContains, FUN=FUN, sep=sep))
}

#' Intended to be used with a wide table that still has ImageChannel and MaskChannel information.
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

#' Example:
#' 
#' # Calculate Geometric Size Ratios
#' ratioData <- calculateChannelLogRatios(x2, comboCol='MaskChannel', valsToPermute=c("1_Lower","1_Upper","3","4","5","6","Cyt","WholeCell"), idCols=c('cId','ImageChannel'), mCols='Geometric.SizeIterable.Total')
#' x2 <- mergeComboData(x2, ratioData, mCol.old='Geometric.SizeIterable.Total', mCol.new='Geometric.SizeIterable.Total.Ratio')
#' rm(ratioData)
#' 
#' # also
#' diffs <- calculateChannelDifferences(x2, comboCol='MaskChannel', idCols=c('cId','ImageChannel'), mCols='Geometric.SizeIterable')
#' x3 <- mergeComboData(x3, diffs, mCol.old='Geometric.SizeIterable', mCol.new='Geometric.SizeIterableDiff')
#' rm(diffs)
#' 
#' The setnames is performed on the combodata just prior to merging with x
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

##### JUNK YARD #####

# getTableList <- function(dir, fileList, class, expt, sampleSize=NULL, cellIds=NULL)
# {
# 	if(!is.null(sampleSize))
# 	{
# 		subSampleSize <- sampleSize / length(fileList)
# 	}
# 	tableList <- list()
# 	for(f in fileList)
# 	{
# 		print(paste0('Reading file: ', file.path(dir, f)))
# 		temp <- fread(file.path(dir, f))
# 		temp$Class <- class
# 		temp$Expt <- expt
# 		if('Z' %in% names(temp))
# 		{
# 			temp[,Z:=NULL]
# 		}
# 		if('A' %in% names(temp))
# 		{
# 			temp[,A:=NULL]
# 		}
# 		if('B' %in% names(temp))
# 		{
# 			temp[,B:=NULL]
# 		}
# 		if(!('ImRow' %in% names(temp)))
# 		{
# 			temp$ImRow <- 1
# 		}
# 		if(!('ImCol' %in% names(temp)))
# 		{
# 			temp$ImCol <- 1
# 		}
# 		if(!('Loc' %in% names(temp)))
# 		{
# 			temp[,Loc:=getLocsFromRCs(ImRow, ImCol, max(ImRow) + 1)]
# 		}
# 		temp$file <- f
# 		temp$cId <- paste(temp$Expt, temp$file, temp$Loc, temp$Id, sep='.')
# 		temp[,ImRow:=NULL]
# 		temp[,ImCol:=NULL]
# 		if(!is.null(cellIds))
# 		{
# 		     rIds <- cellIds
# 		     temp <- temp[cId %in% rIds]
# 		}
# 		if(!is.null(sampleSize))
# 		{
# 			rIds <- trySample(unique(temp$cId), subSampleSize)
# 			temp <- temp[cId %in% rIds]
# 		}
# 		tableList <- append(tableList, list(temp))
# 	}
# 	return(tableList)
# }
# 
# getTableList <- function(dir, featureTable, colocTable, class, expt, sampleSize=NULL, cellIds=NULL)
# {
#      if(!is.null(sampleSize))
#      {
#           subSampleSize <- sampleSize / length(fileList)
#      }
#      tableList <- list()
#      for(f in fileList)
#      {
#           print(paste0('Reading file: ', file.path(dir, f)))
#           temp <- fread(file.path(dir, f))
#           temp$Class <- class
#           temp$Expt <- expt
#           if('Z' %in% names(temp))
#           {
#                temp[,Z:=NULL]
#           }
#           if('A' %in% names(temp))
#           {
#                temp[,A:=NULL]
#           }
#           if('B' %in% names(temp))
#           {
#                temp[,B:=NULL]
#           }
#           if(!('ImRow' %in% names(temp)))
#           {
#                temp$ImRow <- 1
#           }
#           if(!('ImCol' %in% names(temp)))
#           {
#                temp$ImCol <- 1
#           }
#           if(!('Loc' %in% names(temp)))
#           {
#                temp[,Loc:=getLocsFromRCs(ImRow, ImCol, max(ImRow) + 1)]
#           }
#           temp$file <- f
#           temp$cId <- paste(temp$Expt, temp$file, temp$Loc, temp$Id, sep='.')
#           temp[,ImRow:=NULL]
#           temp[,ImCol:=NULL]
#           if(!is.null(cellIds))
#           {
#                rIds <- cellIds
#                temp <- temp[cId %in% rIds]
#           }
#           if(!is.null(sampleSize))
#           {
#                rIds <- trySample(unique(temp$cId), subSampleSize)
#                temp <- temp[cId %in% rIds]
#           }
#           tableList <- append(tableList, list(temp))
#      }
#      return(tableList)
# }

# getXYCSVsAsTableFromDir <- function(dir, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
# {
# 	ret <- list()
# 	fList <- list.files(path = dir, recursive = TRUE)
# 	for(f in fList)
# 	{
# 		if((grepl('x', f) || grepl('y', f)) & grepl('.csv', f))
# 		{
# 			fileName <- strsplit(f, "\\.")[[1]][1]
# 			ret[[fileName]] <- getXYCSVAsTable(dir, f, xName, xExpression, yName, yExpression)
# 		}
# 	}
# 	retTable <- rbindlist(ret)
# 	return(retTable)
# }
# 
# getXYCSVAsTable <- function(dir, file, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
# {
# 	fileName <- strsplit(file, "\\.")[[1]][1]
# 	xy <- strsplit(fileName, "_")[[1]]
# 	y <- as.numeric(substr(xy[1],2,nchar(xy[1])))
# 	x <- as.numeric(substr(xy[2],2,nchar(xy[2])))
# 	xVal <- eval(parse(text=xExpression))
# 	yVal <- eval(parse(text=yExpression))
# 	print(paste0('Reading ', file.path(dir,file), ' as ', xName, '=', xVal, ', ', yName, '=', yVal, '.'))
# 	theTable <- fread(file.path(dir,file))
# 	theTable[,(xName),with=FALSE] <- xVal
# 	theTable[,(yName),with=FALSE] <- yVal
# 	return(theTable)
# }
# 
# getXYArffsAsTableFromDir <- function(dir, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
# {
# 	ret <- list()
# 	fList <- list.files(path = dir, recursive = TRUE)
# 	for(f in fList)
# 	{
# 		if((grepl('x', f) || grepl('y', f)) & grepl('.arff', f))
# 		{
# 			fileName <- strsplit(f, "\\.")[[1]][1]
# 			ret[[fileName]] <- getXYArffAsTable(dir, f, xName, xExpression, yName, yExpression)
# 		}
# 	}
# 	retTable <- rbindlist(ret)
# 	return(retTable)
# }
# 
# getXYArffAsTable <- function(dir, file, xName='SNR', xExpression='(x+1)', yName='BLUR', yExpression='(y+1)*0.05')
# {
# 	fileName <- strsplit(file, "\\.")[[1]][1]
# 	xy <- strsplit(fileName, "_")[[1]]
# 	y <- as.numeric(substr(xy[1],2,nchar(xy[1])))
# 	x <- as.numeric(substr(xy[2],2,nchar(xy[2])))
# 	xVal <- eval(parse(text=xExpression))
# 	yVal <- eval(parse(text=yExpression))
# 	print(paste0('Reading ', file.path(dir,file), ' as ', xName, '=', xVal, ', ', yName, '=', yVal, '.'))
# 	theTable <- read.arff(file.path(dir,file))
# 	theTable[,xName] <- xVal
# 	theTable[,yName] <- yVal
# 	return(data.table(theTable))
# }

#' Depricated
#' StandardizeWideData is preferred
dep_standardizeLongData <- function(x, val.col='Value', by, ignoreMeasurementsContaining=c('ZernikePhase'))
{
	robustScale <- function(x, measurement, ignore)
	{
		if(any(sapply(ignore, grepl, measurement, fixed=T)))
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
	x[,Value:=robustScale(get(val.col),Measurement,ignore=ignoreMeasurementsContaining),by=by]
	return(x)
}

#' Depricated
#' StandardizeWideData is preferred
dep_standardizeLongData2 <- function(x, val.col='Value', by)
{
	robustScale <- function(x)
	{
		m <- median(x, na.rm=TRUE)
		return((x-m)/mad(x, center=m, na.rm=TRUE))
	}
	x[,c(val.col):=robustScale(get(val.col)),by=by]
	return(x)
}

#' Depricated
#' removeNoVarianceCols preferred
dep_removeNoVarianceMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'))
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

#' Depricated
#' removeNoVarianceCols preferred
dep_removeNoMADMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'))
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

#' Depricated
#' readJEXData and readJEXDataTables
dep_getTableListFromDB <- function(db, ds, x, y, objectName, isArff=F, storeFilePath=F, class=NULL, assignClass=T, expt=NULL, repl=NULL, sampleSize=NULL, colsToRemove = c(), cIdCols = c(), fsep='\\\\')
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

#' Depricated
#' use reorganize directly followed by call to sortColsByName(x)
dep_getWideTable <- function(x)
{
	idCols <- getAllColNamesExcept(x, c('Value','Measurement'))
	x <- reorganize(x, idCols)
	x <- sortColsByName(x);
	return(x)
}

#' Depricated
#' Built into readJEXDataTables
dep_trySample <- function(x, n, replace=F, prob=NULL)
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

# Zernike Moment calcualtions

#' JEX outputs an affine invariant Zernike moment (position, scale, & rotation)
#' In order to make it inensity invariant, the magnitude must be divided by the mean.
#' This is meant to be used with long tables.
#' 
#' Not sure yet whether to call this before or after calculating dot products (maybe)
meanNormalizeLongZernikeMoments <- function(x, idCols='cId')
{
	mNames <- getMeasurementNamesContaining(x, 'ZernikeMag')
	for(mName in mNames)
	{
		x <- divideMAbyMBbyRef(x, mName, 'Stats.Mean', idCols=idCols)
	}
	return(x)
}

#' JEX outputs an affine invariant Zernike moment (position, scale, & rotation)
#' In order to make it inensity invariant, the magnitude must be divided by the mean.
#' This is meant to be used with wide tables that HAVEN'T YET INCORPORATED ImageChannel
#' and MaskChannel into the feature name.
#' 
#' Not sure yet whether to call this before or after calculating dot products (maybe)
meanNormalizeWideZernikeMoments <- function(x, meanName='Stats.Mean')
{
	mNames <- getColNamesContaining(x, 'ZernikeMag')
	for(mName in mNames)
	{
		x[, c(mName):=get(mName)/get(meanName)]
	}
	return(x)
}


# Do this AFTER the '_Order_' and '_Rep_' strings have been removed from names
# Do this BEFORE differencing channels
# Do this BEFORE standardizing any data
calculateZernikeDotProduct <- function(x, imageChannelValsToPermute=unique(x$ImageChannel)[!grepl('_',unique(x$ImageChannel,fixed=T))], idCols=c('cId','MaskChannel'))
{
	calculateZernikeReAndIm(x)
	reProducts <- calculateChannelProducts(x, sep='_', comboCol='ImageChannel', valsToPermute=imageChannelValsToPermute, idCols=idCols, mColContains='ZernikeRe')
	imProducts <- calculateChannelProducts(x, sep='_', comboCol='ImageChannel', valsToPermute=imageChannelValsToPermute, idCols=idCols, mColContains='ZernikeIm')
	ret <- addReAndImProducts(reProducts=reProducts, imProducts=imProducts)
	removeColsContaining(x, c('ZernikeRe','ZernikeIm'))
	x <- merge(x, ret, by=c(idCols, 'ImageChannel'), all=T)
	return(x)
}

# This gets called on the entire table to convert Mag and Phase into Re and Im.
calculateZernikeReAndIm <- function(x)
{
	magNames <- getColNamesContaining(x, 'ZernikeMag')
	for(magName in magNames)
	{
		phaseName <- gsub('Mag','Phase',magName)
		if(!(phaseName %in% names(x)))
		{
			next
		}
		reName <- gsub('Mag','Re',magName)
		imName <- gsub('Mag','Im',magName)
		x[, (reName):=get(magName) * cosd(get(phaseName))]
		x[, (imName):=get(magName) * sind(get(phaseName))]
	}
	return(x[])
}

# This is a helper function for calculateZernikeDotProduct.
addReAndImProducts <- function(reProducts, imProducts)
{
	reNames <- getColNamesContaining(reProducts, 'ZernikeRe')
	dotNames <- gsub('Re','Dot',reNames)
	for(reName in reNames)
	{
		imName <- gsub('Re','Im',reName)
		if(!(imName %in% names(x)))
		{
			next
		}
		reProducts[, (reName):=get(reName) + imProducts[[imName]]]
	}
	setnames(reProducts, reNames, dotNames)
	return(reProducts)
	
	# mags <- x[grepl('ZernikeMag',Measurement,fixed=T)]
	# phases <- x[grepl('ZernikePhase',Measurement,fixed=T)]
	# result1 <- copy(mags)
	# result1[,Measurement:=gsub('ZernikeMag','ZernikeRe',Measurement)]
	# result1$Value <- mags$Value * cosd(phases$Value)
	# result2 <- copy(mags)
	# result2[,Measurement:=gsub('ZernikeMag','ZernikeIm',Measurement)]
	# result2$Value <- mags$Value * sind(phases$Value)
	# ret <- rbindlist(list(result1,result2),use.names=TRUE)
	# return(ret)
}

# Do this AFTER the '_Order_' and '_Rep_' strings have been removed from names (part of fixColNames())
getZernikeOrderReps <- function(x)
{
	ret <- unique(tstrsplit(tstrsplit(getColNamesContaining(x, 'ZernikeMag'), 'Mag', keep=2L)[[1]], '_', keep=1L)[[1]])
	return(ret)
}

preprocessData <- function(x,
					  time.col,
					  imageDims,
					  labelDims,
					  doGeometry=F,
					  doHaralick=F,
					  doZernike=F,
					  doDNZandNucwSecOnly=T,
					  doSymmetry=F,
					  channelsToDot=NULL,
					  log.args=NULL, # 			lapply.data.table(x, by=log.args$by, cols=log.args$cols, FUN=logicle, logicle.params=list(), base=getDefault(log.args$base, 10), neg.rm=getDefault(log.args$neg.rm, F), in.place=T)
					  sim.args=NULL, # 			lapply.data.table(x, by=sim.args$by, cols=sim.args$cols, FUN=sim.transform, in.place = T) percentile and 1-percentile are used as upper and lower bounds
					  sim.percentile.args=NULL, # 	lapply.data.table(x, by=sim.percentile.args$by, cols=sim.percentile.args$cols, FUN=sim.transform.percentile, percentile=getDefault(sim.args$percentile, 0.01), in.place = T) percentile and 1-percentile are used as upper and lower bounds
					  logit.args=NULL, #			lapply.data.table(x, by=logit.args$by, cols=logit.args$cols, FUN=logit.transform, base=getDefault(logit.args$base, 10), in.place = T)
					  save.dir=NULL)

{
	library(data.table)
	library(bit64)
	library(MASS)
	library(zoo)
	
	if(!is.null(channelsToDot) && any(!(channelsToDot %in% unique(x$ImageChannel))))
	{
		stop(paste('Some of the channels for calculating Zernike dot products do not exist in the table (', paste(channelsToDot[!(log.args$cols %in% unique(x$ImageChannel))], collapse=','), '). Aborting.'))
	}
	
	x2 <- copy(x)
	setnames(x2, time.col, 'Time')
	
	setorderv(x2, c('e.x', 'e.y', 'cId', 'Time'))
	
	# Avoid standardizing the Id, ImRow, and ImCol etc.
	# lapply.data.table(x, FUN=as.character, in.place=T, cols=c('Id',imageDims,labelDims))
	
	# Get rid of extraneous text in the measurement names that are carried over from JEX/Java
	fixLongTableStringsInCol(x2, 'Measurement')
	
	if(!is.null(log.args) && (is.null(log.args$cols) || any(!(log.args$cols %in% names(x2)))))
	{
		print(uniqueo(x2$Measurement))
		stop(paste('Some of the cols to log transform do not exist in the table (', paste(log.args$cols[!(log.args$cols %in% unique(x2$Measurement))], collapse=','), '). Aborting.'))
	}
	
	if(!is.null(sim.args) && (is.null(sim.args$cols) || any(!(sim.args$cols %in% names(x2)))))
	{
		print(uniqueo(x2$Measurement))
		stop(paste('Some of the cols to log transform do not exist in the table (', paste(sim.args$cols[!(sim.args$cols %in% unique(x2$Measurement))], collapse=','), '). Aborting.'))
	}
	
	if(!is.null(logit.args) && (is.null(logit.args$cols) || any(!(logit.args$cols %in% names(x2)))))
	{
		print(uniqueo(x2$Measurement))
		stop(paste('Some of the cols to log transform do not exist in the table (', paste(logit.args$cols[!(logit.args$cols %in% unique(x2$Measurement))], collapse=','), '). Aborting.'))
	}
	
	##### CONVERT TO WIDE TABLE WITH IMAGEHCANNEL AND MASKCHANNEL #####
	
	# dups <- duplicated(x, by=getAllColNamesExcept(x, c('Value','Measurement')))
	# duh <- x[dups]
	# View(duh)
	
	# First reorganize, keeping ImageChannel and MaskChannel as columns for various calculations
	x2 <- reorganize(x2, idCols=c('cId','ImageChannel','MaskChannel','Time',labelDims), measurementCols=c('Measurement'), valueCols = c('Value'), sep='.')
	x2 <- removeCols(x2,c('Geometric.X','Geometric.Y')) # This data set was quantified with a previous version of FeatureExtraction that doesn't allow 'tracking' X and Y. The position is relative.
	
	# # Normalize standard deviation (PER CELL) and remove redundant metric of variance
	# x2[, ':='(Stats.StdDev=Stats.StdDev/abs(Stats.Mean), Stats.Variance=NULL)]
	
	# Summarize the geometric data that currently has data for each subregion of each mask.
	if(doGeometry)
	{
		x2 <- summarizeGeometry(x2, cellIdCols=c('cId','Time',labelDims), removeXY=F) # Expect a warning from duplication of values (these duplicates are then removed)
	}
	
	# Calculate final rotationally 'invariant' Haralick measures
	if(doHaralick)
	{
		x2 <- calculateRMSofHaralick(x2, removeOriginalHaralickMeasures=T)
	}
	
	# Normalize the Zernike mangnitudes to make them intensity invariant (Should be ok if the mean is neg here atlhough 0 is tricky :-/)
	if(doZernike)
	{
		x2 <- meanNormalizeWideZernikeMoments(x2, meanName='Stats.Mean')
	}
	
	# Summarize the symmetry data that currently has raw symmetry measurements
	if(doSymmetry)
	{
		summarizeSymmetryData(x2, sim.trans=T, logit.trans=T)
	}
	
	# Calculate Zernike Dot Products as a way to measure co-localization
	# Just keep DNZ Zernike features and then also NUCwSEC Zernike features
	if(doZernike && doDNZandNucwSecOnly)
	{
		zernikeCols <- getColNamesContaining(x2, 'Zernike')
		zernikeCols <- zernikeCols[!grepl('NUCwSEC',zernikeCols) & grepl('_',zernikeCols)]
		removeCols(x2, zernikeCols)
	}
	if(!is.null(channelsToDot))
	{
		x2 <- calculateZernikeDotProduct(x2, imageChannelValsToPermute=channelsToDot, idCols=c('cId','MaskChannel',labelDims))
	}
	
	# Calculate Geometric Size Ratios
	# ratioData <- calculateChannelLogRatios(x2, comboCol='MaskChannel', valsToPermute=c("1_Lower","1_Upper","3","4","5","6","Cyt","WholeCell"), idCols=c('cId','ImageChannel'), mCols='Geometric.SizeIterable.Total')
	# x2 <- mergeComboData(x2, ratioData, mCol.old='Geometric.SizeIterable.Total', mCol.new='Geometric.SizeIterable.Total.Ratio')
	# rm(ratioData)
	
	# ##### CLEAN UP WEIRD DATA #####
	# # Remove bad features or features that are no longer necessary
	# # removeColsContaining(x2, c('ImageMoments.Moment','ImageMoments.HuMoment','ImageMoments.NormalizedCentralMoment'))#, getColNamesContaining(x, 'ImageMoments.CentralMoment'))
	# # removeColsContaining(x2, 'HistogramBin')
	# removeColsContaining(x2, c('Circle','Phase'))
	# # filterLBPCodes(x2, nSigma=-2)
	
	##### TRANFORMATIONS #####
	# # Perform and transformations of the data if needed such as log transformations of intensity data.
	# colsToTransform <- c('Stats.Sum', 'Stats.SumOfInverses', 'Stats.SumOfLogs', 'Stats.SumOfSquares', 'Stats.Max', 'Stats.Mean', 'Stats.Min', 'Stats.Median')
	if(!is.null(log.args))
	{
		lapply.data.table(x2, by=log.args$by, cols=log.args$cols, FUN=logicle, logicle.params=list(), base=getDefault(log.args$base, 10), neg.rm=getDefault(log.args$neg.rm, F), in.place=T)
	}
	if(!is.null(sim.args))
	{
		lapply.data.table(x2, by=sim.args$by, cols=sim.args$cols, FUN=sim.transform, in.place = T)
	}
	if(!is.null(sim.percentile.args))
	{
		lapply.data.table(x2, by=sim.percentile.args$by, cols=sim.percentile.args$cols, FUN=sim.transform.percentile, percentile=getDefault(sim.args$percentile, 0.01), in.place = T)
	}
	if(!is.null(logit.args))
	{
		lapply.data.table(x2, by=logit.args$by, cols=logit.args$cols, FUN=logit.transform, base=getDefault(logit.args$base, 10), in.place = T)
	}
	
	fwrite(file=file.path(save.dir, 'x2.csv'), x2, row.names=F)
	x3 <- copy(x2)
	
	##### TRANSPOSE #####
	# Transpose data so that each row contains all data for an individual cell
	x3 <- refactor(x3)
	x3 <- reorganize(x3, idCols=c('cId','Time',labelDims), measurementCols=c('ImageChannel','MaskChannel'), valueCols=getAllColNamesExcept(x3, c('cId','Time','ImageChannel','MaskChannel',labelDims)), sep='.')
	
	# Remove all columns that are ALL NON-FINITE
	removeColsMatching(x3, col.test=function(a){all(!is.finite(a))})
	###########################
	
	fwrite(file=file.path(save.dir, 'x3.csv'), x3, row.names=F)
	return(x3)
}

calculateDrugSensitivityMetrics <- function(x,
								    imageDims,
								    labelDims,
								    maskChan,
								    nucChan,
								    phaseChan,
								    gap.max = 2,
								    initialCellsOnly = T,
								    P.string = '.P',
								    Neg.string = '.Neg',
								    oldTxs=NULL,
								    newTxs=NULL)
{
	library(data.table)
	library(bit64)
	library(MASS)
	library(zoo)
	
	if(length(oldTxs) != length(newTxs))
	{
		stop('Old and new treatement names need to be vectors of the same length. Aborting.')
	}
	
	if(length(oldTxs > 0))
	{
		if(!all(oldTxs %in% unique(x$Tx)))
		{
			warning(paste0('Some of the old Tx names do not exist in the table (', paste(oldTxs[!(oldTxs %in% unique(x$Tx))], collapse=',')))
		}
	}
	
	x3 <- copy(x)
	
	if(length(oldTxs) > 0)
	{
		for(i in 1:length(oldTxs))
		{
			if(oldTxs[i] %in% unique(x3$Tx))
			{
				x3[Tx==oldTxs[i], Tx:=newTxs[i]]
			}
		}
	}
	
	# Make Time real
	if(min(x3$Time) == 0)
	{
		x3[, Time:=Time+1]
	}
	x3[, Time:= 0.5 + (Time-1)*0.5]
	max.time <- max(x3$Time)
	breaks <- merge.vectors(seq(0,max.time,2), seq(0,max.time+1,2))
	x3[, Period.1:=as.numeric(as.character(cut(Time, breaks=breaks, labels=breaks[2:length(breaks)], right=T)))]
	x3[, Period.2:=as.numeric(cut(Time, breaks=c(0,8,16,24,48,96,120,max(Time)), right=T))]
	if(max(x3$Time) <= 120)
	{
		x3[, Period.2:=factor(Period.2, levels=1:6, labels=c('0-8','8-16','16-24','24-48','48-96','96-120'))]
		
	}
	else
	{
		x3[, Period.2:=factor(Period.2, levels=1:7, labels=c('0-8','8-16','16-24','24-48','48-96','96-120',paste0('120-',max(Time))))]
	}
	
	# Filter out cells with time gaps
	x3[, maxGap:=max(getDeltas(sort(Time))), by=c('cId')]
	x3 <- x3[maxGap <= gap.max]
	
	if(initialCellsOnly)
	{
		x3[, minTime:=min(Time), by=c('cId')]
		x3 <- x3[minTime == min(unique(Time))]
	}
	
	# Other Calcs
	# setnames(x3, old=getAllColNamesExcept(x3, c('cId','Time','Period.2','Period.1',labelDims)), new=paste0('Stats.Mean.',getAllColNamesExcept(x3,c('cId','Time','Period.2','Period.1',labelDims))))
	setorderv(x3, c(labelDims,'cId','Time'))
	suppressWarnings(x3[, c('Nuc','Movement','Death','Phase','Death2'):=NULL])
	suppressWarnings(x3[, Nuc:=log10(get(paste0('Stats.Mean.', nucChan, '.', maskChan)))])
	suppressWarnings(x3[, Phase:=log10(get(paste0('Stats.Mean.', phaseChan, '.', maskChan)))])
	suppressWarnings(x3[, Phase:=log10(get(paste0('Stats.Mean.', phaseChan, '.', maskChan)))])
	# x3[, Phase2:=log10(roll.mean(Stats.Mean.1.1, na.rm=T, win.width=3)), by='cId']
	
	setorderv(x3, c(labelDims, 'cId', 'Time'))
	
	# Create Plasma, Neg, and Drug columns if necessary
	x3[, P:=grepl('.P',Tx,fixed=T)]
	x3[, Neg:=grepl('.Neg',Tx,fixed=T)]
	
	return(x3)
}

getEventStatus <- function(time, LD, n.true=3, suffix='')
{
	if(!is.logical(LD))
	{
		stop('This function assumes LD is a logical with 0 = live and 1=dead')
	}
	if(length(time) != length(LD))
	{
		stop('This function assumes time and LD are same length.')
	}
	LD <- rollapply(LD, width=n.true, FUN=sum, na.rm=T, partial=T, align='left')
	event <- which(LD==n.true)
	if(length(event)==0)
	{
		# The event never happened
		ret <- list(LD.time=time[length(time)], LD.status=0)
		names(ret) <- paste0(names(ret), suffix)
		return(ret)
	}
	ret <- list(LD.time=time[event[1]], LD.status=1)
	names(ret) <- paste0(names(ret), suffix)
	return(ret)
}
# LIVE / DEAD
getEventRecoveryStatus <- function(time, LD, n.true=3, suffix='')
{
	if(!is.logical(LD))
	{
		stop('This function assumes LD is a logical with 0 = live and 1=dead')
	}
	if(length(time) != length(LD))
	{
		stop('This function assumes time and LD are same length.')
	}
	LD <- rollapply(LD, width=n.true, FUN=sum, na.rm=T, partial=T, align='left')
	DL <- rollapply(!LD, width=n.true, FUN=sum, na.rm=T, partial=T, align='left')
	deathEvent <- which(LD==n.true)
	lifeEvent <- which(DL==n.true)
	if(length(deathEvent)==0 | length(lifeEvent)==0)
	{
		# The event never happened
		ret <- list(LD.time=time[length(time)], LD.status=0)
		names(ret) <- paste0(names(ret), suffix)
		return(ret)
	}
	else
	{
		recoveryEvent <- lifeEvent[lifeEvent > deathEvent[1]]
		if(length(recoveryEvent)==0)
		{
			# The event never happened
			ret <- list(LD.time=time[length(time)], LD.status=0)
			names(ret) <- paste0(names(ret), suffix)
			return(ret)
		}
	}
	ret <- list(LD.time=time[recoveryEvent[1]], LD.status=1)
	names(ret) <- paste0(names(ret), suffix)
	return(ret)
}

normalizePhase <- function(x)
{
	peaks <- x[Period.1 %in% uniqueo(Period.1)[1:2], getDensityPeaks(copy(Phase), peak.min.sd=10, neighlim=10, n=-1, make.plot=T, plot.args=list(ylim=c(0,5), xlim=getPercentileValues(copy(Phase), c(0.005, 0.995)), main=copy(.BY[[1]]))), by='Tx']
	data.table.plot.all(x3[Period.1 %in% uniqueo(Period.1)[1:2]], percentile.limits=c(0.005, 0.995, 0, 1), xcol='Phase', type='d', by='Tx', density.args=list(draw.area=F, lwd=2), legend.args=list(lty=2), main='Phase Normalization (Before)')
	x[, Phase.norm:=10*(Phase/peaks[peaks$Tx==.BY[[1]]]$peak.x-1), by=c('Tx')]
	data.table.plot.all(x[Period.1 %in% uniqueo(Period.1)[1:2]], percentile.limits=c(0.005, 0.995, 0, 1), xcol='Phase.norm', type='d', by='Tx', density.args=list(draw.area=F), alpha=1, main='Phase Normalization (After)')
}

# Standardize Nuclear Signal at each timepoint to account for drift in labeling intensity.
normalizeNuc <- function(x,
					nuc.lo,
					nuc.hi,
					to.plot=getDefault(uniqueo(x$Tx)[grepl('Veh', uniqueo(x$Tx), fixed=T)][1], uniqueo(x$Tx)[1], test=function(blah){is.null(blah) || !is.finite(blah)}),
					times=NULL,
					sample.size=5000)
{
	data.table.plot.all(x[Tx %in% to.plot], sample.size=sample.size, xcol='Nuc', percentile.limits=c(0.005,0.995), cumulative=F, type='d', by='Period.2', alpha=0.2, density.args=list(draw.area=F), legend.args=list(lty=2), main='Nuc Standardization (Before)')
	standardizeWideData(x, suffix='norm', data.cols=c('Nuc'), na.rm.no.variance.cols=T, col.use.percentiles=T, percentiles=c(nuc.lo,nuc.hi), by=c('Time'))
	if(!is.null(times))
	{
		data.table.plot.all(x[Time %in% times & Tx %in% to.plot], sample.size=5000, percentile.limits=c(0.005,0.995), cumulative=F, xcol='Nuc.norm', type='d', by='Time', alpha=1, density.args=list(draw.area=F), main='Nuc Standardization (After)')
	}
	else
	{
		data.table.plot.all(x[Tx %in% to.plot], sample.size=5000, percentile.limits=c(0.005,0.995), cumulative=F, xcol='Nuc.norm', type='d', by='Period.2', alpha=1, density.args=list(draw.area=F), main='Nuc Standardization (After)')
	}
}

makeDrugSensitivityHistograms <- function(x, PhaseThresh, NucThresh, DeathThresh=0, makePhase=T, makeNuc=T, makeDeath=T, save.dir, type='cmd')
{
	dir.create(file.path(save.dir, 'MovieFrames'), showWarnings=F, recursive = T)
	dir.create(file.path(save.dir, 'Movies'), showWarnings=F, recursive = T)
	x[, Tx:=factor(Tx, levels=uniqueo(Tx))]
	setkey(x, Tx)
	.use.lightFont()
	if(makePhase)
	{
		data.table.plot.all(x, xcol='Phase.norm', by=c('Tx'), type='d', save.plot=T, save.file=file.path(save.dir, 'MovieFrames/PhaseDist_'), v=PhaseThresh, plot.by='Period.1', xlim=c(-3,1), density.args=list(draw.area=F, lwd=2))#, xlim=c(-2.5,2.5))
		makeMovie(full.dir.path=save.dir, in.filename='MovieFrames/PhaseDist_%d.png', out.filename = 'Movies/Phase Histograms.mp4', frame.rate = 2, type=type)
	}
	if(makeNuc)
	{
		data.table.plot.all(x, xcol='Nuc.norm', by=c('Tx'), type='d', save.plot=T, save.file=file.path(save.dir, 'MovieFrames/NucDist_'), v=NucThresh, plot.by=c('Period.1'), xlim=c(-3,3), density.args=list(draw.area=F, lwd=2))
		makeMovie(full.dir.path=save.dir, in.filename='MovieFrames/NucDist_%d.png', out.filename = 'Movies/Nuc Histograms.mp4', frame.rate = 2, type=type)	
	}
	if(makeDeath)
	{
		data.table.plot.all(x, xcol='Death', by=c('Tx'), type='d', save.plot=T, save.file=file.path(save.dir,'MovieFrames/DeathDist_'), v=DeathThresh, plot.by=c('Period.1'), xlim=c(-1,1.3), density.args=list(draw.area=F, lwd=2))
		makeMovie(full.dir.path=save.dir, in.filename='MovieFrames/DeathDist_%d.png', out.filename = 'Movies/Death Histograms.mp4', frame.rate = 2, type=type)
	}
}

plotSurvivalCurve <- function(x.surv, x, Txs=uniqueo(x$Tx), ylab, flip=F, save.plot=T, save.file='Survival Plot.png', ylim=c(0,1), viability.y=0.5, pval.y=0.2, xlim=c(0,50), save.dir, width=4, height=4, Tx.Labels.New=NULL, Tx.Colors=NULL)
{
	
	pval.coord = c(0, ylim[2]*pval.y)
	
	temp <- x.surv[!(LD.time==min(LD.time) & LD.status==1)]
	temp <- temp[Tx %in% Txs]
	inits <- x[Time == uniqueo(Time)[2], list(L=sum(Nuc.norm < NucThresh, na.rm=T), N=.N, init.viability=round(100*sum(Phase.norm >= PhaseThresh, na.rm=T)/.N, 0)), by=c('Tx','Time')]
	paste.cols(inits, cols=c('Tx','init.viability'), name='txt', sep=':')
	inits[, txt2:=paste(txt,'%', sep='')]
	if(!is.null(Tx.Labels.New))
	{
		if(length(Tx.Labels.New) > 0 && (length(Tx.Labels.New) == length(Txs)))
		{
			temp[, Tx:=factor(Tx, levels=Txs, labels=Tx.Labels.New)]
		}
		else
		{
			stop('New and old Tx Label vectors should be the same length as the Txs vector. Aborting.')
		}
	}
	else
	{
		temp[, Tx:=factor(Tx, levels=Txs)]
	}
	
	setkey(temp, Tx)
	fit <- survfit(Surv(LD.time, LD.status) ~ Tx, data = temp)
	labs <- getUniqueCombos(temp, c('Tx'))
	makeComplexId(labs, c('Tx'))
	if(!is.null(Tx.Colors))
	{
		labs$cols <- Tx.Colors
	}
	else
	{
		labs$cols <- (loopingPastels(1:length(labs$cId), a=1, max.k=length(labs$cId)))
	}
	
	.use.lightFont()
	
	if(flip)
	{
		daPlot <- ggsurvplot(fit, fun='event', pval.coord=pval.coord, data = temp, censor=F, pval=T, palette=labs$cols, xlim=xlim, ylim=ylim, conf.int = T, ylab=ylab, xlab='Time [h]', legend.title='', legend.labs=labs$cId, ggtheme = theme_classic2(base_family = "Open Sans Light", base_size = 14))
	}
	else
	{
		daPlot <- ggsurvplot(fit, pval.coord=pval.coord, data = temp, censor=F, pval=T, palette=labs$cols, xlim=xlim, ylim=ylim, conf.int = T, ylab=ylab, xlab='Time [h]', legend.title='', legend.labs=labs$cId, ggtheme = theme_classic2(base_family = "Open Sans Light", base_size = 14))
	}
	
	daPlot$plot <- daPlot$plot+ 
		ggplot2::annotate("text", 
					   x = min(x$Time), y = ylim[1]+(ylim[2]-ylim[1])*viability.y, # x and y coordinates of the text
					   label = paste0("Init. Viability:\n", paste(inits[Tx %in% Txs]$txt2, collapse="\n")),
					   size = 2.5,
					   hjust=0)+
		guides(colour = guide_legend(nrow = (length(Txs) %/% 4) + 1))
	stats <- as.data.table(pairwise_survdiff(Surv(LD.time, LD.status) ~ Tx, data = temp)$p.value, keep.rownames=T)
	stats.sym <- lapply.data.table(stats, getPSymbol, cols=getAllColNamesExcept(stats, 'rn'), by=c('rn'))
	if(save.plot)
	{
		fwrite(stats, file=file.path(save.dir, paste0(tools::file_path_sans_ext(save.file),' - pvalues.csv')))
		png(file.path(save.dir, save.file), res=300, units='in', width=width, height=height)
		print(daPlot)
		dev.off()
	}
	else
	{
		print(daPlot)
		print(stats)
	}
	return(list(stats=stats, stats.sym=stats.sym, daPlot=daPlot))
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

