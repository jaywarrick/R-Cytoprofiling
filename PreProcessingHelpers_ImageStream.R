library(data.table)
library(foreign)

##### Visualization #####

browseShinyData <- function()
{
	sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/ui.R')
	sourceGitHubFile(user='jaywarrick', repo='R-General', branch='master', file='DataClassBrowser/server.R')
	shinyApp(ui=myUI, server=myServer)
}

plotHist <- function(x, feature)
{
	breaks=c(-1000, seq(-4,4,0.5), 1000)
	wt <- x[Class == 'WT'][[feature]]
	mt <- x[Class == 'MT'][[feature]]
	cmt <- rgb(0,0,1,0.8)
	cwt <- rgb(1,0,0,0.8)
	wtd <- density(wt, from=-4, to=4)
	mtd <- density(mt, from=-4, to=4)
	if(max(wtd$y) > max(mtd$y))
	{
		plot(wtd, col='red', xlim=c(-4,4), main='', xlab=feature)
		lines(mtd, col='blue')
	}
	else
	{
		plot(mtd, col='blue', xlim=c(-4,4), main='', xlab=feature)
		lines(wtd, col='red')
	}
	legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)
}

##### General #####

resample <- function(x, ...)
{
	x[sample.int(length(x), ...)]
}

getLocsFromRCs <- function(r, c, numRows, zeroIndexing=true)
{
	if(zeroIndexing)
	{
		r + max(numRows) * c
	}
	else
	{
		(r-1) + max(numRows) * (c-1)
	}
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

getTableList <- function(dir1, fileList, class, expt, sampleSize=NULL, cIds=NULL)
{
	if(!is.null(sampleSize))
	{
		subSampleSize <- sampleSize / length(fileList)
	}
	tableList <- list()
	for(f in fileList)
	{
		path <- file.path(dir1, f)
		print(paste0('Reading file: ', path))
		temp <- fread(input=path)
		temp$Class <- class
		temp$Expt <- expt
		fileXY <- unlist(strsplit(f, 'x', fixed=T))
		temp$File <- paste0('x', fileXY[length(fileXY)])
		temp$cId <- paste(temp$Expt, temp$File, temp$Id, sep='.')
		if(!is.null(sampleSize))
		{
			rIds <- trySample(unique(temp$cId), subSampleSize)
			temp <- temp[cId %in% rIds]
		}
		else if(!is.null(cIds))
		{
			temp <- temp[cId %in% cIds]
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

applySimilarityTransform <- function(x, measurements='Stats.PearsonsCorrelationCoefficient', bounds=c(-1,1), removeOld=TRUE)
{
	for(measure in measurements)
	{
		x[, c(paste0(measure, '.Dilate')) := log((abs(bounds[1]) + get(measure))/(abs(bounds[2]) - get(measure)))]
		if(removeOld)
		{

			x[, c(measure):=NULL]
		}
	}
}

summarizeGeometry <- function(x, idCols)
{
	x[, MaskChannel2 := tstrsplit(MaskChannel, '.p', fixed=T, keep=1L)]
	x[is.finite(Geometric.SizeIterable), ':='(weights=Geometric.SizeIterable/(sum(Geometric.SizeIterable, na.rm=T)), countWeights=Geometric.SizeIterable/(max(Geometric.SizeIterable, na.rm=T))), by=c(idCols, 'MaskChannel2')]
	x[, ':='(Geometric.FeretsAspectRatio = Geometric.MaximumFeretsDiameter/Geometric.MinimumFeretsDiameter, Geometric.EllipseAspectRatio = Geometric.MajorAxis/Geometric.MinorAxis)]
	x[, ':='(Geometric.MaximumFeretsDiameter=NULL, Geometric.MinimumFeretsDiameter=NULL, Geometric.MajorAxis=NULL, Geometric.MinorAxis=NULL)]

	# Decide how to combine different geometric features
	#geomFeatures <- c('Convexity', 'Solidity', 'SizeIterable', 'BoundarySize', 'MainElongation', 'Circularity',
	#'Boxivity', 'Eccentricity', 'MajorAxis', 'MaximumFeretsDiameter', 'MinimumFeretsDiameter',
	#'MinorAxis', 'Roundness', 'X', 'Y')
	geomFeatures_Total <- c('Geometric.SizeIterable', 'Geometric.BoundarySize')
	geomFeatures_SizeWeightedMean <- c('Geometric.SizeIterable', 'Geometric.Convexity', 'Geometric.Solidity', 'Geometric.MainElongation', 'Geometric.Circularity', 'Geometric.Boxivity', 'Geometric.Eccentricity', 'Geometric.BoundarySize','Geometric.FeretsAspectRatio','Geometric.EllipseAspectRatio','Geometric.Roundness')

	for(feature in geomFeatures_Total)
	{
		x[is.finite(get(feature)), (paste0(feature, '.Total')):=sum(get(feature), na.rm=TRUE), by=c(idCols, 'MaskChannel2')]
	}
	for(feature in geomFeatures_SizeWeightedMean)
	{
		x[is.finite(get(feature)), (feature):=sum(get(feature)*weights, na.rm=TRUE), by=c(idCols, 'MaskChannel2')]
	}
	x[is.finite(countWeights), ':='(N=as.double(.N), weightedN=sum(countWeights)), by=c(idCols, 'MaskChannel2')]

	getPointStats <- function(x, y, weights)
	{
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
		ret <- lapply(list(Intensity=intensity(pp), WeightedIntensity=intensity(pp, weights=weights), ConvexArea=area(hull), ConvexPerimeter=perimeter(hull), Diameter=circ$rad), as.double)
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
	x[is.finite(Geometric.X), c('Geometric.SubRegionIntensity', 'Geometric.SubRegionWeightedIntensity', 'Geometric.SubRegionConvexArea', 'Geometric.SubRegionConvexPerimeter', 'Geometric.SubRegionRadius', 'Geometric.SubRegionCircularity'):=getPointStats(Geometric.X, Geometric.Y, weights), by=c(idCols, 'MaskChannel2')]
	x[, ':='(MaskChannel=MaskChannel2, MaskChannel2=NULL, weights=NULL, countWeights=NULL, Geometric.X=NULL, Geometric.Y=NULL)]
	return(x)
}

labeledBoxPlot <- function(data, labelCols=c('name'), ord=c(1), valueCol='MeanDecreaseAccuracy', col=c('blue'), bottom=10, left=4, top=1, right=1, num.offset=1, axlab.offset=5, ...)
{
	library(data.table)
	myFormula <- reformulate(termlabels=labelCols, response=valueCol)

	setorderv(rfImp, labelCols, ord)
	labels <- unique(data[[labelCols[length(labelCols)]]])
	at <- (1:length(labels))

	print(labels)
	par(mar=c(bottom,left,top,right), mgp=c(axlab.offset,num.offset,0))
	duh <- boxplot(myFormula, data=data, las=2, mar=c(6,6,6,6), col=col, xaxt='n', ...)
	axis(1, labels = labels, at=at, las=2)
}

# getSubRegionSizeWeightedMeans2 <- function(x, weights, features, featureResultNames=features, idCols=getAllColNamesExcept(x, c('MaskChannel','Value')))
# {
# 	subWeights <- function(weights, this.cId, this.MaskChannel2)
# 	{
# 		return(weights[cId==this.cId & MaskChannel2==this.MaskChannel2]$Value)
# 	}
# 	maskChannel2Index <- match('MaskChannel2', idCols)
# 	cIdIndex <- match('cId', idCols)
# 	results <- x[Measurement %in% features, list(MaskChannel=.BY[[maskChannel2Index]], Value=sum(Value * subWeights(weights=weights, this.cId=.BY[[cIdIndex]], this.MaskChannel2=.BY[[maskChannel2Index]]), na.rm=T)), by=idCols]
# 	results[, MaskChannel2:=NULL]
# 	replaceMeasurementNames(results, features, featureResultNames)
# 	return(results)
# }
#
# getSubRegionTotals <- function(x, features, featureResultNames=paste0(features, '.Total'), idCols=getAllColNamesExcept(x, c('MaskChannel','Value')))
# {
# 	# Might need to do na.rm=T for sum(...)
# 	maskChannel2Index <- match('MaskChannel2', idCols)
# 	results <- x[Measurement %in% features, list(MaskChannel=.BY[[maskChannel2Index]], Value=sum(Value, na.rm=T)), by=idCols]
# 	results[, MaskChannel2:=NULL]
# 	replaceMeasurementNames(results, features, featureResultNames)
# 	return(results)
# }
#
# getSubRegionTotals2 <- function(x, features, featureResultNames=paste0(features, '.Total'), idCols=getAllColNamesExcept(x, c('MaskChannel','Value')))
# {
# 	# Might need to do na.rm=T for sum(...)
# 	maskChannel2Index <- match('MaskChannel2', idCols)
# 	results <- x[Measurement %in% features, list(MaskChannel=.BY[[maskChannel2Index]], Value=sum(Value, na.rm=T)), by=idCols]
# 	results[, MaskChannel2:=NULL]
# 	replaceMeasurementNames(results, features, featureResultNames)
# 	return(results)
# }
#
# getSubRegionCountsAndSizeWeightedCounts <- function(x, idCols=getAllColNamesExcept(x, c('MaskChannel','Value')))
# {
# 	# rename MaskChannel - Also rename the MaskChannel to eliminate the '.p#' designation
# 	# weightedN - Use SizeIterable to make this weighted count since we need it for weighting anyway
# 	# N - The number of subregions for this cell, Array.X, Array.Y, and MaskChannel
# 	results <- x[Measurement == 'Geometric.SizeIterable', list(MaskChannel=MaskChannel2, Value=.N), by=idCols]
# 	results$Measurement <- 'N'
# 	results2 <- x[Measurement == 'Geometric.SizeIterable', list(MaskChannel=MaskChannel2, Value=sum(Value/max(Value))), by=idCols]
# 	results2$Measurement <- 'weightedN'
# 	results <- rbindlist(list(results, results2), use.names = T)
# 	results[, MaskChannel2:=NULL]
# 	return(results)
# }
#
# getSubRegionCountsAndSizeWeightedCounts2 <- function(x, idCols=getAllColNamesExcept(x, c('MaskChannel','Value')))
# {
# 	# rename MaskChannel - Also rename the MaskChannel to eliminate the '.p#' designation
# 	# weightedN - Use SizeIterable to make this weighted count since we need it for weighting anyway
# 	# N - The number of subregions for this cell, Array.X, Array.Y, and MaskChannel
# 	results <- x[Measurement == 'Geometric.SizeIterable', list(MaskChannel=MaskChannel2, Value=.N), by=idCols]
# 	results$Measurement <- 'N'
# 	results2 <- x[Measurement == 'Geometric.SizeIterable', list(MaskChannel=MaskChannel2, Value=sum(Value/max(Value))), by=idCols]
# 	results2$Measurement <- 'weightedN'
# 	results <- rbindlist(list(results, results2), use.names = T)
# 	results[, MaskChannel2:=NULL]
# 	return(results)
# }

findDuplicateValuesInWideTableFaultyResult <- function(x)
{
	troubleNames <- names(x5)[sapply(x, function(x){if(is.numeric(x)){return(max(x) > 1)}else{return(TRUE)}})]
	return(x[,troubleNames,with=F])
}

removeExtraneousColumns <- function(x)
{
	dumbCols <- c(getColNamesContaining(x, 'Phase'), getColNamesContaining(x, 'ImageMoments.Moment'), getColNamesContaining(x, 'ImageMoments.HuMoment'), getColNamesContaining(x, 'ImageMoments.NormalizedCentralMoment'))#, getColNamesContaining(x, 'ImageMoments.CentralMoment'))
	dumbCols <- unique(dumbCols)
	print('Removing the following extraneous columns of information...')
	for(colName in dumbCols)
	{
		print(colName)
	}
	x[,(dumbCols):=NULL]
	return(x)
}

divideColAByColB <- function(x, colA, colB)
{
	x[get(colB)==0,(colA):=NA]
	x[get(colB)!=0,(colA):=get(colA)/get(colB)]
	return(x)
}

removeColsWithNonfiniteVals <- function(x)
{
	duh <- x[,lapply(.SD, function(y){length(which(!is.finite(y))) > 0}), .SDcols=getNumericCols(x)]
	duh2 <- getNumericCols(x)[as.logical(as.vector(duh))]
	if(length(duh2 > 0))
	{
		print("Removing cols with infinite values...")
		for(col in duh2)
		{
			print(col)
			x[,(col):=NULL]
		}
	}
	else
	{
		print("No non-finite data found in any column. Yay!")
	}
}

removeCellsWithNonfiniteVals <- function(x, colDeletionThreshold=0.05)
{
	duh <- x[,lapply(.SD, function(y){as.double(length(which(!is.finite(y))))/as.double(length(y)) > colDeletionThreshold}), .SDcols=getNumericCols(x)]
	duh1 <- getNumericCols(x)[as.logical(as.vector(duh))]
	for(col in duh1)
	{
		print(paste0("Removing column ", col, " as more than 5% of cells have NAs/Infs for this feature"))
		x[, c(col):=NULL]
	}
	duh <- x[,lapply(.SD, function(y){length(which(!is.finite(y))) > 0}), .SDcols=getNumericCols(x)]
	duh2 <- getNumericCols(x)[as.logical(as.vector(duh))]
	if(length(duh2 > 0))
	{
		print("Finding cells to remove...")
		print(paste0('Checking column: ', duh2[1]))
		cIds <- x[!is.finite(get(duh2[1]))]$cId
		if(length(cIds) > 0)
		{
			ret <- x[is.finite(get(duh2[1]))]
			cIds <- unique(cIds)
			print('Removing the following ids')
			print(paste(cIds, collapse=', '))

			# print(ret$cId)

			# Call recursively
			return(removeCellsWithNonfiniteVals(ret))
		}
		else
		{
			print("Shouldn't ever get here!!!")
		}
	}
	else
	{
		print("No more non-finite data was found in any column. Yay!")
		return(x)
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

removeColsContainingNames <- function(x, namesToMatch)
{
	colsToRemove <- getColNamesContaining(x, namesToMatch[1])
	print(paste0("Removing colums with names containing..."))
	for(nameToMatch in namesToMatch)
	{
		print(nameToMatch)
		colsToRemove <- colsToRemove[colsToRemove %in% getColNamesContaining(x, nameToMatch)]
	}
	for(colToRemove in unique(colsToRemove))
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
	x <- reorganize(x, idCols=idCols)
	x <- sortColsByName(x);
	return(x)
}

sortColsByName <- function(x)
{
	setcolorder(x, sort(names(x)))
}

# standardizeWideData <- function(x)
# {
# 	removeNoVarianceCols(x)
# 	robustScale <- function(x)
# 	{
# 		m <- median(x, na.rm=TRUE)
# 		return((x-m)/mad(x, center=m, na.rm=TRUE))
# 	}
# 	x[,lapply(.SD, function(x){if(is.numeric(x)){return(robustScale(x))}else{return(x)}})]
# }

# removeNoVarianceCols <- function(x)
# {
# 	namesToRemove <- getNoVarianceCols(x)
# 	if(length(namesToRemove) > 0)
# 	{
# 		print("Removing cols with a variance of zero...")
# 		for(name in namesToRemove)
# 		{
# 			print(name)
# 			x[,(name):=NULL]
# 		}
# 	}
# }

# getNoVarianceCols <- function(x)
# {
# 	tempSD <- function(y){sd(y, na.rm = TRUE)}
# 	tempNames <- x[,lapply(.SD, tempSD), .SDcols=getNumericCols(x)]
# 	return(names(tempNames)[as.numeric(as.vector(tempNames))==0])
# }

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

intIntensityNormalizeCentralMoments <- function(x)
{
	mNames <- getMeasurementNamesContaining(x, 'ImageMoments.CentralMoment')
	for(mName in mNames)
	{
		x <- divideMAbyMBbyRef(x, mName, 'Stats.Sum')
	}
	return(x)
}

meanNormalizeZernikeMoments <- function(x)
{
	mNames <- getMeasurementNamesContaining(x, 'ZernikeMag')
	for(mName in mNames)
	{
		x <- divideMAbyMBbyRef(x, mName, 'Stats.Mean')
	}
	return(x)
}

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
	names <- names[names %in% unique(x$Measurement)]
	for(name in names)
	{
		print(name)
	}
	x <- x[!(Measurement %in% names)]
	return(x)
}

replaceMeasurementNames <- function(x, oldNames, newNames)
{
	print("Replacing the following Measurement names...")
	oldNames.temp <- oldNames[oldNames %in% unique(x$Measurement)]
	for(name in oldNames.temp)
	{
		print(name)
	}

	print("with...")
	newNames.temp <- newNames[oldNames %in% unique(x$Measurement)]
	for(name in newNames.temp)
	{
		print(name)
	}

	x <- x[Measurement %in% oldNames.temp, Measurement:=newNames.temp[match(Measurement,oldNames.temp)]]
	return(x)
}

standardizeLongData <- function(x, by=c('MaskChannel','ImageChannel','Measurement','Expt'), use.mad=TRUE)
{
	robustStandardize <- function(x, measurement, use.mad)
	{
		if(substr(measurement,1,12) == 'ZernikePhase')
		{
			return(x)
		}
		else
		{
			if(use.mad)
			{
				m <- median(x, na.rm=TRUE)
				return((x-m)/mad(x, center=m, na.rm=TRUE))
			}
			else
			{
				m <- mean(x, na.rm=TRUE)
				return((x-m)/sd(x, center=m, na.rm=TRUE))
			}
		}
	}
	x <- removeNoVarCombos(x, by=by, use.mad=use.mad)
	x[,Value:=robustStandardize(Value,Measurement,use.mad=use.mad),by=by]
	return(x)
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

removeNonFiniteVarCombos <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'))
{
	# If any of the experiments saw 0 variance (see how unique call does not include 'Expt'), then we should eliminate that measure (combo of Measurement, ImageChannel, and MaskChannel) from all experiments
	toRemove <- unique(x[!is.finite(Value),c('Measurement','ImageChannel','MaskChannel')])
	toRemove[, combo:=paste(Measurement,ImageChannel,MaskChannel,sep=' ')]

	if(nrow(toRemove)>0)
	{
		print("Removing measurements with non-finite values...")
		maxOption <- getOption('max.print')
		options(max.print=nrow(toRemove) + 10)
		print(toRemove$combo, nrows=nrow(toRemove))
		options(max.print=maxOption)
		y <- x[!(paste(Measurement,ImageChannel,MaskChannel,sep=' ') %in% toRemove$combo)]
		x[, VAR:=NULL]
		y[, VAR:=NULL]
		return(y)
	}else
	{
		return(x)
	}
}

removeNoVarCombos <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'), use.mad=TRUE)
{
	# For each experiment calculate the variance/mad seen
	if(use.mad)
	{
		x[,VAR:=mad(get(val), na.rm=TRUE), by=by]
	}
	else
	{
		x[,VAR:=var(get(val), na.rm=TRUE), by=by]
	}

	# If any of the experiments saw 0 variance (see how unique call does not include 'Expt'), then we should eliminate that measure (combo of Measurement, ImageChannel, and MaskChannel) from all experiments
	toRemove <- unique(x[VAR == 0,c('Measurement','ImageChannel','MaskChannel')])
	toRemove[, combo:=paste(Measurement,MaskChannel,ImageChannel,sep=' ')]

	if(nrow(toRemove)>0)
	{
		print("Removing measurements with 0 VAR/MAD...")
		maxOption <- getOption('max.print')
		options(max.print=nrow(toRemove) + 10)
		print(toRemove$combo, nrows=nrow(toRemove))
		options(max.print=maxOption)
		y <- x[!(paste(Measurement,MaskChannel,ImageChannel,sep=' ') %in% toRemove$combo)]
		x[, VAR:=NULL]
		y[, VAR:=NULL]
		return(y)
	}else
	{
		x[, VAR:=NULL]
		return(x)
	}
}

# removeNoMADMeasurements <- function(x, val='Value', by=c('MaskChannel','ImageChannel','Measurement','Expt'))
# {
# 	# Tempororarily add a column in the table with stdev in it
# 	x[,MAD:=mad(get(val), na.rm=TRUE), by=by]
# 	toRemove <- unique(x[MAD == 0]$Measurement)
# 	if(length(toRemove)>0)
# 	{
# 		print("Removing measurements with 0 MAD...")
# 		for(m in toRemove)
# 		{
# 			print(m)
# 		}
# 		y <- x[!(Measurement %in% toRemove)]
# 		x[, MAD:=NULL]
# 		y[, MAD:=NULL]
# 		return(y)
# 	}else
# 	{
# 		x[, MAD:=NULL]
# 		return(x)
# 	}
# }

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
	replaceSubStringInAllRowsOfCol(x,' ','',col)
	replaceSubStringInAllRowsOfCol(x,':','_',col)
}

##### Feature Calculations #####

unmergeChannelNames <- function(channelString)
{
	temp <- unlist(strsplit(channelString,'_minus_',fixed=TRUE))
	return(list(channel1=temp[1], channel2=temp[2]))
}

calculateChannelDifferences <- function(x)
{
	if(length(unique(x$ImageChannel)) > 1)
	{
		# Calculate differences between channels for each Cell and Measurement (but keep other column information too so include other cols in 'by')
		idCols <- getAllColNamesExcept(x, c('Value','ImageChannel'))
		return(x[ImageChannel != 'None' & !grepl('_',ImageChannel,fixed=T),list(ImageChannel=getComboNames(ImageChannel), Value=getComboDifferences(Value)), by=idCols])

	}else
	{
		# return an empty table with the same columns as provided
		return(x[FALSE])
	}
}

# Meant to be called on a subset of the main table
calculateChannelProducts <- function(x)
{
	if(length(unique(x$ImageChannel)) > 1)
	{
		# Calculate differences between channels for each Cell and Measurement (but keep other column information too so include other cols in 'by')
		idCols <- getAllColNamesExcept(x, c('Value','ImageChannel'))
		x2 <- x[ImageChannel != 'None',list(ImageChannel=getComboNames(ImageChannel, '_times_'), Value=getComboProducts(Value)), by=idCols]
	}else
	{
		# return an empty table with the same columns as provided
		return(x[FALSE])
	}
}

getComboNames <- function(x, operation='_minus_')
{
	if(length(x) < 2)
	{
		return(NULL)
	}
	temp <- combn(x, 2)
	#print(temp)
	temp <- paste0(temp[1,],operation,temp[2,])
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

