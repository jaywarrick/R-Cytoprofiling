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

getTableList <- function(dir, fileList, isArff=F, storeFilePath=F, class=NULL, assignClass=T, expt=NULL, repl=NULL, sampleSize=NULL, colsToRemove = c(), cIdCols = c())
{
	if(!is.null(sampleSize))
	{
		subSampleSize <- sampleSize / length(fileList)
	}
	tableList <- list()
	
	# For each file in the fileList
	for(f in fileList)
	{
		# Read the file in
		print(paste0('Reading file: ', file.path(dir, f)))
		
		if(isArff)
		{
			library(foreign)
			temp <- data.table(read.arff(file.path(dir, f)))
		}
		else
		{
			temp <- fread(file.path(dir, f))
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
		if(length(idColsFound) != length(cIdCols))
		{
			warning(cat('The specified cIdCols (', cIdCols[!(cIdCols %in% names(temp))], 'is/are not column names of the table being retrieved... (', names(temp), ')'))
		}
		temp[,c('cId'):=paste(mapply(function(x){unique(as.character(x))}, mget(idColsFound)), collapse='.'), by=idColsFound]
		
		print(temp[cId == '118.11.1.HS5'])
		
		# put the complex Id first and the class column last
		setcolorder(temp, c('cId', names(temp)[names(temp) != 'cId']))
		
		# Put the 'Class' column as the last column of the table
		setcolorder(temp, c(names(temp)[names(temp) != 'Class'], 'Class'))
		
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

removeColsContainingAny <- function(x, colNames)
{
	dumbCols <- c()
	for(dumbCol in colNames)
	{
		dumbCols <- c(dumbCols, getColNamesContaining(x, dumbCol))
	}
	dumbCols <- unique(dumbCols)
	print('Removing the following extraneous columns of information...')
	for(dumbCol in dumbCols)
	{
		print(dumbCol)
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

removeColsContaining <- function(x, name)
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
	replaceStringInColNames(x,'_Order_','')
	replaceStringInColNames(x,'_Rep_','')
	replaceStringInColNames(x,'$','.')
	replaceStringInColNames(x,'net.imagej.ops.Ops.','')
	replaceStringInColNames(x,' ','')
	replaceStringInColNames(x,':','_')
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

removeIncompleteRows <- function(x)
{
	valid <- NULL
	for(colName in names(x))
	{
		#print(colName)
		if(is.numeric(x[, colName, with=F][[1]][1]))
		{
			temp <- is.finite(x[,colName,with=F][[1]])
			if(is.null(valid))
			{
				valid <- temp
			}
			else
			{
				valid <- valid & temp
			}
		}
	}
	cat('Removing rows... ', which(!valid), sep=',')
	return(x[valid])
}

calculateLogRatiosOfColsContainingName <- function(x, name)
{
	mNames <- getColNamesContaining(x, name)
	combos <- combn(mNames,2)
	for(j in seq_along(combos[1,]))
	{
		combo <- combos[,j]
		ending1 <- substring(combo[1], first=nchar(name) + 2)
		ending2 <- substring(combo[2], first=nchar(name) + 2)
		x[,c(paste0(name, "LR.", ending1, ".", ending2)) := log(get(combo[1]) / get(combo[2]))]
	}
	return(x)
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

integratedIntensityNormalizeCentralMoments <- function(x)
{
	# mNames <- getMeasurementNamesContaining(x, 'ImageMoments.CentralMoment')
	# for(mName in mNames)
	# {
	# 	x <- divideMAbyMBbyRef(x, mName, 'Stats.Sum')
	# }
	# return(x)
	mNames <- getColNamesContaining(x, 'ImageMoments.CentralMoment')
	newMNames <- paste(mNames, '.M00Normalized', sep='')
	for(mName in mNames)
	{
		x[,c(mName) := get(mName)/Stats.Sum]
	}
	setnames(x, mNames, newMNames)
	return(x)
}

meanNormalizeZernikeMoments <- function(x)
{
	# mNames <- getMeasurementNamesContaining(x, 'ZernikeMag')
	# for(mName in mNames)
	# {
	# 	x <- divideMAbyMBbyRef(x, mName, 'Stats.Mean')
	# }
	# return(x)
	mNames <- getColNamesContaining(x, 'ZernikeMag')
	for(mName in mNames)
	{
		x[,c(mName) := get(mName)/Stats.Mean]
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
		return(x[ImageChannel != 'None' & !grepl('_dot_',ImageChannel,fixed=T),list(ImageChannel=getComboNames(ImageChannel), Value=getComboDifferences(Value)), by=idCols])

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

