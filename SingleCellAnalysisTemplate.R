rm(list=ls())
gc()
# If needed, run the following script to install necessary packages before running the rest of the script.
# source(curl_download(url='https://raw.githubusercontent.com/jaywarrick/R-General/master/InstallPackages.R', tempfile()))
library(curl)
library(data.table)
library(bit64)
source(curl_download(url='https://raw.githubusercontent.com/jaywarrick/R-General/master/.Rprofile', tempfile(fileext = '.R')))
source(curl_download(url='https://raw.githubusercontent.com/jaywarrick/R-Cytoprofiling/master/PreProcessingHelpers.R', tempfile(fileext = '.R')))
library(MASS)
library(zoo)
library(survminer)
library(survival)

#### INSTRUCTIONS ####
# Use Ctrl+Enter (or CMD+Enter on a mac) to step through each lines of code as necessary
# Most of the code always stays the same but some parts need to be altered to match the imaging parameters (e.g., channel names etc)
# There are 6 points in the code that should be revisted for each dataset to ensure proper numerical analysis (STEP A, B, C, D, E, and F)
# Simply run through the code and read the information associated with each of the 6 steps to perform the analysis.

#### STEP A ####
# Alter these parameters to match your imaging setup, then run the code up to the "STEP B" header
save.dir <- '<your folder path>' # !!!CHANGE THIS!!! to match the RStudio project directory
maskChan <- 'WholeCell' # Mask Channel. Shouldn't change really, but needs to match the channel name from JEX
nucChan <- 3 # Nuc Fluor. Channel (shouldn't really change), but needs to match channel name from JEX
phaseChan <- 2 # Phase Channel (shouldn' really change)
gap.max <- 2 # How many timepoints in a row is a cell allowed to be lost during tracking 
initialCellsOnly <- F # Should we only analyze cells that were present at t=0?
P.string <- '.Pos' # Can leave, only applies to studies with and without Neg fraction
Neg.string <- '.Neg' # Can leave, only applies to studies with and without Neg fraction
time.completeness <- 0.1 # Filter out cells that cannot be tracked for longer than this fraction of the entire timelapse (e.g., 0.1 = 10% of timelapse)
sample.size <- -1 # Leave as -1 to sample ALL cells. Set as a positive number to just sample that many cells randomly for each well
lines.with <- c('Id','Stats$Mean','Stats$Sum','Geometric$COMX','Geometric$COMX','Geometric$SizeIterable','Geometric$Circularity') # Usually don't change but can add to list if you want. In the big csv master table, we are only going to grab lines that have these strings in them. This dramatically reduces how much memory we use.
oldTxs='Dara0' # Usually don't change. Use this in conjuction with the 'newTxs' variable to do a 'Find and replace' while you read the data in. It only applies to 'Tx' labels (see the 'labels' argument in the 'readJexData' call below)
newTxs='Veh' # Usually don't change
####

#### STEP B ####
# Create a data structure that represents the JEX database file/table objects that should be imported. This only loads the "metadata" of the table, but not the table itself
# Simply replace the JEX database path and object information to point R to the appropriate database objects in JEX
# dbPath is the path to the folder with the JEX database.
# ds is the name of the "Dataset" within the JEX database where you want to get a database file object from
# e.x and e.y are the X and Y location in the array of wells (index starting at 0) within the dataset where you want to get a database file object
# Type is always 'File' for this type of analysis although this can be used to get the metadata associated with any type of JEX database object
# Name is the name of the database object to get (i.e., the csv file of raw image calculations from JEX)
# Labels is set to whatever "name=<value>" pairs (i.e., metadata) that you would like to assign to these file objects. This makes it easier for plotting later (e.g., Tx=Drug vs Tx=Control)
# Valid can be used to set the information from a particular well as valid or invalid. Generally TRUE, however, if set to false, the data won't be imported in the next steps (i.e., will be ignored)
# The index within jData doesn't matter, but just needs to be unique (i.e., jData[['0']] could be jData[['Yay']] but each readJEXData should be a separate item in the jData list)
jData <- list()
jData[['0']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=0, e.y=0, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='WT.MM',Valid='true'))
jData[['1']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=0, e.y=1, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='WT.MM.bort',Valid='true'))
jData[['2']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=0, e.y=2, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='shH1.MM',Valid='true'))
jData[['3']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=0, e.y=3, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='shH1.MM.bort',Valid='true'))
jData[['4']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=1, e.y=0, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='WT.MM',Valid='true'))
jData[['5']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=1, e.y=1, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='WT.MM.bort',Valid='true'))
jData[['6']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=1, e.y=2, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='shH1.MM',Valid='true'))
jData[['7']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=1, e.y=3, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='shH1.MM.bort',Valid='true'))
jData[['8']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=2, e.y=0, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='shMMP2.MM',Valid='true'))
jData[['9']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=2, e.y=1, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='shMMP2.MM.bort',Valid='true'))
jData[['10']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=2, e.y=2, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='MM',Valid='true'))
jData[['11']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=2, e.y=3, type='File', name='Feature CSV Table', labels=list(Rep='1',Tx='MM.bort',Valid='true'))
jData[['12']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=3, e.y=0, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='shMMP2.MM',Valid='true'))
jData[['13']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=3, e.y=1, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='shMMP2.MM.bort',Valid='true'))
jData[['14']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=3, e.y=2, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='MM',Valid='true'))
jData[['15']] <- readJEXData(dbPath='//smph.drive.wisc.edu/smph-data/private/oncology/smlab/Tina/JEX Databases/2020.11.10-Pt902shH1MMP2bort', ds='Plate', e.x=3, e.y=3, type='File', name='Feature CSV Table', labels=list(Rep='2',Tx='MM.bort',Valid='true'))
jData <- rbindlist(jData)

# Use the data structure to read in the actual tables of data and associated "metadata"
l(x, time.col, idCols, imageDims, labelDims) %=% readJEXDataTables(jData=jData,
													  time.completeness=time.completeness,
													  sample.size=sample.size,
													  lines.with=lines.with)

# Perform preprocessing of the data to calculate additional metrics.
x3 <- preprocessData(x,
				 time.col=time.col,
				 imageDims=imageDims,
				 labelDims=labelDims,
				 doGeometry=T,
				 save.dir=save.dir)
gc()

# Save a copy of the data at this stage in case steps need to be rerun for any reason
x3b <- copy(x3)
x3 <- copy(x3b)

# Calculate drug sensitivity metrics
x3 <- calculateDrugSensitivityMetrics(x3,
							   imageDims = imageDims,
							   labelDims = labelDims,
							   maskChan = maskChan,
							   nucChan = nucChan,
							   phaseChan = phaseChan,
							   gap.max = gap.max,
							   initialCellsOnly = initialCellsOnly,
							   P.string = P.string,
							   Neg.string = Neg.string,
							   oldTxs = oldTxs,
							   newTxs = newTxs)

# Filter debris from analysis based upon circularity and size
x3[, Circ:=Geometric.Circularity.None.WholeCell]
x3.s <- x3[, list(Circ.med=median(Circ, na.rm=T), Size.med=median(Geometric.SizeIterable.Total.None.WholeCell, na.rm=T)), by=c('cId','Tx','Rep')]
data.table.plot.all(x3.s, xcol='Circ.med', ycol='Size.med', alpha=0.2, cex=0.6, randomize=T, sample.size=2000, legend.args=list(bty='n', cex=0.8), log='y', type='p', by=c('Tx'))
x3.s <- x3.s[Circ.med > 0.6 & Size.med > 1000 & Size.med < 4*1000]# TYPICALLY 0.6, 1000, and 4000 for RPMI. SEEMS 0.6, 500, 2000 is more appropriate for small primary cells.
x3 <- x3[cId %in% x3.s$cId]
data.table.plot.all(x3, xcol='Circ', ycol='Geometric.SizeIterable.Total.None.WholeCell', alpha=0.2, cex=0.6, randomize=T, sample.size=2000, legend.args=list(bty='n', cex=0.8), log='y', type='p', by=c('Tx'))

# Account for overall differences brightfield imaging between wells
normalizePhase(x3, pooled=T)

# Account for overall differences in nuclear staining intensity between wells
alignDistributions(x3,
			    col='Nuc',
			    align.by=c('Tx','Period.1'),
			    line.color.by = c('Period.1'),
			    plot.by='Tx',
			    norm.by=c('Tx'),
			    align.filter.fun=NULL,
			    plot.filter.fun=NULL,
			    data.filter.fun=NULL,
			    adj=c(0.1, 1),
			    bias=-1,
			    two.pass=F,
			    norm.method='subtraction',
			    norm.percentile = 0.5) 

#### STEP B ####
# These following 5 settings can be altered for a particular imaging setup but then no longer really need to be changed once you solidfy how image datasets are gathered.
# Histograms are analyzed for certain ranges of times to set thresholds automatically. Periods are 2-hour chunks of time (i.e., 4 timepoints if imaged at 30 min intervals).
# Less than that and histograms are less smooth making it more difficult to perform automatic threshold determination. 
periods.clust <- 1:10 # TYPICALLY 1:10. What Period.1 values should be used for finding thresholds (want early where there is at least some death but not everybody dead)
periods.rel <- 2:4 # TYPICALLY 2:4. What Period.1 values should be used to calculate the initial value of each cell's intensity (More robust than just using t=0 timepoint only) (Periods are 2 hour chunks of time consiting of 4 timepoints when imaged at 0.5 hour intervals)
sample.size <- 10000 # TYPICALLY 10000. Use this many cells to draw histograms for finding thresholds (speeds things way up compared to using ALL cells)
clust.filter.fun <- function(x){x$Nuc > 0} # TYPICALLY function(x){x$Nuc > 0}. What filter should be applied to the data when plotting histograms to find thresholds (If Nuc < 0, then the cell was probably a red blood cell, so typically want to exclude them)
p.rel <- 1-0.001^(1/3) # TYPICALLY 1-0.001^(1/3). p-value threshold for watching for relative changes in cellular phase and fluorescence. 99.9% chance of staying within threshold for 3 timepoints in a row assuming measurement error alone

#### STEP C ####
# Find NucThreshold and calculate Nuc.norm.rel (this occurs directly within x3)
# Three possible thresholds are calculated (2 absolue and 1 relative), we only choose between the 1st two absolute thresholds.
# The first represents the threshold between histogram clusters.
# The parameter 'starts' represents how many and generally where we would expect cluster means to be (i.e., an initial guess), which is then refined by the clustering algorithm.
# NucThres.clust holds the values of the thresholds between each cluster. If two clusters, then one threshold is calculated. If three clusters, then 2 thresholds are calculated, etc. and returned as a vector, thus the need to specify an index when saving the value (e.g., NucThresh.clust[1])
# NucThresh.p is used when cluster definition is difficult. A median and MAD estimate of standard deviation are used to provide an estimate of the mean and standard deviation of the population.
# The median and MAD are used because they are robust to outliers or "tails" of the distribution. 
# Then, a threshold "p-value" associated with a normal distribution is used to calculate a threshold relative to the mean, called NucThresh.p
# Therefore, if p.clust is set to 0.95, ~ 95 percent of the data would be below the calculated threshold.
# Both NucThresh.clust and NucThresh.p are calculated if needed but generally the clustering threshold is all that is necessary for the nuclear channel
l(NucThresh.clust, NucThresh.p, NucThresh.rel) %=% getThresholds(x=x3[Tx %in% c('WT.MM','WT.MM.bort','shH1.RPMI','shH1.MM.bort','shMMP2.MM','shMMP2.MM.bort','MM','MM.bort')], 
													col='Nuc.norm',
													starts=c(0,5.5),# TYPICALLY c(0, 5.5). Change this to help algorithm find appropriate cluster threshold
													seed=1223, 
													p.clust=0.95, # TYPICALLY 0.95
													p.rel=p.rel,
													periods.clust=periods.clust, 
													periods.rel=periods.rel, 
													clust.filter.fun=clust.filter.fun, 
													sample.size=sample.size)
NucThresh <- NucThresh.clust[1] # TYPICALLY NucThresh.clust[1] because we usually suggest finding 2 clusters. If 3 clusters it could be [1] or [2] based on the data.

#### STEP D ####
# Find PhaseThreshold and calculate Phase.norm.rel (this occurs directly within x3)
# PhaseThresh is generally best determined using the PhaseThresh.p threshold estimate.
l(PhaseThresh.clust, PhaseThresh.p, PhaseThresh.rel) %=% getThresholds(x=x3, 
														 col='Phase.norm',
														 starts=c(-4.0,0), # TYPICALLY c(-4.0,0). Change this to help algorithm find appropriate cluster threshold
														 seed=1223, 
														 p.clust=0.05, # TYPICALLY 0.05.
														 p.rel=p.rel,
														 periods.clust=periods.clust, 
														 periods.rel=periods.rel, 
														 clust.filter.fun=clust.filter.fun, 
														 sample.size=sample.size)
PhaseThresh <- PhaseThresh.p # TYPICALLY PhaseThresh.p

#### STEP E ####
# Calculate Death, Death.rel (this occurs directly within x3), and DeathThresh
suppressWarnings(x3[, Death:=log10(exp(Nuc.norm-NucThresh)/exp(Phase.norm-PhaseThresh)), by='cId'])
DeathThresh <- 0 # ALWAYS LEAVE AT 0 and change Phase thresh or Nuc thresh instead. The Death metric is based upon the other two relative to the Phase and Nuc thresholds
# The next 2 lines of code is only run if you would like view the distributions and where the 0 threshold lies within it.
# l(DeathThresh.clust, DeathThresh.p, DeathThresh.rel) %=% getThresholds(x=x3, 
# 														 col='Death',
# 														 starts=c(-1,1), # Change this to help algorithm find appropriate cluster threshold
# 														 seed=1223, 
# 														 p.clust=0.05, 
# 														 p.rel=p.rel,
# 														 periods.clust=periods.clust, 
# 														 periods.rel=periods.rel, 
# 														 clust.filter.fun=clust.filter.fun, 
# 														 sample.size=sample.size)
# abline(v=DeathThresh, lwd=2)

#### Data Visualization Templates ###

#### Plotting Normalized Cell Number vs Time ####
setorder(x3, Tx, cId, Time)
x3.n <- x3[, list(N=.N, Phase=sum(Phase.norm >= PhaseThresh), Nuc=sum(Nuc.norm <= NucThresh), Death=sum(Death <= DeathThresh)), by=c('Tx','Rep','Time')]
setorder(x3.n, Tx, Time)
x3.n[, ':='(N0=N[1], Phase0=Phase[1], Nuc0=Nuc[1], Death0=Death[1]), by=c('Tx','Rep')]
x3.n[, ':='(N.norm=N/N0, Phase.N=Phase/N, Nuc.N=Nuc/N, Death.N=Death/N), by=c('Tx','Rep')]
# Total Cells
data.table.plot.all(x3.n, xcol='Time', ycol='N.norm', type='l', by=c('Tx','Rep'), ylim=c(0,3), legend.args=list(x='bottomleft', bty='n', cex=0.8))
# Number of cells above the Phase threshold (i.e., ~living) according to the Phase channel
data.table.plot.all(x3.n, xcol='Time', ycol='Phase.N', type='l', by='Tx', ylim=c(0,3), legend.args=list(x='bottomleft', bty='n', cex=0.8))
# Number of cells below the Nuc threshold (i.e., ~living) according to nuclear staining
data.table.plot.all(x3.n, xcol='Time', ycol='Nuc.N', type='l', by='Tx', ylim=c(0,3), legend.args=list(x='bottomleft', bty='n', cex=0.8))
# Number of cells below the Death threshold (i.e., ~living) according to combined phase and nuclear staining information
data.table.plot.all(x3.n, xcol='Time', ycol='Death.N', type='l', by='Tx', ylim=c(0,3), legend.args=list(x='bottomleft', bty='n', cex=0.8))

#### Generating histograms of variabls over time ####
makeDrugSensitivityHistograms(x3, makePhase=T, makeNuc=F, makeDeath=F, PhaseThresh=PhaseThresh, NucThresh=NucThresh, DeathThresh=DeathThresh, save.dir=save.dir)
makeDrugSensitivityHistograms(x3, makeNuc=T, makePhase=F, makeDeath=F, PhaseThresh=PhaseThresh, NucThresh=NucThresh, DeathThresh=DeathThresh, save.dir=save.dir)
makeDrugSensitivityHistograms(x3, makeDeath=T, makePhase=F, makeNuc=F, PhaseThresh=PhaseThresh, NucThresh=NucThresh, DeathThresh=DeathThresh, save.dir=save.dir)

#### Plotting cell survival curves based on different cell death criteria ####
# Two functions are provided for assessing cell events
# getEventStatus - Determines when a cell meets the specified criteria for at least n.true number of timepoints in a row
# getEventRecoveryStatus - Determines if a cell is first observed to meet the criteria for n.true number of timepoints, THEN (i.e., afterward) is observed to NOT meet the criteria for n.true number of timepoints (i.e., the cell "recovers" from the initial event) (i.e., a recovery event is observed)
library(survminer)
library(survival)
# Order the data
setorderv(x3, c('cId',labelDims,'Time'))
# Print out the 'Tx' conditions within the dataset for reference
paste("'", paste(uniqueo(x3$Tx), collapse="','"), "'", sep='') # Print out a comma separated list of all the Tx conditions and use below
# Create a variable with the Tx's that you want to plot (CHECK SPELLING AND CAPITALIZATION!!!)
Txs <- c('WT.MM','WT.MM.bort')

# Death - Nuc Metric - Identify cell death events (i.e., when Nuc.norm > NucThresh) and generate an associated survival probability curve
x.surv <- x.m[, getEventStatus(Time, Nuc.norm > NucThresh, n.true=3, suffix=''), by=c('Tx','cId')]
plotSurvivalCurve(x.surv=x.surv,
			   x=x3,
			   ylab='Survival Probability (Nuc Metric)',
			   save.file='Survival Probability (Nuc Metric).png',
			   flip=F,
			   save.plot=T,
			   Txs=Txs,
			   viability.y=-1,
			   pval.y=0.2,
			   xlim=c(0,50),
			   ylim=c(0,1),
			   legend.cex=0.6,
			   width=4,
			   height=4,
			   save.dir=save.dir)

# Lysis - Nuc Metric - Identify cell lysis events (i.e., when Nuc.norm > NucThresh as the cell condenses then Nuc.norm <= NucThresh when nuclear material is released into the surrounding media) and generate an associated lysis probability curve
# The 'flip' parameter is used to flip (1-prob) the y axis of the probability plot (i.e., prob=1 changes to a prob=(1-1)=0 or 0.75 changes to (1-0.75)=0.25)
x.surv <- x.m[, getEventRecoveryStatus(Time, Nuc.norm > NucThresh, n.true=3, suffix=''), by=c('Tx','cId')]
plotSurvivalCurve(x.surv=x.surv,
			   x=x3,
			   ylab='Recovery Probability (Nuc Metric) [Lysis]',
			   save.file='Recovery Probability (Nuc Metric) [Lysis].png',
			   flip=T,
			   save.plot=T,
			   Txs=Txs,
			   viability.y=-1,
			   pval.y=1,
			   xlim=c(0,50),
			   ylim=c(0,1),
			   legend.cex=0.8,
			   width=4,
			   height=4,
			   save.dir=save.dir)
# Death - Phase Metric - Identify cell death events based solely on phase information (i.e., Phase.norm < PhaseThresh) and generate an associated survival probability curve
x.surv <- x.m[, getEventStatus(Time, Phase.norm < PhaseThresh, n.true=3, suffix=''), by=c('Tx','cId')]
plotSurvivalCurve(x.surv=x.surv,
			   x=x3,
			   ylab='Survival Probability (Phase Metric)',
			   save.file='Survival Probability (Phase Metric).png',
			   flip=F,
			   save.plot=T,
			   Txs=Txs,
			   viability.y=0.5,
			   legend.cex=0.6,
			   width=4,
			   height=4,
			   pval.y=0.2,
			   xlim=c(0,50),
			   ylim=c(0,1),
			   legend.cex=0.6,
			   width=4,
			   height=4,
			   save.dir=save.dir)

# Recovery - Phase Metric - Identify cell recovery events based solely on phase information (i.e., Phase.norm < PhaseThresh then later Phase.norm >= PhaseThresh) and generate an associated recovery probability curve
x.surv <- x.m[, getEventRecoveryStatus(Time, Phase.norm < PhaseThresh, n.true=3, suffix=''), by=c('Tx','cId')]
plotSurvivalCurve(x.surv=x.surv,
			   x=x3,
			   ylab='Recovery Probability (Phase Metric)',
			   save.file='Recovery Probability (Phase Metric).png',
			   flip=T,
			   save.plot=T,
			   Txs=Txs,
			   viability.y=0.5,
			   pval.y=0.2,
			   xlim=c(0,50),
			   ylim=c(0,1),
			   legend.cex=0.6,
			   width=4,
			   height=4,
			   save.dir=save.dir)

# Death - Death Metric - Identify cell death events based on both phase and nuclear information (i.e., Death > DeathThresh) and generate an associated survival probability curve
x.surv <- x.m[, getEventStatus(Time, Death > DeathThresh, n.true=3, suffix=''), by=c('Tx','cId')]
plotSurvivalCurve(x.surv=x.surv,
			   x=x3,
			   ylab='Survival Probability (Death Metric)',
			   save.file='Survival Probability (Death Metric).png',
			   flip=F,
			   save.plot=T,
			   Txs=Txs,
			   viability.y=-1,
			   pval.y=0.2,
			   xlim=c(0,50),
			   ylim=c(0,1),
			   legend.cex=0.8,
			   width=4,
			   height=4,
			   save.dir=save.dir)

# Recovery - Death Metric - Identify cell recovery events based on both phase and nuclear information (i.e., Death < DeathThresh then later Death <= DeathThresh) and generate an associated recovery probability curve
x.surv <- x.m[, getEventRecoveryStatus(Time, Death > DeathThresh, n.true=3, suffix=''), by=c('Tx','cId')]
plotSurvivalCurve(x.surv=x.surv,
			   x=x3,
			   ylab='Recovery Probability (Death Metric)',
			   save.file='Recovery Probability (Death Metric).png',
			   flip=T,
			   save.plot=T,
			   Txs=Txs,
			   viability.y=-1,
			   pval.y=1,
			   xlim=c(0,50),
			   ylim=c(0,1),
			   legend.cex=0.8,
			   width=4,
			   height=4,
			   save.dir=save.dir)
