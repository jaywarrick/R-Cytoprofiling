rm(list=ls())

source('D:/GitHub/R-General/.Rprofile')
source('D:/GitHub/R-Cytoprofiling/20161013_ImageFlow/PreprocessingHelpers.R')

tList <- getTableListFromDB(db='G:/JEX Databases/R01 Renewal', ds='20161216 - Bead Staining', x=0, y=0, objectName='Output CSV Table Intensity', isArff=F, storeFilePath=F, class=NULL, assignClass=F, expt=NULL, repl=NULL, sampleSize=NULL, colsToRemove = c(), cIdCols = c(), fsep='\\\\')

x <- rbindlist(tList, use.names = T)

x <- replaceSubStringInAllRowsOfCol(x, '$', '.', 'Measurement')
x <- replaceSubStringInAllRowsOfCol(x, '390 X 440', 'Blue', 'ImageChannel')
x <- replaceSubStringInAllRowsOfCol(x, '648 X 684', 'Red', 'ImageChannel')

x <- x[Measurement == 'net.imagej.ops.Ops.Stats.Sum' & ImageChannel != 'DIC']
x[, Measurement:='Mean']
x1 <- reorganize(x, measurementCols = 'Measurement,ImageChannel,MaskChannel')
setorder(x1, Beads, Tx, Well, Rep, Id)

x1 <- removeIncompleteRows(x1)

x1 <- x1[Rep=='A']
x1[,Rep:=NULL]

x1[,NucRatio:=Mean_Red_Nuc/(Mean_Red_WholeCell)]

duh <- x1[,list(meanNucRatio=mean(NucRatio)), by=.(Beads,Tx)]

# NoMin <- duh[Beads=='No' & Tx=='TNF']$meanNucRatio
# NoRange <- duh[Beads=='No' & Tx=='GST']$meanNucRatio - NoMin
# x1[Beads=='No', NucRatio:=(NucRatio - NoMin) / NoRange]
# 
# YesMin <- duh[Beads=='Yes' & Tx=='TNF']$meanNucRatio
# YesRange <- duh[Beads=='Yes' & Tx=='GST']$meanNucRatio - YesMin
# x1[Beads=='Yes', NucRatio:=(NucRatio - YesMin) / YesRange]

##### PTR1
pdf('G:/R01 Renewal/DensityPlot.pdf', width=6, height=4)
y <- x1[Beads=='Yes' & Tx=='PTR1']
n <- x1[Beads=='No' & Tx=='PTR1']
toPlot <- density(y$NucRatio, n=40)
toPlot2 <- density(n$NucRatio, n=40)
plot(toPlot$x-0.0285, toPlot$y, xlab='Nuclear Fraction', ylab='Histogram Density', lwd=2, type='l')
lines(toPlot2$x, toPlot2$y, lwd=2, col='red', lty=1)
dev.off()
wilcox.test(x=y$NucRatio, y=n$NucRatio)
wilcox.test(x=y$Mean_Red_Nuc, y=n$Mean_Red_Nuc)
wilcox.test(x=y$Mean_Red_Cyt, y=n$Mean_Red_Cyt)
plot(y$Mean_Red_Nuc, y$Mean_Red_Cyt, pch=20, col=adjustcolor('red', alpha.f=0.5))
points(n$Mean_Red_Nuc, n$Mean_Red_Cyt, pch=20, col=adjustcolor('blue', alpha.f=0.5))

##### TNF
y <- x1[Beads=='Yes' & Tx=='TNF']
n <- x1[Beads=='No' & Tx=='TNF']
plot(density(y$NucRatio, breaks=40))
lines(density(n$NucRatio, breaks=40), col='red')
wilcox.test(x=y$NucRatio, y=n$NucRatio)
wilcox.test(x=y$Mean_Red_Nuc, y=n$Mean_Red_Nuc)
wilcox.test(x=y$Mean_Red_Cyt, y=n$Mean_Red_Cyt)
plot(y$Mean_Red_Nuc, y$Mean_Red_Cyt, pch=20, col=adjustcolor('red', alpha.f=0.5))
points(n$Mean_Red_Nuc, n$Mean_Red_Cyt, pch=20, col=adjustcolor('blue', alpha.f=0.5))

##### GST
y <- x1[Beads=='Yes' & Tx=='GST']
n <- x1[Beads=='No' & Tx=='GST']
plot(density(y$NucRatio, breaks=40))
lines(density(n$NucRatio, breaks=40), col='red')
wilcox.test(x=y$NucRatio, y=n$NucRatio)
wilcox.test(x=y$Mean_Red_Nuc, y=n$Mean_Red_Nuc)
wilcox.test(x=y$Mean_Red_Cyt, y=n$Mean_Red_Cyt)
plot(y$Mean_Red_Nuc, y$Mean_Red_Cyt, pch=20, col=adjustcolor('red', alpha.f=0.5))
points(n$Mean_Red_Nuc, n$Mean_Red_Cyt, pch=20, col=adjustcolor('blue', alpha.f=0.5))

##### PTR1 vs TNF
y <- x1[Beads=='No' & Tx=='PTR1']
n <- x1[Beads=='No' & Tx=='TNF']
plot(density(y$NucRatio, breaks=40))
lines(density(n$NucRatio, breaks=40), col='red')
wilcox.test(x=y$NucRatio, y=n$NucRatio)
wilcox.test(x=y$Mean_Red_Nuc, y=n$Mean_Red_Nuc)
wilcox.test(x=y$Mean_Red_Cyt, y=n$Mean_Red_Cyt)
plot(y$Mean_Red_Nuc, y$Mean_Red_Cyt, pch=20, col=adjustcolor('red', alpha.f=0.5))
points(n$Mean_Red_Nuc, n$Mean_Red_Cyt, pch=20, col=adjustcolor('blue', alpha.f=0.5))

##### PTR1 vs TNF With Beads
y <- x1[Beads=='Yes' & Tx=='PTR1']
n <- x1[Beads=='Yes' & Tx=='TNF']
plot(density(y$NucRatio, breaks=40))
lines(density(n$NucRatio, breaks=40), col='red')
wilcox.test(x=y$NucRatio, y=n$NucRatio)
wilcox.test(x=y$Mean_Red_Nuc, y=n$Mean_Red_Nuc)
wilcox.test(x=y$Mean_Red_Cyt, y=n$Mean_Red_Cyt)
plot(y$Mean_Red_Nuc, y$Mean_Red_Cyt, pch=20, col=adjustcolor('red', alpha.f=0.5))
points(n$Mean_Red_Nuc, n$Mean_Red_Cyt, pch=20, col=adjustcolor('blue', alpha.f=0.5))


##### Cell insets #####
x1[Beads=='No' & Id=='17' & Tx=='PTR1' & Well==2]
x1[Beads=='No' & Id=='22' & Tx=='PTR1' & Well==2]

x1[Beads=='No' & Tx=='PTR1' & Well==2, c('Id','NucRatio'), with=F]
