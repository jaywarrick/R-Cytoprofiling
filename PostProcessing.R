###########
# Cell Area
###########

breaks=c(-1000, seq(-4,4,0.5), 1000)
wt <- shinyData[Class == 'WT']$Geometric.SizeIterable_WholeCell_None
mt <- shinyData[Class == 'MT']$Geometric.SizeIterable_WholeCell_None
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
t.test(mt, wt, na.rm=T)
(median(mt)-median(wt))
hist(wt, col=cwt, breaks=breaks, xlim=c(-3,3), freq=TRUE, main='', xlab='Standardized Cell Area')
hist(mt, col=cmt, breaks=breaks, xlim=c(-3,3), freq=TRUE, add=T)
wtd <- density(wt, n=1000)
mtd <- density(mt, n=1000)
plot(mtd, col='blue', xlim=c(-3,3), main='', xlab='Standardized Cell Area')
lines(wtd, col='red')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

###########
# Integrate Intensity
###########

breaks=c(-1000, seq(-4,4,0.5), 1000)
wt <- shinyData[Class == 'WT']$Stats.Sum_WholeCell_395X455M
mt <- shinyData[Class == 'MT']$Stats.Sum_WholeCell_395X455M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtd <- density(wt, n=1000)
mtd <- density(mt, n=1000)
plot(wtd, col='red', xlim=c(-3,3), main='', xlab='Standardized Cell Area')
lines(mtd, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)


###########
# Central Moment (2,0) & (0,2)
###########
wtX <- shinyData[Class == 'WT']$ImageMoments.CentralMoment20_WholeCell_395X455M_minus_650X705M
mtX <- shinyData[Class == 'MT']$ImageMoments.CentralMoment20_WholeCell_395X455M_minus_650X705M
wtY <- shinyData[Class == 'WT']$ImageMoments.CentralMoment02_WholeCell_395X455M_minus_650X705M
mtY <- shinyData[Class == 'MT']$ImageMoments.CentralMoment02_WholeCell_395X455M_minus_650X705M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdX <- density(wtX, n=1000)
mtdX <- density(mtX, n=3000)
wtdY <- density(wtY, n=1000)
mtdY <- density(mtY, n=3000)
plot(wtdX, col='red', xlim=c(-3,3), main='', xlab='CMX [Nuc-p65]')
lines(mtdX, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)
plot(wtdY, col='red', xlim=c(-3,3), main='', xlab='CMY [Nuc-p65]')
lines(mtdY, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

###########
# Hu Moment 1
###########
wtHu1 <- shinyData[Class == 'WT']$RHu1_WholeCell_395X455M_minus_650X705M
mtHu1 <- shinyData[Class == 'MT']$RHu1_WholeCell_395X455M_minus_650X705M
headshinyData$ImageMoments.CentralMoment03_WholeCell_395X455M_minus_650X705M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=1000)
mtdHu1 <- density(mtHu1, n=3000)
plot(wtdHu1, col='red', xlim=c(-3,3), main='', xlab='Hu1')
lines(mtdHu1, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

wtHu1 <- shinyData[Class == 'WT']$RHu1_WholeCell_395X455M
mtHu1 <- shinyData[Class == 'MT']$RHu1_WholeCell_395X455M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=1000)
mtdHu1 <- density(mtHu1, n=3000)
plot(wtdHu1, col='red', xlim=c(-3,3), main='', xlab='Hu1')
lines(mtdHu1, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

wtHu1 <- shinyData[Class == 'WT']$RHu1_WholeCell_650X705M
mtHu1 <- shinyData[Class == 'MT']$RHu1_WholeCell_650X705M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=1000)
mtdHu1 <- density(mtHu1, n=3000)
plot(wtdHu1, col='red', xlim=c(-3,3), main='', xlab='Hu1')
lines(mtdHu1, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

###########
# Sum
###########
wtHu1 <- shinyData[Class == 'WT']$Stats.Sum_WholeCell_395X455M
mtHu1 <- shinyData[Class == 'MT']$Stats.Sum_WholeCell_395X455M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=1000)
mtdHu1 <- density(mtHu1, n=3000)
plot(wtdHu1, col='red', xlim=c(-3,3), main='', xlab='Hu1')
lines(mtdHu1, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

wtHu1 <- shinyData[Class == 'WT']$Stats.Sum_WholeCell_650X705M
mtHu1 <- shinyData[Class == 'MT']$Stats.Sum_WholeCell_650X705M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=1000)
mtdHu1 <- density(mtHu1, n=3000)
plot(wtdHu1, col='red', xlim=c(-3,3), main='', xlab='Hu1')
lines(mtdHu1, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

###########
# Sum of Logs
###########
wtHu1 <- shinyData[Class == 'WT']$Stats.SumOfLogs_WholeCell_395X455M
mtHu1 <- shinyData[Class == 'MT']$Stats.SumOfLogs_WholeCell_395X455M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=1000)
mtdHu1 <- density(mtHu1, n=3000)
plot(wtdHu1, col='red', xlim=c(-3,3), main='', xlab='Hu1')
lines(mtdHu1, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)

wtHu1 <- shinyData[Class == 'WT']$Stats.SumOfLogs_WholeCell_650X705M
mtHu1 <- shinyData[Class == 'MT']$Stats.SumOfLogs_WholeCell_650X705M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=1000)
mtdHu1 <- density(mtHu1, n=3000)
plot(wtdHu1, col='red', xlim=c(-3,3), main='', xlab='Hu1')
lines(mtdHu1, col='blue')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)


###########
# ZDot
###########
wtHu1 <- shinyData[Class == 'WT']$ZernikeDot33_THISwFIXED_WholeCell_485X525M_dot_650X705M
mtHu1 <- shinyData[Class == 'MT']$ZernikeDot33_THISwFIXED_WholeCell_485X525M_dot_650X705M
cmt <- rgb(0,0,1,0.8)
cwt <- rgb(1,0,0,0.8)
wtdHu1 <- density(wtHu1, n=15000)
mtdHu1 <- density(mtHu1, n=1300000)
plot(mtdHu1, col='blue', xlim=c(-6,10), main='', xlab='Hu1')
lines(wtdHu1, col='red')
legend('topright', legend=c('MT','WT'), col=c('blue','red'), lty=1)












cmt <- rgb(0,0,1,0.4)
cwt <- rgb(1,0,0,0.4)
x4$Measurement <- paste(x4$Measurement, x4$MaskChannel, x4$ImageChannel, sep='_')
x <- x2b[MaskChannel=='WholeCell' & net.imagej.ops.Ops.Stats.Mean > 500, c('Id','ImRow','ImCol','Class','ImageChannel','net.imagej.ops.Ops.Stats.Mean'), with=FALSE]
x$col <- cwt
x[Class=='MT']$col <- cmt
x <- x[sample.int(nrow(x))]
x$cId <- paste0(x$Id, ' RCClass[', x$ImRow, ',', x$ImCol, ',', x$Class, ']')
x$cId <- paste0(x$Id, ' RCClass[', x$ImRow, ',', x$ImCol, ',', x$Class, ']')
missingA <- x[ImageChannel=='485X525M',][!(x[ImageChannel=='485X525M',]$cId %in% x[ImageChannel=='650X705M',]$cId),]$cId
missingB <- x[ImageChannel=='650X705M',][!(x[ImageChannel=='650X705M',]$cId %in% x[ImageChannel=='485X525M',]$cId),]$cId
plot(x[!(cId %in% c(missingA, missingB)) & ImageChannel=='485X525M']$net.imagej.ops.Ops.Stats.Mean-500, x[!(cId %in% c(missingA, missingB)) & ImageChannel=='650X705M']$net.imagej.ops.Ops.Stats.Mean-500, log='xy', cex=0.5, pch=20, col=x$col, xlim=c(1,1500),ylim=c(1,50), xlab='CMFDA Intensity [au]', ylab='p65 Staining Intensity [au]')
