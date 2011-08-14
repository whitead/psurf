
source("SQLShareLib.R")
#plot ecoli distribution with error bars

efracs <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE, username="wenjunh")
ems <- apply(efracs, MARGIN=2, mean)
evs <- apply(efracs, MARGIN=2, var)


tiff("ecolidist.tiff", width=3.3*500, height=2.5*500, pointsize=8, res=500)
par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
yy <- matrix(ems, byrow=T, nrow=1)
ee <- matrix(sqrt(evs), byrow=T, nrow=1)
barx <- barplot(yy, col=c("dark gray"), main="", xlab="Amino Acid", ylab="Density", beside=F,
                ylim=c(0,ee[which.max(yy)] + max(yy)), names.arg=aalist.sh)
error.bar(barx, yy, ee, length=0.02)
legend("topright", col=c("dark gray"), legend=c("E. Coli"), ncol=1, pch=15)
graphics.off()

#plot open/close with error bars


sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))

tiff("groELDist.tiff", width=3.3*500, height=2.5*500, pointsize=8, res=500)
par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
barplot(as.matrix(groDist[1:2,]), col=c("dark gray", "light gray"), main="", xlab="Amino Acid", ylab="Count", beside=T, names.arg=aalist.sh)
legend("topright", col=c("dark gray", "light gray"), legend=c("GroEL Close", "GroEL Open"), ncol=1, pch=15)
graphics.off()
