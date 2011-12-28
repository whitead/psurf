
source("SQLShareLib.R")
#plot ecoli distribution with error bars

efracs <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE, username="wenjunh")
ufracs <- fetchAllSurfResidues("ecoli", -1., normalize=TRUE, username="wenjunh")
ems <- apply(efracs, MARGIN=2, mean)
eci <- apply(efracs, MARGIN=2, FUN=function(x) matrix(c(quantile(x, c(0.975)) - median(x), median(x) - quantile(x, c(0.025))), byrow=T, nrow=2))
ums <- apply(ufracs, MARGIN=2, mean)
uci <- apply(ufracs, MARGIN=2, FUN=function(x) matrix(c(quantile(x, c(0.975)) - median(x), median(x) - quantile(x, c(0.025))), byrow=T, nrow=2))


tiff("ecolidist_ci.tiff", width=3.3*500, height=2.5*500, pointsize=8, res=500)
par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
yy <- matrix(c(ems,ums), byrow=T, nrow=2)
print(yy)
ee.up <- matrix(c(eci[1,], uci[1,]), byrow=T, nrow=2)
ee.low <- matrix(c(eci[2,], uci[2,]), byrow=T, nrow=2)
barx <- barplot(yy, col=c("gray90", "gray30"), main="", xlab="Amino Acid", ylab="Density", beside=T, names.arg=aalist.sh,ylim = c(0.00,0.14))#ylim=c(0,ee.up[which.max(yy)] + max(yy)))
error.bar(barx, yy, upper=ee.up, lower=ee.low,length=0.01)
legend("topright", col=c("gray90", "gray30"), legend=c("E. Coli Surface", "E. Coli All"), ncol=1, pch=15)
graphics.off()

evs <- apply(efracs, MARGIN=2, var)
uvs <- apply(ufracs, MARGIN=2, var)



tiff("ecolidist.tiff", width=3.3*500, height=2.5*500, pointsize=8, res=500)
par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
yy <- matrix(c(ems,ums), byrow=T, nrow=2)
ee <- matrix(c(sqrt(evs / nrow(efracs)), sqrt(uvs / nrow(ufracs))), byrow=T, nrow=2)
print(ee)
barx <- barplot(yy, col=c("gray90", "gray30"), main="", xlab="Amino Acid", ylab="Density", beside=T, names.arg=aalist.sh,ylim = c(0.00,0.14))#ylim=c(0,ee[which.max(yy)] + max(yy)))
error.bar(barx, yy, ee, length=0.01)
legend("topright", col=c("gray90", "gray30"), legend=c("E. Coli Surface", "E. Coli All"), ncol=1, pch=15)
graphics.off()

#plot open/close with error bars


#sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
#rawData <- fetchdata(sql)
#groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
#groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)
#
#groDist[1, ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
#groDist[2, ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))
#
#tiff("groELDist.tiff", width=3.3*500, height=2.5*500, pointsize=8, res=500)
#par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
#barplot(as.matrix(groDist[c(2,1),]), col=c("gray90", "gray25"), main="", xlab="Amino Acid", ylab="Density", beside=T, names.arg=aalist.sh)
#legend("topright", col=c("gray90", "gray25"), legend=c("GroEL Open", "GroEL Close"), ncol=1, pch=15)
#graphics.off()
