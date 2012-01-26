
source("SQLShareLib.R")
#plot ecoli distribution with error bars

#efracs <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE, username="wenjunh")
#ufracs <- fetchAllSurfResidues("ecoli", -1., normalize=TRUE, username="wenjunh")
#ems <- apply(efracs, MARGIN=2, mean)
 #eci <- apply(efracs, MARGIN=2, FUN=function(x) matrix(c(quantile(x, c(0.975)) - median(x), median(x) - quantile(x, c(0.025))), byrow=T, nrow=2))
#ums <- apply(ufracs, MARGIN=2, mean)
#uci <- apply(ufracs, MARGIN=2, FUN=function(x) matrix(c(quantile(x, c(0.975)) - median(x), median(x) - quantile(x, c(0.025))), byrow=T, nrow=2))


#tiff("ecolidist_ci.tiff", width=3.3*500, height=2.5*500, pointsize=8, res=500)
#par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
#yy <- matrix(c(ems,ums), byrow=T, nrow=2)
#print(yy)
#ee.up <- matrix(c(eci[1,], uci[1,]), byrow=T, nrow=2)
#ee.low <- matrix(c(eci[2,], uci[2,]), byrow=T, nrow=2)
#barx <- barplot(yy, col=c("gray90", "gray30"), main="", xlab="Amino Acid", ylab="Density", beside=T, names.arg=aalist.sh,ylim = c(0.00,0.14))#ylim=c(0,ee.up[which.max(yy)] + max(yy)))
#error.bar(barx, yy, upper=ee.up, lower=ee.low,length=0.01)
#legend("topright", col=c("gray90", "gray30"), legend=c("E. Coli Surface", "E. Coli All"), ncol=1, pch=15)
#graphics.off()

#evs <- apply(efracs, MARGIN=2, var)
#uvs <- apply(ufracs, MARGIN=2, var)



tiff("ecolidist.tiff", width=3.3*500, height=2.5*500, pointsize=8, res=500)
#par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
#yy <- matrix(c(ems,ums), byrow=T, nrow=2)
#ee <- matrix(c(sqrt(evs / nrow(efracs)), sqrt(uvs / nrow(ufracs))), byrow=T, nrow=2)
#print(ee)
#barx <- barplot(yy, col=c("gray90", "gray30"), main="", xlab="Amino Acid", ylab="Density", beside=T, names.arg=aalist.sh,ylim = c(0.00,0.14))#ylim=c(0,ee[which.max(yy)] + max(yy)))
#error.bar(barx, yy, ee, length=0.01)
#legend("topright", col=c("gray90", "gray30"), legend=c("E. Coli Surface", "E. Coli All"), ncol=1, pch=15)
#graphics.off()




sql <- paste("select * FROM [whitead@washington.edu].[GroEL_counts.csv]")
rawData <- fetchdata(sql)


groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist[1, ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist[2, ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))

b <- array(0,20)

cairo_pdf('groDist_SE.pdf', width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.6)
barx <- barplot(as.matrix(groDist), col=c("gray15","gray75"), main="", xlab="Amino Acid", ylab="", beside=T, names.arg=aalist.sh, ylim = c(0.00,0.20))
legend("topright", col=c("gray15","gray75"), legend=c("GroEL-GroES (Closed)","GroEL (Open)"), pch=15)
graphics.off()


#plot2 The plot of 5 species

hspRawCount <- fetchHSPResCount("hsp","in",username="whitead")

for (i in 1:nrow(hspRawCount)) {
  hspRawCount[i,] <- hspRawCount[i,] / sum(hspRawCount[i,])
}
hspRawCount <- as.matrix(hspRawCount)

hspDist <- matrix(0,nrow(hspRawCount)+1, ncol(hspRawCount))
rownames(hspDist) <- c("E.Coli","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic CCT")
colnames(hspDist) <- colnames(hspRawCount)

hspDist[2,] <- hspRawCount["1SX4",]
hspDist[3,] <- hspRawCount["1WE3",]
hspDist[4,] <- hspRawCount["3KFB",]
hspDist[5,] <- hspRawCount["2CG9",]
hspDist[6,] <- hspRawCount["3P9D",]

#various statistics
print(sum(hspRawCount["3KFB", c("ASP", "GLU", "ARG", "HIS")]))

cutoff <- 0.3
psurf <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE, "wenjunh")
hspDist[1,] <- apply(psurf, MARGIN=2, FUN=mean)

cairo_pdf('HSP_Protein_Dist.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
barplot(hspDist, beside=TRUE, col=c('blue','gray75','gray60','gray45','gray30','gray15'), ylim=c(0.0,0.3), names.arg=aalist.sh,)
legend("topright", col=c('blue','gray75','gray60','gray45','gray30','gray15'), legend=c("E.Coli","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic CCT"), pch=rep(15,6), cex=0.8)
graphics.off()
