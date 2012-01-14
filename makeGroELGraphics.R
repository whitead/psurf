
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
barx <- barplot(yy, col=c("gray90", "gray30"), main="", xlab="Amino Acid", ylab="Density", beside=T, names.arg=aalist.sh,ylim = c(0.00,0.14))#ylim=c(0,ee[which.max(yy)] + max(yy)))
error.bar(barx, yy, ee, length=0.01)
legend("topright", col=c("gray90", "gray30"), legend=c("E. Coli Surface", "E. Coli All"), ncol=1, pch=15)
graphics.off()


#plot1(a) ecoli distribution with error bars

#efracs <- fetchAllSurfResidues("ecoli40", cutoff, normalize=TRUE, username="wenjunh")
#
#ems <- apply(efracs, MARGIN=2, mean)
#eci <- apply(efracs, MARGIN=2, FUN=function(x) matrix(c(quantile(x, c(0.975)) - median(x), median(x) - quantile(x, c(0.025))), byrow=T, nrow=2))
#ese <- apply(efracs, MARGIN=2, FUN=function(x) matrix(c(median(x) + sd(x) / sqrt(nrow(efracs)) , median(x) - sd(x) / sqrt(nrow(efracs))), byrow=T, nrow=2))
#
#sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
#rawData <- fetchdata(sql)
#sql2 <- paste("select * FROM [wenjunh@washington.edu].[groel_outsurfres_count.csv]")
#rawData2 <- fetchdata(sql2)
#
#groDist <- empty.df(rawData[[1]][[1]][-1], c('ecoli', rawData[[1]][[2]][1], rawData[[1]][[3]][1], rawData2[[1]][[2]][1]))
#groDist[1,] <- ems
#groDist[2:3,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)
#groDist[4,] <- matrix(unlist(lapply(rawData2[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=1) 
#
#groDist[2, ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
#groDist[3, ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))
#groDist[4, ] <- as.double(groDist["GroEL_Out", ]) / sum(as.double(groDist["GroEL_Out", ]))
#
#b <- array(0,20)
#ese.up <- matrix(c(ese[1,],b,b,b), byrow=T, nrow=4)
#ese.low <- matrix(c(ese[2,],b,b,b), byrow=T, nrow=4)
#
#cairo_pdf('groDist_SE.pdf', width=3.42, height=2.58)
##tiff("groDist_SE.tiff", width=3.42, height=2.58) 
#par(family="LMSans10", cex.axis=0.65)
#barx <- barplot(as.matrix(groDist), col=c("gray90","gray75","gray60","gray45"), main="", xlab="Amino Acid", ylab="Density", beside=T, names.arg=aalist.sh, ylim = c(0.00,0.20))
##error.bar(barx, groDist, upper=ese.up, lower=ese.low)#, length=0.01)
##legend("topright", col=c("gray90","gray75","gray60","gray45", legend=c("E.Coli Surface","GroEL Close","GroEL Open","GroEL Out")))
#graphics.off()


#plot2 The plot of 5 species

hspRawCount <- fetchHSPResCount("hsp","In","wenjunh")

for (i in 1:nrow(hspRawCount)) {
  hspRawCount[i,] <- hspRawCount[i,] / sum(hspRawCount[i,])
}
hspRawCount <- as.matrix(hspRawCount)

hspDist <- matrix(0,nrow(hspRawCount)+1, ncol(hspRawCount))
rownames(hspDist) <- c("E.Coli","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic HSP")
colnames(hspDist) <- colnames(hspRawCount)

hspDist[2,] <- hspRawCount["1SX4",]
hspDist[3,] <- hspRawCount["1WE3",]
hspDist[4,] <- hspRawCount["3KFB",]
hspDist[5,] <- hspRawCount["2CG9",]
hspDist[6,] <- hspRawCount["3P9D",]

cutoff <- 0.3
psurf <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE, "wenjunh")
hspDist[1,] <- apply(psurf, MARGIN=2, FUN=mean)

cairo_pdf('HSP_Protein_Dist.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
barplot(hspDist, beside=TRUE, col=c('gray90','gray75','gray60','gray45','gray30','gray15'), ylim=c(0.0,0.3), names.arg=aalist.sh,ylab="Fraction")
legend("topright", col=c('gray90','gray75','gray60','gray45','gray30','gray15'), legend=c("E.Coli","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic HSP"), pch=rep(15,6), cex=0.8)
graphics.off()


#plot ##--supplement: plot of E.Coli distribution with different surface cutoffs
SurfCutOff <- c(0,0.1,0.2,0.3,0.4,0.5)

psurfplot <- matrix(0,6,20)
rownames(psurfplot) <- c("All","10","20","30","40","50")
colnames(psurfplot) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

for (i in 1:nrow(psurfplot)) {
  temp <- fetchAllSurfResidues("ecoli",SurfCutOff[i],normalize=TRUE,"wenjunh")
  for (j in 1:ncol(temp)) {
    psurfplot[i,j] <- mean(temp[,j])
  }
}

cairo_pdf("ecoli_cutoff.pdf",width=8,height=5)
barplot(psurfplot,beside=T,col=c('blue','gray75','gray60','gray45','gray30','gray15'),ylim=c(0,0.14),ylab="Fraction")
legend("topright",col=c('blue','gray75','gray60','gray45','gray30','gray15'),legend=c("All","0.1 Cutoff","0.2 Cutoff","0.3 Cutoff","0.4 Cutoff","0.5 Cutoff"),pch=15)
dev.off()

