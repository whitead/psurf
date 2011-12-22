
source("SQLShareLib.R")

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
  hspRawCount[i,] <- hspRawCount[i,] / max(hspRawCount[i,])
}
hspRawCount <- as.matrix(hspRawCount)

hspDist <- matrix(0,nrow(hspRawCount)+1, ncol(hspRawCount))
rownames(hspDist) <- c("E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic HSP","E.Coli")
colnames(hspDist) <- colnames(hspRawCount)

hspDist[1,] <- hspRawCount["1SX4",]
hspDist[2,] <- hspRawCount["1WE3",]
hspDist[3,] <- hspRawCount["3KFB",]
hspDist[4,] <- hspRawCount["2CG9",]
hspDist[5,] <- hspRawCount["3P9D",]

cutoff <- 0.3
psurf <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE, "wenjunh")
hspDist[6,] <- apply(psurf, MARGIN=2, FUN=mean)

cairo_pdf('HSP_Protein_Dist.pdf',width=6.84, height=2.58)
par(family='LMSans10', cex.axis=0.65, ps=9)
barplot(hspDist, beside=TRUE, col=c('gray90','gray75','gray60','gray45','gray30','gray15'), ylim=c(-1.0,2.0))
legend("topright", col=c('gray90','gray75','gray60','gray45','gray30','gray15'), legend=c("E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic HSP","E.Coli"),)
graphics.off()


