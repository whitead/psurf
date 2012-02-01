
source("SQLShareLib.R")

#get the interaction energies between residue types
#getInteractionEnergy <- function(dataset, username=NULL, glycine=FALSE) {
#
#  dataset <- paste(paste(dataset,"total","contacts", sep="_"), "csv" ,sep=".")
# 
#  #get the data
#  if(is.null(username)) {
#    countMatrix <- fetchContacts(dataset)
#  }
#  else{
#    countMatrix <- fetchContacts(dataset, username)
#  }
#
#  #re-order it
#  anames <- colnames(countMatrix[-c(1,2)])
#  aanum <- length(anames)
#  anames.ord <- order(anames)
#
#  countMatrix <- countMatrix[c(1,2,anames.ord + 2)]
#
#  #Including categories, i.e. aromatic, hydrophobic, polar, charged
#  aromatic <- c("PHE", "TYR", "TRP")
#  polar <- c("SER", "THR", "PRO", "ASN", "GLN", "GLY")
#  charged <- c("GLU", "LYS", "ARG", "ASP", "HIS")
#  hydrophobic <- c("ALA", "VAL", "LEU", "ILE", "MET")
#
#  
#  #sum it
#  contactMatrix <- sampleContacts(countMatrix, random=FALSE)
#  contactMatrix[1:aanum,] <- contactMatrix[order(rownames(contactMatrix)[1:aanum]),]
#  rownames(contactMatrix) <- c(sort(rownames(contactMatrix)[1:aanum]), rownames(contactMatrix)[-(1:aanum)])
#  contactMatrix[1:aanum, 1:aanum] <- contactMatrix[1:aanum,1:aanum] + t(contactMatrix[1:aanum, 1:aanum])
#
#
#
#  
#  #normalize it to the effects remove amounts of each amino acid
#  #add the effect of the free residues
#  for(i in 1:aanum) {
#    contactMatrix[1:aanum, i] <- contactMatrix[1:aanum, i] / sum(contactMatrix[1:aanum, i])
#    contactMatrix[1:aanum, i] <- contactMatrix[1:aanum, i] * (1 - contactMatrix["FREE", i] / contactMatrix["TOTAL", i])
#    contactMatrix["FREE", i] <- contactMatrix["FREE", i] / contactMatrix["TOTAL", i]
#  }
#
#  
#
#  #normalize it so all events sum to 1
#  normRows <- c(1:aanum, which(rownames(contactMatrix) == "FREE"))
#  contactMatrix[normRows, 1:aanum] <- contactMatrix[normRows, 1:aanum] / sum(contactMatrix[normRows, 1:aanum])
#  
#  #make it relative to being a free residue
#  mat <- matrix(rep(0, aanum**2), nrow=aanum)
#  rownames(mat) <- sort(anames)
#  colnames(mat) <- sort(anames)
#  for(i in 1:aanum) {
#    for(j in 1:aanum) {
#      mat[i,j] <- contactMatrix[i,j] * contactMatrix[j,i] / (contactMatrix["FREE",i] * contactMatrix["FREE", j])
#    }
#  }
#  
#  mat <- -log(mat)
#
#  #make glycine 0, if wanted
#  if(!glycine) {
#    mat["GLY", ] <- 0
#    mat[,"GLY"] <- 0
#  }
#  
#  return(mat)
#}
#
#
#
#
#sql <- paste("select * FROM [whitead@washington.edu].[GroEL_counts.csv]")
#rawData <- fetchdata(sql)
#
#
#groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
#groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)
#
#groDist[1, ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
#groDist[2, ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))
#
#b <- array(0,20)
#
#cairo_pdf('groDist_SE.pdf', width=3.42, height=2.58, pointsize=8)
#par(family="LMSans10", cex.axis=0.6)
#barx <- barplot(as.matrix(groDist), col=c("gray15","gray75"), main="", xlab="Amino Acid", ylab="", beside=T, names.arg=aalist.sh, ylim = c(0.00,0.20))
#legend("topright", col=c("gray15","gray75"), cex=0.8, legend=c("GroEL-GroES (Closed)","GroEL (Open)"), pch=15)
#graphics.off()


#plot2 The plot of 5 species

hspRawCount <- fetchHSPResCount("hsp","in",username="whitead")

for (i in 1:nrow(hspRawCount)) {
  hspRawCount[i,] <- hspRawCount[i,] / sum(hspRawCount[i,])
}
hspRawCount <- as.matrix(hspRawCount)

hspDist <- matrix(0,nrow(hspRawCount)+1, ncol(hspRawCount))
rownames(hspDist) <- c("E.Coli","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Yeast CCT")
colnames(hspDist) <- colnames(hspRawCount)

hspDist[2,] <- hspRawCount["1SX4",]
hspDist[3,] <- hspRawCount["1WE3",]
hspDist[4,] <- hspRawCount["3KFB",]
hspDist[5,] <- hspRawCount["2CG9",]
hspDist[6,] <- hspRawCount["3P9D",]

#various statistics
print(sum(hspRawCount["1WE3", c("ASP", "GLU", "LYS", "ARG", "HIS")]))
print(sum(hspRawCount["1SX4", c("ASP", "GLU", "LYS", "ARG", "HIS")]))

cutoff <- 0.3
psurf <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE, "wenjunh")

ResSum <- rep(0, nrow(psurf))

for (i in 1:nrow(psurf)) {
  ResSum[i] <- sum(psurf[i,c("LYS","ARG","GLU","ASP","HIS")])
}
quantile(ResSum)

hspDist[1,] <- apply(psurf, MARGIN=2, FUN=mean)

up <- rep(0,ncol(psurf))
low <- rep(0,ncol(psurf))

for (i in ncol(psurf)) {
  up[i] <- quantile(psurf[,i],0.95)
  low[i] <- quantile(psurf[,i],0.05)
}

cairo_pdf('HSP_Protein_Dist.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
barplot(hspDist, beside=TRUE, col=c('blue','gray75','gray60','gray45','gray30','gray15'), ylim=c(0.0,0.3), names.arg=aalist.sh,)
legend("topright", col=c('blue','gray75','gray60','gray45','gray30','gray15'), legend=c("E.Coli","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic CCT"), pch=rep(15,6), cex=0.8)
graphics.off()
