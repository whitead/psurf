source("SQLShareLib.R")

#get the interaction energies between residue types
getInteractionEnergy <- function(dataset, username=NULL, glycine=FALSE) {

 dataset <- paste(paste(dataset,"backbone","contacts", sep="_"), "csv" ,sep=".")

  #get the data
  if(is.null(username)) {
    countMatrix <- fetchContacts(dataset)
  }
  else{
    countMatrix <- fetchContacts(dataset, username)
  }

  #re-order it
  anames <- colnames(countMatrix[-c(1,2)])
  aanum <- length(anames)
  anames.ord <- order(anames)

  countMatrix <- countMatrix[c(1,2,anames.ord + 2)]

  #sum it
  contactMatrix <- sampleContacts(countMatrix, random=FALSE)
  contactMatrix[1:aanum,] <- contactMatrix[order(rownames(contactMatrix)[1:aanum]),]
  rownames(contactMatrix) <- c(sort(rownames(contactMatrix)[1:aanum]), rownames(contactMatrix)[-(1:aanum)])

  
  #Adjust sums, so that pairs which were double counted or multiple pairings are fixed
  for(i in 1:aanum) {
    contactMatrix[1:aanum,i] <- contactMatrix[1:aanum,i] * (contactMatrix["TOTAL", i] - contactMatrix["FREE", i]) / sum(contactMatrix[1:aanum,i])
  }

  #Make it symmetric
  contactMatrix[1:aanum, 1:aanum] <- (contactMatrix[1:aanum,1:aanum] + t(contactMatrix[1:aanum, 1:aanum])) / 2


  
  #normalize it so all events sum to 1
  normRows <- c(1:aanum, which(rownames(contactMatrix) == "FREE"))
  contactMatrix[normRows, 1:aanum] <- contactMatrix[normRows, 1:aanum] / sum(contactMatrix[normRows, 1:aanum])




  #Get the categories to consider 
  aromatic <- c("PHE", "TYR", "TRP")
  polar <- c("SER", "THR", "PRO", "ASN", "GLN", "GLY")
  charged <- c("GLU", "LYS", "ARG", "ASP", "HIS")
  hydrophobic <- c("ALA", "VAL", "LEU", "ILE", "MET")
  categories <- list(Aromatic=aromatic, Polar=polar, Charged=charged, Hydrophobic=hydrophobic)

  #make it relative to being a free residue

  
  mat <- matrix(rep(0, aanum * (aanum + length(categories))), nrow=aanum + length(categories))

  rownames(mat) <- c(sort(anames), names(categories))
  colnames(mat) <- sort(anames)
  for(i in 1:aanum) {
    for(j in 1:aanum) {
      mat[i,j] <- contactMatrix[i,j]  / (sum(contactMatrix[1:aanum,i]) * sum( contactMatrix[1:aanum, j]))
    }
  }

  for(i in 1:length(categories)) {

    for(j in 1:aanum) { 
      mat[names(categories)[i] ,j] <- sum(contactMatrix[j,categories[[i]]]) / (sum(sapply(categories[[i]], function(x) {contactMatrix[1:aanum,x]})) * sum( contactMatrix[1:aanum, j]))
      
    }
    
  }

  mat <- -log(mat)

  #make glycine 0, if wanted
  if(!glycine) {
    mat["GLY", ] <- 0
    mat[,"GLY"] <- 0
  }
  
  return(mat)
}



#Interaction Table
chis <- getInteractionEnergy("ecoli", "wenjunh")

interactionTable <- chis

surface <- fetchAllSurfResidues("ecoli", cutoff, normalize=FALSE, "wenjunh")
interior <- fetchAllBuriedResidues("ecoli", cutoff, normalize=FALSE, "wenjunh")
all <- fetchAllSurfResidues("ecoli", -1, normalize=FALSE, "wenjunh")

distributions <- list(Proteins=all, Surface=surface, Interior=interior, Diff=(surface - all))

for(i in 1:length(distributions)) {

  dist <- distributions[[i]]
  energies <- matrix(0, nrow=nrow(dist), ncol=20)

  for(j in 1:nrow(dist)) {

    energies[j,] <- as.matrix(dist[j,]) %*% chis[1:20, 1:20]
    
  }

  interactionTable <- rbind(interactionTable, apply(energies, 2, median))
  
}

rownames(interactionTable) <- c(rownames(chis), "Proteins", "Surface", "Interior", "Surface - Proteins")

#remove glycine
gindex <- which(colnames(interactionTable) == "GLY")

interactionTable <- interactionTable[-gindex, -gindex]

write.table(round(interactionTable, 2), file="Interaction_Table.txt", quote=FALSE, sep=" & ", eol="\\\\\n")


sql <- paste("select * FROM [whitead@washington.edu].[GroEL_counts.csv]")
rawData <- fetchdata(sql)



groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist[1, ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroL_Close", ]))
groDist[2, ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroE_Open", ]))


b <- array(0,20)

cairo_pdf('groDist_SE.pdf', width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.6)

barx <- barplot(as.matrix(groDist), col=c("gray15","gray75"), main="", xlab="Amino Acid", ylab="", beside=T, names.arg=aalist.sh, ylim = c(0.00,0.20))
legend("topright", col=c("gray15","gray75"), cex=0.8, legend=c("GroEL-GroES (Closed)","GroEL (Open)"), pch=15)
graphics.off()


#plot2 The plot of 5 species

hspRawCount <- fetchHSPResCount("hsp","in",username="whitead")

for (i in 1:nrow(hspRawCount)) {
  hspRawCount[i,] <- hspRawCount[i,] / sum(hspRawCount[i,])
}
hspRawCount <- as.matrix(hspRawCount)

hspDist <- matrix(0,nrow(hspRawCount)+1, ncol(hspRawCount))
rownames(hspDist) <- c("Thermo Thermopilius","E.Coli GroEL","Thermo GroEL","Group II","HSP90","Yeast CCT")
colnames(hspDist) <- colnames(hspRawCount)

hspDist[2,] <- hspRawCount["1SX4",]
hspDist[3,] <- hspRawCount["1WE3",]
hspDist[4,] <- hspRawCount["3KFB",]
hspDist[5,] <- hspRawCount["2CG9",]
hspDist[6,] <- hspRawCount["3P9D",]

#various statistics
print(sum(hspRawCount["1WE3", c("ASP", "GLU", "LYS", "ARG", "HIS")]))
print(sum(hspRawCount["1SX4", c("ASP", "GLU", "LYS", "ARG", "HIS")]))
print(sum(hspRawCount["3KFB", c("ASP", "GLU")]))
print(sum(hspRawCount["3KFB", c("LYS", "HIS", "ARG")]))

cutoff <- 0.3
psurf <- fetchAllSurfResidues("thermus_nogaps", cutoff, normalize=TRUE, "wenjunh")

chargeSum <- rep(0, nrow(psurf))

for (i in 1:nrow(psurf)) {
  chargeSum[i] <- sum(psurf[i,c("LYS","ARG","GLU","ASP","HIS")])
}
quantile(chargeSum, seq(0,1,0.01))

hspDist[1,] <- apply(psurf, MARGIN=2, FUN=median)

up <- matrix(0,6,ncol(psurf))
low <- matrix(0,6,ncol(psurf))

for (i in 1:ncol(psurf)) {
  up[1,i] <- quantile(psurf[,i],0.95) - median(psurf[,i])
  low[1,i] <- median(psurf[,i]) - quantile(psurf[,i],0.05)
}

cairo_pdf('HSP_Protein_Dist_Thermo.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
barx <- barplot(hspDist, beside=TRUE, col=c('blue','gray75','gray60','gray45','gray30','gray15'), ylim=c(0.0,0.3), names.arg=aalist.sh,)
error.bar(barx,hspDist,lower=low,upper=up,length=0.01)
legend("topright", col=c('blue','gray75','gray60','gray45','gray30','gray15'), legend=c("Thermo Thermopilius","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic CCT"), pch=rep(15,6), cex=0.8)
graphics.off()

#plot3 The plot of E.Coli, Human, and Thermophilius
cutoff <- 0.3
psurf_ecoli <- fetchAllSurfResidues("ecoli_nogaps", cutoff, normalize=TRUE, "wenjunh")
psurf_thermus <- fetchAllSurfResidues("thermus_nogaps", cutoff, normalize=TRUE, "wenjunh")
psurf_human <- fetchAllSurfResidues("h2_w_nogaps", cutoff, normalize=TRUE, "whitead")

protein_dist <- matrix(0, 3, 20)
rownames(protein_dist) <- c("E.Coli","Thermus Thermopilius","Human")
colnames(protein_dist) <- colnames(psurf_ecoli)

protein_dist[1,] <- apply(psurf_ecoli, MARGIN=2, FUN=median)
protein_dist[2,] <- apply(psurf_thermus, MARGIN=2, FUN=median)
protein_dist[3,] <- apply(psurf_human, MARGIN=2, FUN=median)

up <- matrix(0,3,ncol(psurf_ecoli))
low <- matrix(0,3,ncol(psurf_ecoli))

for (i in 1:ncol(psurf_ecoli)) {
  up[1,i] <- quantile(psurf_ecoli[,i],0.95) - median(psurf_ecoli[,i])
  low[1,i] <- median(psurf_ecoli[,i]) - quantile(psurf_ecoli[,i],0.05)
}

for (i in 1:ncol(psurf_thermus)) {
  up[2,i] <- quantile(psurf_thermus[,i],0.95) - median(psurf_thermus[,i])
  low[2,i] <- median(psurf_thermus[,i]) - quantile(psurf_thermus[,i],0.05)
}

for (i in 1:ncol(psurf_human)) {
  up[3,i] <- quantile(psurf_human[,i],0.95) - median(psurf_human[,i])
  low[3,i] <- median(psurf_human[,i]) - quantile(psurf_human[,i],0.05)
}

cairo_pdf('Three_Protein_SurfRes.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
barx <- barplot(protein_dist, beside=TRUE, col=c('gray75','gray50','gray25'), ylim=c(0.0,0.3), names.arg=aalist.sh,)
error.bar(barx,protein_dist,lower=low,upper=up,length=0.01)
legend("topright", col=c('gray75','gray50','gray25'), legend=c("E.Coli","Thermus Thermopilius","Human"), pch=rep(15,6), cex=0.8)
graphics.off()
