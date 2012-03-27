
source("SQLShareLib.R")

#get the interaction energies between residue types
getInteractionEnergy <- function(username=NULL, glycine=FALSE) {

  dataset <- paste(paste("h2","ionic","sufrace", "contacts", sep="_"), "csv" ,sep=".")
 
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
  contactMatrix <- sampleContacts(countMatrix, random=TRUE)

  
  contactMatrix[1:aanum,] <- contactMatrix[order(rownames(contactMatrix)[1:aanum]),]
  rownames(contactMatrix) <- c(sort(rownames(contactMatrix)[1:aanum]), rownames(contactMatrix)[-(1:aanum)])

  #Adjust sums, so that pairs which were double counted or multiple pairings are fixed
  for(i in 1:aanum) {
    contactMatrix[1:aanum,i] <- (contactMatrix[1:aanum,i] + 1) * (contactMatrix["TOTAL", i] - contactMatrix["FREE", i]) / (sum(contactMatrix[1:aanum,i]) + aanum)
  }

  #Make it symmetric
  contactMatrix[1:aanum, 1:aanum] <- (contactMatrix[1:aanum,1:aanum] + t(contactMatrix[1:aanum, 1:aanum])) / 2

  #normalize it so all events sum to 1
  normRows <- c(1:aanum, which(rownames(contactMatrix) == "FREE"))
  contactMatrix[normRows, 1:aanum] <- contactMatrix[normRows, 1:aanum] / sum(contactMatrix[normRows, 1:aanum])




  mat <- matrix(0, nrow=aanum, ncol=aanum)

  rownames(mat) <- sort(anames)
  colnames(mat) <- sort(anames)
  for(i in 1:aanum) {
    for(j in 1:aanum) {
      mat[i,j] <- contactMatrix[i,j]  / (sum(contactMatrix[1:aanum,i]) * sum( contactMatrix[1:aanum, j]))
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

kt <- 2.494
bootstrap <- 3
pairs <- empty.df(cnames=c("RE", "RD", "KE", "KD"), rnames=1:bootstrap)
names(pairs) <- c("Arg-Glu", "Arg-Asp", "Lys-Glu", "Lys-Asp")
for(i in 1:bootstrap){
  cat(paste("Bootrap:", i, "\n"))
  intMat <- getInteractionEnergy(username="whitead")
  pairs[i,] <- kt*c(intMat["ARG", "GLU"], intMat["ARG", "ASP"], intMat["LYS", "GLU"], intMat["LYS", "ASP"])
}

ims <- apply(pairs, MARGIN=2, FUN=median)
ilower <- ims - apply(pairs, MARGIN=2, FUN=function(x) {quantile(x, c(0.025))} )
iupper <- apply(pairs, MARGIN=2, FUN=function(x) {quantile(x, c(0.975))} ) - ims

ims <- ims - ims[1]
ilower <- ilower
iupper <- iupper

color <- c("blue", "red")


#Free energies from simulation techniques
deltaG <- c(1.7,7.7, 3.6, 4.2)
deltaG <- deltaG - deltaG[1]
deltaG.upper <- c(0.4, 0.3, 0.6, 0.5)
deltaG.lower <- c(-0.4, -0.3, -0.6, -0.5)


iyy <- matrix(c(ims, deltaG), byrow=T, nrow=2)
ill <- matrix(c(ilower, deltaG.lower), byrow=T, nrow=2)
iuu <- matrix(c(iupper, deltaG.upper), byrow=T, nrow=2)

cairo_pdf("salt_interaction.pdf", width=4.5, height=3.3, pointsize=12)
par(family="LMSans10", cex.axis=0.8)
barx <- barplot(iyy, xlab="Pair Type", axis.lty=1, ylab=expression(paste(Delta * G, " [kJ/mol]")), names.arg=c("Arg-Glu", "Arg-Asp", "Lys-Glu", "Lys-Asp"), beside=T, col=color, ylim=c(min(c(iyy - ill - 0.06)), max(c(iyy + iuu + 0.06))))
error.bar(barx, iyy, lower=ill, upper=iuu, length=0.05)
legend("topright", c("Bioinformatics", "Metadynamics"), pch=15, ncol=1, col=color)
graphics.off()

print(ims - ilower)
print(ims + iupper)
print(ims)
