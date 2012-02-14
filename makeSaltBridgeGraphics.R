
source("SQLShareLib.R")

#get the interaction energies between residue types
getInteractionEnergy <- function(username=NULL, glycine=FALSE) {

  dataset <- paste(paste("h2","ionic","contacts", sep="_"), "csv" ,sep=".")
 
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
bootstrap <- 1000
pairs <- empty.df(cnames=c("RE", "RD", "KE", "KD"), rnames=1:bootstrap)
names(pairs) <- c("Arg-Glu", "Arg-Asp", "Lys-Glu", "Lys-Asp")
for(i in 1:bootstrap){
  print(paste("Bootrap:", i))
  intMat <- getInteractionEnergy(username="whitead")
  pairs[i,] <- kt*c(intMat["ARG", "GLU"], intMat["ARG", "ASP"], intMat["LYS", "GLU"], intMat["LYS", "ASP"])
}

ims <- apply(pairs, MARGIN=2, FUN=median)
ilower <- ims - apply(pairs, MARGIN=2, FUN=function(x) {quantile(x, c(0.025))} )
iupper <- apply(pairs, MARGIN=2, FUN=function(x) {quantile(x, c(0.975))} ) - ims

color <- c("dark gray")

cairo_pdf("salt_interaction.pdf", width=4.5, height=3.3, pointsize=12)
par(family="LMSans10", cex.axis=0.8)
barx <- barplot(ims, xlab="Pair Type", space=0.5, ylab=expression(paste(Delta * G, " [kJ/mol]")), names.arg=c("Arg-Glu", "Arg-Asp", "Lys-Glu", "Lys-Asp"), ylim=c(min(c(ims - ilower - 0.06)), 0), col=color)
error.bar(barx, ims, lower=ilower, upper=iupper, length=0.05)
graphics.off()

print(ims - ilower)
print(ims + iupper)
