#!/usr/bin/env Rscript
source("SQLShareLib.R")

#get the interaction energies between residue types
getInteractionEnergy <- function(dataset, username=FALSE, random=FALSE, glycine=FALSE, countMatrix=fetchContacts(paste(paste(dataset,"backbone","contacts", sep="_"), "csv" ,sep="."), username)) {

  #re-order it
  anames <- colnames(countMatrix[-c(1,2)])
  aanum <- length(anames)
  anames.ord <- order(anames)

  countMatrix <- countMatrix[c(1,2,anames.ord + 2)]

  #sum it
  contactMatrix <- sampleContacts(countMatrix, random=random)
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
