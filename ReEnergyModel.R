source("SQLShareLib.R")

groDist <- fetchFrame("greoldist.csv", username="wenjunh")

energyCycle <- function(protein, distribution) {
  countMatrix <- sampleCounts(distribution)

  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues(protein, cutoff, normalize=TRUE)
  pidsfold <- fetchPDBIDs(protein)

  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues(protein, cutoff, normalize=TRUE)
  pidsunfold <- fetchPDBIDs(protein)


  pfracs <- list(pfracsfold, pfracsunfold)

  Gclose <- 0
  Gopen <- 0
  ener <- c(Gclose,Gopen)
  names(ener) <- c("Gclose","Gopen")

  ener["Gclose"] <- freeEnergyModel(countMatrix, groDist["Close", ], pfracs)

  return(ener)
}

freeEnergyModel <- function(countMatrix,yDist,pfracs) {

  Num <- distriNum

  rownames(countMatrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","WATER")
  colnames(countMatrix) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

  countMatrix <- countMatrix[order(rownames(countMatrix)),order(colnames(countMatrix))]
  water <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    water[l] <- countMatrix["WATER",l] / sum(countMatrix["WATER",])
  }
  names(water) <- colnames(countMatrix)
#  print(water)
#  print(countMatrix)
  
  deltaGOpen <- array(0,1000)
  deltaGClose <- array(0,1000)


  sfracGrEL <- list(sfracGrELClose_1C, sfracGrELOpen_1C)

  sigmaXln <- array(0,2)
  sigmaX <- array(0,2)
 
  for (i in 1:2) { 
    X <- array(0, ncol(pfracs[[i]]))
    for (j in 1:ncol(pfracs[[i]])) {
       Y <- array(0, length(yDist))
       for (k in 1:length(yDist)) {
         Y[k] <- countMatrix[which(rownames(countMatrix) == colnames(pfracs[[i]][j])), which(colnames(countMatrix) == names(yDist[k]))] * yDist[k] * (1-countMatrix[which(rownames(countMatrix) == "WATER"), which(colnames(countMatrix) == names(yDist[k]))]) + water[which(names(water) == colnames(pfracs[[i]][j]))] * countMatrix[which(rownames(countMatrix) == "WATER"), which(colnames(countMatrix) == names(yDist[k]))] * yDist[k] 
       }
#       print(Y)
       sigmaY <- log(sum(Y), base = exp(1))
#       print(sigmaY)
       X[j] <- mean(pfracs[[i]][,j]) * sigmaY
     }
    sigmaXln[i] <- sum(X)
    sigmaX[i] <- exp(sum(X))
  }
  deltaG <- array(0,2)
  names(deltaG) <- c("deltaGln","deltaG")
  deltaG[1] <- sigmaXln[1]-sigmaXln[2]
  deltaG[2] <- sigmaX[1]-sigmaX[2]
  return (deltaG)

}

energybootstrap <- function(bootstrap,protein,distribution) {
  bootstrap <- bootstrap
  source('obtainSurfaceResidues.R')
  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues(protein, cutoff, normalize=TRUE)
  pidsfold <- fetchPDBIDs(protein)

  colnames(pfracsfold) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  pfracsfold <- pfracsfold[,order(colnames(pfracsfold))]
  
  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues(protein, cutoff, normalize=TRUE)
  pidsunfold <- fetchPDBIDs(protein)

  colnames(pfracsunfold) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  pfracsunfold <- pfracsunfold[,order(colnames(pfracsunfold))]

  pfracs <- list(pfracsfold, pfracsunfold)

  source('countMatrixGenerator.R')
  Gclose <- matrix(0,bootstrap,2)
  colnames(Gclose) <- c("G_ln","G")
  Gopen <- matrix(0,bootstrap,2)
  colnames(Gopen) <- c("G_ln","G")
  ener <- list(Gclose,Gopen)
  names(ener) <- c("Gclose","Gopen")
  for (i in 1:bootstrap) {    
    countMatrix <- matrixGenerator(distribution)
    for (distriNum in 1:2) {
      ener[[distriNum]][i,] <- freeEnergyModel(countMatrix, distriNum, pfracs)
    }
    cat(paste(i,"/",bootstrap))
  }
  return(ener)
}




