load('sfracGrELClose_1C.Rdata')
load('sfracGrELOpen_1C.Rdata')

#load('allcount.Rdata')

energyCycle <- function(counts,protein) {
  countMatrix <- counts
  
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

  Gclose <- 0
  Gopen <- 0
  ener <- list(Gclose,Gopen)
  names(ener) <- c("Gclose","Gopen")

  for (distriNum in 1:2) {
    ener[[distriNum]] <- freeEnergyModel(countMatrix, distriNum, pfracs)
  }
  return(ener)
}

freeEnergyModel <- function(countMatrix,distriNum,pfracs) {
  InterfaceCount <- countMatrix
  Num <- distriNum
  pfracs <- pfracs

  rownames(InterfaceCount) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","FREE")
  colnames(InterfaceCount) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

  InterfaceCount <- InterfaceCount[order(rownames(InterfaceCount)),order(colnames(InterfaceCount))]
 

  deltaGOpen <- array(0,1000)
  deltaGClose <- array(0,1000)


  sfracGrEL <- list(sfracGrELClose_1C, sfracGrELOpen_1C)

  sigmaXln <- array(0,2)
  sigmaX <- array(0,2)
 
  for (i in 1:2) { 
    X <- array(0, ncol(pfracs[[i]]))
    for (j in 1:ncol(pfracs[[i]])) {
       Y <- array(0, length(sfracGrEL[[Num]]))
       for (k in 1:length(sfracGrEL[[Num]])) {
         Y[k] <- InterfaceCount[which(rownames(InterfaceCount) == colnames(pfracs[[i]][j])), which(colnames(InterfaceCount) == names(sfracGrEL[[Num]][k]))] * sfracGrEL[[Num]][k]
       }
       sigmaY <- log(sum(Y), base = exp(1))
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

energybootstrap <- function(bootstrap,protein) {
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
    countMatrix <- matrixGenerator('surf.cal')
    for (distriNum in 1:2) {
      ener[[distriNum]][i,] <- freeEnergyModel(countMatrix, distriNum, pfracs)
    }
    cat(paste(i,"/",bootstrap))
  }
  return(ener)
}




