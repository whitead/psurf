load('sfracGrELClose_1C.Rdata')
load('sfracGrELOpen_1C.Rdata')

#load('allcount.Rdata')
#countMatrix <- count
#countMatrix <- t(countMatrix)

freeEnergyModel <- function(countMatrix,distriNum) {
  InterfaceCount <- countMatrix
  Num <- distriNum
#  read in counting matrix
#  InterfaceCount <- read.table('H1_Counts_table_ppairs.txt',header=T)
#  InterfaceCount <- t(InterfaceCount)
  rownames(InterfaceCount) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","FREE")
  colnames(InterfaceCount) <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

  InterfaceCount <- InterfaceCount[order(rownames(InterfaceCount)),order(colnames(InterfaceCount))]
  source('countSurfaceResidues.R')

  deltaGOpen <- array(0,1000)
  deltaGClose <- array(0,1000)

#  obtain surface residues for both E.Coli fold and E.Coli unfold
  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE)
  pidsfold <- fetchPDBIDs("ecoli")

  colnames(pfracsfold) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  pfracsfold <- pfracsfold[,order(colnames(pfracsfold))]
  
  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues("ecoli", cutoff, normalize=TRUE)
  pidsunfold <- fetchPDBIDs("ecoli")

  colnames(pfracsunfold) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  pfracsunfold <- pfracsunfold[,order(colnames(pfracsunfold))]

  pfracs <- list(pfracsfold, pfracsunfold)
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
       #print(Y)
       sigmaY <- log(sum(Y), base = exp(1))
       #print(sigmaY)
       X[j] <- mean(pfracs[[i]][,j]) * sigmaY
    }
    #print(X)
    sigmaXln[i] <- sum(X)
    #print(sigmaXln[i])
    sigmaX[i] <- exp(sum(X))
  }
  deltaG <- array(0,2)
  names(deltaG) <- c("deltaGln","deltaG")
  deltaG[1] <- sigmaXln[1]-sigmaXln[2]
  deltaG[2] <- sigmaX[1]-sigmaX[2]
  #print(deltaG)
  return (deltaG)

}

#countMatrix <- read.table('H1_Counts_table_ppairs.txt',header=T)
#load('normCys.Rdata')
#countMatrix <- countnormCys
#load('normCysFree.Rdata')
#countMatrix <- countnormCysFree
#countMatrix <- t(countMatrix)
load('allcount.Rdata')
countMatrix <- count


Gclose <- 0
Gopen <- 0
ener <- list(Gclose,Gopen)
names(ener) <- c("Gclose","Gopen")

for (distriNum in 1:2) {
  ener[[distriNum]] <- freeEnergyModel(countMatrix, distriNum)
}
print(ener)
