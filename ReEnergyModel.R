source("SQLShareLib.R")

sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.integer(groDist["GroEL_Close", ]) / sum(as.integer(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.integer(groDist["GroEL_Open", ]) / sum(as.integer(groDist["GroEL_Open", ]))

energyCycle <- function(dataset, username=myUsername, contacts=fetchContacts(paste(dataset, "_surface_contacts.csv",sep=""), username)) {
  countMatrix <- sampleContacts(contacts)

  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsfold <- fetchPDBIDs(dataset, username)

  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsunfold <- fetchPDBIDs(dataset, username)


  pfracs <- list(pfracsfold, pfracsunfold)

  Gclose <- 0
  Gopen <- 0
  ener <- c(Gclose,Gopen)
  names(ener) <- c("Gclose","Gopen")

  ener[1] <- freeEnergyModel(countMatrix, groDist["GroEL_Close", ], pfracs)
  ener[2] <- freeEnergyModel(countMatrix, groDist["GroEL_Open", ], pfracs)

  return(ener)
}

freeEnergyModel <- function(countMatrix,yDist,pfracs,bootstrap=1) {

  water <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    water[l] <- countMatrix["WATER",l] / sum(countMatrix["WATER",])
  }
  names(water) <- colnames(countMatrix)

  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )
  
  print(water)
  print(countMatrix)
  
  deltaGOpen <- array(0,1000)
  deltaGClose <- array(0,1000)

  sigmaXln <- array(0,2)
  sigmaX <- array(0,2)
 
  for (i in 1:2) { 
    X <- array(0, ncol(pfracs[[i]]))
    for (j in 1:ncol(pfracs[[i]])) {
       Y <- array(0, length(yDist))
       print(Y)
       for (k in 1:length(yDist)) {
         print(Y[k])
         Y[k] <- countMatrix[colnames(pfracs[[i]][j]), names(yDist[k])] * yDist[k] * (1 - countMatrix["WATER", names(yDist[k])]) +
           water[colnames(pfracs[[i]][j])] * countMatrix["WATER", names(yDist[k])] * yDist[k]
         print(Y[k])
       }
       print(Y)
       sigmaY <- log(sum(Y))
       print(sigmaY)
       X[j] <- mean(pfracs[[i]][,j]) * sigmaY
     }
    sigmaXln[i] <- sum(X)
    sigmaX[i] <- sum(exp(X))
  }
  deltaG <- array(0,2)
  names(deltaG) <- c("deltaGln","deltaG")
  deltaG[1] <- sigmaXln[1]-sigmaXln[2]
  deltaG[2] <- sigmaX[1]/sigmaX[2]
  return (deltaG)

}

energybootstrap <- function(bootstrap,dataset,contacts) {

  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE)
  pidsfold <- fetchPDBIDs(dataset)

  colnames(pfracsfold) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  pfracsfold <- pfracsfold[,order(colnames(pfracsfold))]
  
  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE)
  pidsunfold <- fetchPDBIDs(dataset)

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
    countMatrix <- matrixGenerator(contacts)
    for (distriNum in 1:2) {
      ener[[distriNum]][i,] <- freeEnergyModel(countMatrix, distriNum, pfracs)
    }
    cat(paste(i,"/",bootstrap))
  }
  return(ener)
}




energyCycle("ecoli", "wenjunh")
