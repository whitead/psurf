source("SQLShareLib.R")

sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))

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

  Gclose <- c(0,0)
  Gopen <- c(0,0)
  ener <- list(Gclose,Gopen)
  names(ener) <- list("Gclose","Gopen")

  ener[[1]] <- freeEnergyModel(countMatrix, groDist["GroEL_Close", ], pfracs)
  ener[[2]] <- freeEnergyModel(countMatrix, groDist["GroEL_Open", ], pfracs)

  return(ener)
}

freeEnergyModel <- function(countMatrix,yDist,pfracs) {

  water <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    water[l] <- countMatrix["WATER",l] / sum(countMatrix["WATER",])
  }
  names(water) <- colnames(countMatrix)

  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )  

  deltaGOpen <- array(0,1000)
  deltaGClose <- array(0,1000)

  sigmaXln <- array(0,2)
  sigmaX <- array(0,2)
 
  for (i in 1:2) { 
    X <- array(0, ncol(pfracs[[i]]))
    for (j in 1:ncol(pfracs[[i]])) {
       Y <- rep(0, length(yDist))
       for (k in 1:length(yDist)) {
         t1 <- countMatrix[colnames(pfracs[[i]][j]), names(yDist[k])] * yDist[k] * (1 - countMatrix["WATER", names(yDist[k])])
         t2 <- water[colnames(pfracs[[i]][k])] * countMatrix["WATER", names(yDist[j])] * yDist[k]
         Y[k] <- as.double(t1) + as.double(t2)
       }
       sigmaY <- log(sum(Y))
       X[j] <- median(pfracs[[i]][,j]) * sigmaY
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

energyBootstrap <- function(bootstrap,dataset,username=myUsername,
                            contacts=fetchContacts(paste(dataset, "_surface_contacts.csv",sep=""), username=username)) {

  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsfold <- fetchPDBIDs(dataset, username)

  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsunfold <- fetchPDBIDs(dataset, username)


  pfracs <- list(pfracsfold, pfracsunfold)

  Gclose <- matrix(0,nrow=bootstrap,ncol=2)
  colnames(Gclose) <- c("deltaGln","deltaG")
  Gopen <- matrix(0,nrow=bootstrap,ncol=2)
  colnames(Gopen) <- c("deltaGln","deltaG")
  ener <- list(Gclose,Gopen)
  names(ener) <- c("Gclose","Gopen")
  for (i in 1:bootstrap) {    
    countMatrix <- sampleContacts(contacts)
    for (distNum in 1:2) {
      ener[[distNum]][i,] <- freeEnergyModel(countMatrix, groDist[distNum,], pfracs)
    }
    cat(paste(i,"/",bootstrap))
  }
  return(ener)
}

minimizeEnergy <- function(dataset, username=myUsername,  contacts=fetchContacts(paste(dataset, "_surface_contacts.csv",sep=""), username=username)) {


  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsfold <- fetchPDBIDs(dataset, username)

  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsunfold <- fetchPDBIDs(dataset, username)

  deltaPs <- apply(pfracsfold, 2, median) - apply(pfracsunfold, 2, median)
  countMatrix <- sampleContacts(contacts)

  water <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    water[l] <- countMatrix["WATER",l] / sum(countMatrix["WATER",])
  }
  names(water) <- colnames(countMatrix)

  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )

  deltaPs <- deltaPs[match(colnames(countMatrix),names(deltaPs))]
  
  g <- function(dist) {
    dist <- dist[match(colnames(countMatrix),names(dist))]
    sum <- 0.
    for(i in 1:length(deltaPs)) {
      pxy <- 0.
      for(j in 1:length(dist)) {
        t1 <- countMatrix[i,j] * dist[j] * (1 - countMatrix["WATER", j])
        t2 <- water[j] * countMatrix["WATER", i] * dist[j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sum <- sum + log(pxy) * deltaPs[i]
    }
    return(sum)
  }

  print(g(groDist["GroEL_Close",]))
  print(g(groDist["GroEL_Open",]))

  devI <- function(dist) {
    dist <- dist[match(colnames(countMatrix),names(dist))]
    devI <- array(0,length(dist))
    names(devI) <- names(dist)
    #print(countMatrix)
    #print(deltaPs)
    #print(dist)
    #print(water)
    for (i in 1:length(deltaPs)) {
      for (j in 1:length(dist)) {
        pxy <- 0.
        for (k in 1:length(dist)) {
          t1 <- countMatrix[i,k] * dist[k] * (1 - countMatrix["WATER", k])
          t2 <- water[i] * countMatrix["WATER", k] * dist[k]
          pxy <- pxy + as.double(t1) + as.double(t2)
        }
        d1 <- deltaPs[i] / pxy
        d2 <- countMatrix[i,j] * (1 - countMatrix["WATER", j])
        d3 <- water[i] * countMatrix["WATER", j]
      }
      devI[i] <- d1 * (d2 + d3)
    }
    return (devI)
  }
  #print(devI(groDist["GroEL_Close",]))
}
energyCycle("ecoli", "wenjunh")
#energyBootstrap(1000, "ecoli", username="wenjunh")
minimizeEnergy("ecoli", "wenjunh")
