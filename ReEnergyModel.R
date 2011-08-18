source("SQLShareLib_Wenjun.R")

sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))

energyCycle <- function(dataset, username=myUsername, contacts=fetchContacts(paste(dataset, "_surface_contacts.csv",sep=""), username), lambdaf=1, lambdau=1) {
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

  ener[[1]] <- freeEnergyModel(countMatrix, groDist["GroEL_Close", ], pfracs, lambdaf, lambdau)
  ener[[2]] <- freeEnergyModel(countMatrix, groDist["GroEL_Open", ], pfracs, lambdaf, lambdau)
  print(ener)

  return(ener)
}

freeEnergyModel <- function(countMatrix,yDist,pfracs,lambdaf,lambdau) {

  hydration <- array(0,ncol(countMatrix))
  contact <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    hydration[l] <- countMatrix["WATER",l] / sum(countMatrix[c(1:20,which(rownames(countMatrix) == "WATER")),l]) * (1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]) + countMatrix["FREE",l] / countMatrix["TOTAL",l]
    contact[l] <- 1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]
  }
  names(hydration) <- colnames(countMatrix)
  names(contact) <- colnames(countMatrix)

  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )  

  sigmaXln <- array(0,2)
  sigmaX <- array(0,2)
 
  for (i in 1:2) { 
    X <- array(0, ncol(pfracs[[i]]))
    for (j in 1:ncol(pfracs[[i]])) {
       Y <- rep(0, length(yDist))
       for (k in 1:length(yDist)) {
         t1 <- countMatrix[colnames(pfracs[[i]][j]), names(yDist[k])] * yDist[k] * contact[names(yDist)[k]]
         t2 <- hydration[colnames(pfracs[[i]][j])] / sum(hydration) * hydration[names(yDist)[k]] * yDist[k]
         Y[k] <- as.double(t1) + as.double(t2)
       }
       sigmaY <- log(sum(Y))
       X[j] <- median(pfracs[[i]][,j]) * sigmaY
     }
    sigmaXln[i] <- sum(X)
  }
  deltaG = -(sigmaXln[1] * lambdaf - sigmaXln[2] * lambdau)
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

minimizeEnergy <- function(dataset, username=myUsername, lambdaf=1,lambdau=1, contacts=fetchContacts(paste(dataset, "_surface_contacts.csv",sep=""), username=username)) {


  cutoff <- 0.3
  pfracsfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsfold <- fetchPDBIDs(dataset, username)

  cutoff <- -1.0
  pfracsunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
  pidsunfold <- fetchPDBIDs(dataset, username)

  deltaPs <- apply(pfracsfold, 2, median) * lambdaf - apply(pfracsunfold, 2, median) * lambdau
  countMatrix <- sampleContacts(contacts)

  hydration <- array(0,ncol(countMatrix))
  contact <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    hydration[l] <- countMatrix["WATER",l] / sum(countMatrix[c(1:20,which(rownames(countMatrix) == "WATER")),l]) * (1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]) + countMatrix["FREE",l] / countMatrix["TOTAL",l]
    contact[l] <- 1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]
  }
  names(hydration) <- colnames(countMatrix)
  names(contact) <- colnames(countMatrix)
  
  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )

  deltaPs <- deltaPs[match(colnames(countMatrix),names(deltaPs))]
  
  g <- function(dist) {
    dist <- dist[match(colnames(countMatrix),names(dist))]
    sum <- 0.
    for(i in 1:length(deltaPs)) {
      pxy <- 0.
      for(j in 1:length(dist)) {
        t1 <- countMatrix[i,j] * dist[j] * contact[j]
        t2 <- hydration[i] / sum(hydration) * hydration[j] * dist[j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sum <- sum + log(pxy) * deltaPs[i]
    }
    return(sum)
  }

#  print(g(groDist["GroEL_Close",]))
#  print(g(groDist["GroEL_Open",]))

  devI <- function(dist) {
    dist <- dist[match(colnames(countMatrix),names(dist))]
    devI <- array(0,length(dist))
    names(devI) <- names(dist)
    #print(countMatrix)
    #print(deltaPs)
    #print(dist)
    #print(water)
    for (i in 1:length(dist)) {
      sum <- 0.
      for (j in 1:length(deltaPs)) {
        pxy <- 0.
        for (k in 1:length(dist)) {
          t1 <- countMatrix[j,k] * dist[k] * contact[k]
          t2 <- hydration[j] / sum(hydration) * hydration[k] * dist[k]
          pxy <- pxy + as.double(t1) + as.double(t2)
        }
        d1 <- deltaPs[j] / pxy
        d2 <- countMatrix[j,i] * contact[i]
        d3 <- hydration[j] * hydration[i] / sum(hydration)
        sum <- sum + d1 * (d2 + d3)
      }
      devI[i] <- sum
    }
    return (devI)
  }

  ini <- c(0.1,0.05,0.01,0.01,0,0,0.2,0.11,0.01,0.01,0.02,0.1,0,0.03,0.04,0.07,0.06,0,0.042,0.016)
  names(ini) <- names(groDist["GroEL_Close",])
  
  opt <- optim(ini, g, devI, method = "CG", control = list(trace = T, fnscale = -1))
  print(opt[1])
  print(opt[2])
  print(opt[3])
  print(groDist["GroEL_Close",])

#  print(devI(groDist["GroEL_Close",]))

  numDev <- function(dist,delta) {
    dist <- dist[match(colnames(countMatrix),names(dist))]
    numDev <- array(0,length(dist))
    names(numDev) <- names(dist)
    #print(countMatrix)
    #print(deltaPs)
    #print(dist)
    #print(water)
    for (i in 1:length(dist)) {
      sum1 <- 0.
      sum2 <- 0.
      for(j in 1:length(deltaPs)) {
        pxy <- 0.
        for(k in 1:length(dist)) {
          if (k == i) {
            t1 <- countMatrix[j,k] * (dist[k] + delta) * contact[k]
            t2 <- hydration[j] / sum(hydration) * hydration[k] * (dist[k] + delta)
            pxy <- pxy + as.double(t1) + as.double(t2)}
          else {
            t1 <- countMatrix[j,k] * (dist[k]) * contact[k]
            t2 <- hydration[j] / sum(hydration) * hydration[k] * (dist[k])
            pxy <- pxy + as.double(t1) + as.double(t2)}
        }
      sum1 <- sum1 + log(pxy) * deltaPs[j]
      }

      for(j in 1:length(deltaPs)) {
        pxy <- 0.
        for(k in 1:length(dist)) {
          if (k == i) {
            t1 <- countMatrix[j,k] * (dist[k] - delta) * contact[k]
            t2 <- hydration[j] / sum(hydration) * hydration[k] * (dist[k] - delta)
            pxy <- pxy + as.double(t1) + as.double(t2)}
          else {
            t1 <- countMatrix[j,k] * (dist[k]) * contact[k]
            t2 <- hydration[j] / sum(hydration) * hydration[k] * (dist[k])
            pxy <- pxy + as.double(t1) + as.double(t2)}
        }
      sum2 <- sum2 + log(pxy) * deltaPs[j]
      }
      numDev[i] <- (sum1 - sum2) / (2 * delta)
    }
    return(numDev)
  }

#  print(numDev(groDist["GroEL_Close",], 0.00001))
}

#energyCycle("ecoli", "wenjunh")
#energyCycle("assist", username="wenjunh", contacts=fetchContacts("ecoli_surface_contacts.csv", "wenjunh"))
#energyBootstrap(1000, "ecoli", username="wenjunh")
minimizeEnergy("ecoli", "wenjunh")
