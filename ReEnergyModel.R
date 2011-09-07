library(alabama)
library(KernSmooth)
source("SQLShareLib.R")

sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))



energyCycle <- function(dataset, username=myUsername, contacts=fetchContacts(paste(dataset, "_surface_contacts.csv",sep=""), username), lambdaf=1, lambdau=1, split=FALSE) {
  countMatrix <- sampleContacts(contacts)

  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  cat("Fetching Data ...")
  if (split == FALSE) {
    cutoff <- 0.3
    pfracsfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
    pidsfold <- fetchPDBIDs(dataset, username)

    cutoff <- -1.0
    pfracsunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
    pidsunfold <- fetchPDBIDs(dataset, username)

    pfracs <- list(pfracsfold, pfracsunfold)

  } else {
    
    cutoff <- 0.3
    countsurf <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)
    pidsfold <- fetchPDBIDs(dataset, username)

    cutoff <- -1.0
    countunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)
    pidsunfold <- fetchPDBIDs(dataset, username)

    countburied <- countunfold - countsurf
    surfrac <- sum(countsurf) / (sum(countsurf) + sum(countburied))

    counts <- list(countsurf, countburied, countunfold)
    pfracs <- counts
    names(pfracs) <- c("psurf","pburied","punfold")
    
    for (i in 1:3) {
      for (j in 1:nrow(pfracs[[i]])) {
        pfracs[[i]][j,] <- pfracs[[i]][j,] / sum(pfracs[[i]][j,])
      }
    }
    
  }

  cat("Done! \n")
  ddG <- array(0,2)
  names(ddG) <- c("Close","Open")
  ddG[1] <- freeEnergyModel(countMatrix, groDist["GroEL_Close", ], pfracs, lambdaf, lambdau, split, surfrac)
  ddG[2] <- freeEnergyModel(countMatrix, groDist["GroEL_Open", ], pfracs, lambdaf, lambdau, split, surfrac)

  return(ddG)
}

bootstrapEnergyCycle <- function(dataset, bootstrap=500, username=myUsername, contacts=fetchContacts(paste(dataset, "_surface_contacts.csv",sep=""), username), lambdaf=1, lambdau=1) {

    countMatrix <- sampleContacts(contacts)
    
    cutoff <- 0.3
    pfracsfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)
    pids <- fetchPDBIDs(dataset, username)

    cutoff <- -1.0
    pfracsunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)

    ddG <- empty.df(c("Close", "Open", "Diff"), 1:bootstrap)
    cat("Bootstrapping\n")
  for(i in 1:bootstrap) {
    cat(paste("\rbootstrap:", i, "/", bootstrap))
    pfracs <- list(pfracsfold[pids[sample(length(pids), replace=T)], ], pfracsunfold[pids[sample(length(pids), replace=T)], ])

    ddG[i,"Close"] <- freeEnergyModel(countMatrix, groDist["GroEL_Close", ], pfracs, lambdaf, lambdau, split=FALSE)
    ddG[i,"Open"] <- freeEnergyModel(countMatrix, groDist["GroEL_Open", ], pfracs, lambdaf, lambdau, split=FALSE)
    ddG[i,"Diff"] <- ddG[i,"Close"] - ddG[i,"Open"]
  }
    cat("done\n")
    return(ddG)
}

proteinEnergyCycle <- function(username, dataset1="ecoli", dataset2="assist", contacts=fetchContacts(paste(dataset1, "_surface_contacts.csv",sep=""), username), pidsecoli=NULL, pidsassist=NULL, sample=FALSE) {
  countMatrix <- sampleContacts(contacts)

  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  if(sample == FALSE) {
    pidsecoli <- fetchPDBIDs(dataset1, username)
    pidsassist <- fetchPDBIDs(dataset2, username)
  } else {
    if(is.null(pidsecoli) & is.null(pidsassist)) {
      pidsecoli <- fetchPDBIDs(dataset1, username)
      pidsassist <- fetchPDBIDs(dataset2, username)
      indices <- sample(length(pidsecoli),replace=TRUE)
      indicesassist <- sample(length(pidsassist),replace=TRUE)
      pidsecoli <- pidsecoli[indices]
      pidsassist <- pidsassist[indicesassist]
      save(pidsecoli,file="pidsecoli.txt")
      save(pidsassist,file="pidsassist.txt")
    }
  }
  cat("Fetching Data...")
  
  cutoff <- 0.3
  csurf <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  csurf <- csurf[pidsecoli,]
  csurfassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)
  csurfassist <- csurfassist[pidsassist,]
  
  cutoff <- -1.0
  cunfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  cunfold <- cunfold[pidsecoli,]
  cunfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)
  cunfoldassist <- cunfoldassist[pidsassist,]
  
  cburied <- cunfold - csurf
  cburiedassist <- cunfoldassist - csurfassist
  cat(" Done!\n")
  
  counts <- list(csurf, cburied, cunfold)
  names(counts) <- c("csurf","cburied","cunfold")
  for (i in 1:3) {
    counts[[i]] <- counts[[i]][match(colnames(countMatrix),names(counts[[i]]))]
  }

  countsassist <- list(csurfassist, cburiedassist, cunfoldassist)
  names(countsassist) <- c("csurfassist","cburiedassist","cunfoldassist")
  for (i in 1:3) {
    countsassist[[i]] <- countsassist[[i]][match(colnames(countMatrix),names(countsassist[[i]]))]
  }
  
  pfracs <- counts
  names(pfracs) <- c("psurf","pburied","punfold")
  for (i in 1:3) {
    for (j in 1:nrow(pfracs[[i]])) {
      pfracs[[i]][j,] <- pfracs[[i]][j,] / sum(pfracs[[i]][j,])
    }
  }

  pfracsassist <- countsassist
  names(pfracsassist) <- c("psurfassist","pburiedassist","punfoldassist")
  for (i in 1:3) {
    for (j in 1:nrow(pfracsassist[[i]])) {
      pfracsassist[[i]][j,] <- pfracsassist[[i]][j,] / sum(pfracsassist[[i]][j,])
    }
  }
    
  ddGecoli <- matrix(0,length(pidsecoli),2)
  rownames(ddGecoli) <- pidsecoli
  colnames(ddGecoli) <- c("DDG","length")
  cat("Processing DataSet1...")

  for(i in 1:length(pidsecoli)) {
    if(pidsecoli[i] %in% pidsassist) {
      ddGecoli[i,1] <- 0
      cat(paste("\rSkipping...","         ", i,"/",length(pidsecoli)))
    }
    else {
      ener <- 0.
#      ener <- proteinFreeEnergyModel(pidsecoli[i], countMatrix, groDist["GroEL_Close", ], pfracs, counts, hydration=TRUE)
      ener <- proteinFreeEnergyModel(pidsecoli[i], countMatrix, groDist["GroEL_Open", ], pfracs, counts, hydration=TRUE)
      ddGecoli[i,1] <- -ener
      cat(paste("\rProcessing DataSet1...", i,"/",length(pidsecoli)))
    }
    ddGecoli[i,2] <- sum(counts$cunfold[pidsecoli[i],])
  }
  cat(" \n")
  ddGecoli <- ddGecoli[-which(ddGecoli[,1] == 0),]

  ddGassist <- matrix(0,length(pidsassist),2)
  rownames(ddGassist) <- pidsassist
  colnames(ddGassist) <- c("DDG","length")
  cat("Processing DataSet2...")

  for(i in 1:length(pidsassist)) {
    ener <- 0.
#    ener <- proteinFreeEnergyModel(pidsassist[i], countMatrix, groDist["GroEL_Close", ], pfracsassist, countsassist, hydration=TRUE)
    ener <- proteinFreeEnergyModel(pidsassist[i], countMatrix, groDist["GroEL_Open", ], pfracsassist, countsassist, hydration=TRUE)
    ddGassist[i,1] <- -ener
    cat(paste("\rProcessing DataSet2...", i,"/",length(pidsassist)))
    ddGassist[i,2] <- sum(countsassist$cunfoldassist[pidsassist[i],])
  }
  cat(" \n")
  cat("Completed! \n")
  ddG <- list(ddGecoli, ddGassist)
  names(ddG) <- c("ecoli","assist")

  return(ddG)
}

#This function is written on Aug 18th, 2011, based on hydration factors
freeEnergyModel <- function(countMatrix,yDist,pfracs,lambdaf,lambdau,split, surfrac) {

  hydration <- array(0,ncol(countMatrix))
  contact <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    hydration[l] <- countMatrix["WATER",l] / sum(countMatrix[c(1:20,which(rownames(countMatrix) == "WATER")),l]) * (1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]) + countMatrix["FREE",l] / countMatrix["TOTAL",l]
    contact[l] <- 1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]
  }
  names(hydration) <- colnames(countMatrix)
  names(contact) <- colnames(countMatrix)
  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )

  if (split == FALSE) { 
    sigmaXln <- array(0,2)
    for (i in 1:2) { 
      X <- array(0, ncol(pfracs[[i]]))
      for (j in 1:ncol(pfracs[[i]])) {
         Y <- rep(0, length(yDist))
         for (k in 1:length(yDist)) {
           t1 <- countMatrix[colnames(pfracs[[i]])[j], names(yDist)[k]] * yDist[k] * contact[names(yDist)[k]]
           t2 <- hydration[colnames(pfracs[[i]])[j]] / sum(hydration) * hydration[names(yDist)[k]] * yDist[k]
           Y[k] <- as.double(t1) + as.double(t2)
         }
         sigmaY <- log(sum(Y))
         X[j] <- median(pfracs[[i]][,j]) * sigmaY
       }
      sigmaXln[i] <- sum(X)      
    }
    deltaG = -(sigmaXln[1] * lambdaf - sigmaXln[2] * lambdau)

  } else {

    sigmaXln <- array(0,3)

    X <- array(0, ncol(pfracs$psurf))
    for (i in 1:ncol(pfracs$psurf)) {
      Y <- rep(0, length(yDist))
      for (j in 1:length(yDist)) {
        t1 <- countMatrix[colnames(pfracs$psurf)[i], names(yDist)[j]] * yDist[j] * contact[names(yDist)[j]]
        t2 <- hydration[colnames(pfracs$psurf)[i]] / sum(hydration) * hydration[names(yDist)[j]] * yDist[j]
        Y[j] <- as.double(t1) + as.double(t2)
      }
      sigmaY <- log(sum(Y))
      X[i] <- median(pfracs$psurf[,i]) * sigmaY
    }
    sigmaXln[1] <- sum(X) * surfrac

    X <- array(0, ncol(pfracs$pburied))
    for (i in 1:ncol(pfracs$pburied)) {
      Y <- rep(0, ncol(pfracs$pburied))
      for (j in 1:ncol(pfracs$pburied)) {
        t1 <- countMatrix[colnames(pfracs$pburied)[i], colnames(pfracs$pburied)[j]] * pfracs$pburied[i,j] * contact[colnames(pfracs$pburied)[j]]
        t2 <- hydration[colnames(pfracs$pburied)[i]] / sum(hydration) * hydration[colnames(pfracs$pburied)[j]] * pfracs$pburied[i,j]
        Y[j] <- as.double(t1) + as.double(t2)
      }
      sigmaY <- log(sum(Y))
      X[i] <- median(pfracs$pburied[,i]) * sigmaY
    }
    sigmaXln[2] <- sum(X) * (1-surfrac)

    X <- array(0, ncol(pfracs$punfold))
    for (i in 1:ncol(pfracs$punfold)) {
      Y <- rep(0, length(yDist))
      for (j in 1:length(yDist)) {
        t1 <- countMatrix[colnames(pfracs$punfold)[i], colnames(pfracs$punfold)[j]] * yDist[j] * contact[colnames(pfracs$punfold)[j]]
        t2 <- hydration[colnames(pfracs$punfold)[i]] / sum(hydration) * hydration[colnames(pfracs$punfold)[j]] * yDist[j]
        Y[j] <- as.double(t1) + as.double(t2)
      }
      sigmaY <- log(sum(Y))
      X[i] <- median(pfracs$punfold[,i]) * sigmaY
    }
    sigmaXln[3] <- sum(X)
    deltaG = -((sigmaXln[1]+sigmaXln[2]) * lambdaf - sigmaXln[3] * lambdau)
  }
  return (deltaG)
}

proteinFreeEnergyModel <- function(PDBID,countMatrix,yDist,pfracs,counts,hydration=FALSE) {
  if (hydration == TRUE) {
    
    hydration <- array(0,ncol(countMatrix))
    contact <- array(0,ncol(countMatrix))
    for (l in 1:ncol(countMatrix)) {
      hydration[l] <- countMatrix["WATER",l] / sum(countMatrix[c(1:20,which(rownames(countMatrix) == "WATER")),l]) * (1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]) + countMatrix["FREE",l] / countMatrix["TOTAL",l]
      contact[l] <- 1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]
    }
    names(hydration) <- colnames(countMatrix)
    names(contact) <- colnames(countMatrix)
  
    countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )
    yDist <- yDist[match(colnames(countMatrix),names(yDist))]
  
    ener <- c(0,0)
    names(ener) <- list("Gfold","Gunfold")
    sigmaXln <- array(0,2)
    sigmaX <- array(0,2)
    
    sum <- 0.
    for (i in 1:ncol(counts$csurf)) {
      pxy <- 0.
      for (j in 1:length(yDist)) {
        t1 <- countMatrix[i,j] * yDist[j] * contact[j]
        t2 <- hydration[i] / sum(hydration) * hydration[j] * yDist[j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(counts$csurf[PDBID,i] * sigmaY)
    }
    Xsurf <- sum

    sum <- 0.
    for (i in 1:ncol(counts$cburied)) {
      pxy <- 0.
      for (j in 1:ncol(pfracs$pburied)) {
        t1 <- countMatrix[i,j] * pfracs$pburied[PDBID,j] * contact[j]
        t2 <- hydration[i] / sum(hydration) * hydration[j] * pfracs$pburied[PDBID,j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(counts$cburied[PDBID,i] * sigmaY)
    }
    Xburied <- sum
    ener[1] <- Xsurf + Xburied

    sum <- 0.
    for (i in 1:ncol(counts$cunfold)) {
      pxy <- 0.
      for (j in 1:length(yDist)) {
        t1 <- countMatrix[i,j] * yDist[j] * contact[j]
        t2 <- hydration[i] / sum(hydration) * hydration[j] * yDist[j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(counts$cunfold[PDBID,i] * sigmaY)
    }
    Xunfold <- sum
    ener[2] <- Xunfold
   
  } else {

    countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )
    yDist <- yDist[match(colnames(countMatrix),names(yDist))]
  
    ener <- c(0,0)
    names(ener) <- list("Gfold","Gunfold")
    sigmaXln <- array(0,2)
    sigmaX <- array(0,2)

    sum <- 0.
    for (i in 1:ncol(counts$csurf)) {
      pxy <- 0.
      for (j in 1:length(yDist)) {
        pxy <- pxy + countMatrix[i,j] * yDist[j]
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(counts$csurf[PDBID,i] * sigmaY)
    }
    Xsurf <- sum

    sum <- 0.
    for (i in 1:ncol(counts$cburied)) {
      pxy <- 0.
      for (j in 1:ncol(pfracs$pburied)) {
        pxy <- pxy + countMatrix[i,j] * pfracs$pburied[PDBID,j]
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(counts$cburied[PDBID,i] * sigmaY)
    }
    Xburied <- sum
    ener[1] <- Xsurf + Xburied

    sum <- 0.
    for (i in 1:ncol(counts$cunfold)) {
      pxy <- 0.
      for (j in 1:length(yDist)) {
        pxy <- pxy + countMatrix[i,j] * yDist[j]
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(counts$cunfold[PDBID,i] * sigmaY)
    }
    Xunfold <- sum
    ener[2] <- Xunfold
    
  }

  return(ener[1] - ener[2])
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

hydrationModel <- function(countMatrix, yDist) {
  contact <- array(0,ncol(countMatrix))
  hydration <- array(0,ncol(countMatrix))
  for (l in 1:ncol(countMatrix)) {
    hydration[l] <- countMatrix["WATER",l] / sum(countMatrix[c(1:20,which(rownames(countMatrix) == "WATER")),l]) * (1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]) + countMatrix["FREE",l] / countMatrix["TOTAL",l]
  contact[l] <- 1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]
  }
  names(hydration) <- colnames(countMatrix)
  names(contact) <- colnames(countMatrix)
  yDist <- yDist[match(colnames(countMatrix),names(yDist))]

  g <- 0
  for (i in 1:length(yDist)) {
     g <- g + yDist[i] * log(hydration[i])
   }

  return (g)
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
      if (pxy < 0) {
        pxy <- 1E-199
      }
      sum <- sum + log(pxy) * deltaPs[i]
    }
    return(-sum)
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

  constrEq <- function(dist) {    
    sum <- sum(dist) - 1
    return (sum)
  }

  constrEq.jac <- function(dist) {
    j <- matrix(0, 1 ,length(dist))
    j[1, ] <- array(1, 20)
  }

#  constrIeq <- function(dist) {
#    sum <- sum(dist[which(dist < 0)])
#    return (sum)
#  }

  constrIeq <- function(dist) {
    h <- rep(NA,1)
    listSp <- c(5,6,18)
    list <- c(1,2,3,4,7,8,9,10,11,12,13,14,15,16,17,19,20)
    for (j in list) {
        h[j] <- dist[j]
      }
    for (i in listSp) {
        h[i] <- 0.002 - dist[i]
      }
    return (h)
  }
  
  constrIeq.jac <- function(dist) {
    j <- matrix(0, length(dist), length(dist))
    for (i in 1:length(dist)) {
      j[i,i] <- 1
    }
    return (j)
  }

  ini <- c(0.1, 0.05, 0.01, 0.15, 0.001, 0.001, 0.15, 0.15, 0.01, 0.01, 0.02, 0.15, 0, 0.03, 0.04, 0.07, 0.06, 0.001, 0.04, 0.02)
#  ini <- runif(20,0,0.05)
#  ini[c(5,6,18)] <- 0.001
#  ini <- array(0.1,20)
#  ini <- as.numeric(groDist["GroEL_Close",])
  names(ini) <- names(groDist["GroEL_Close",])
  
#  opt <- optim(ini, g, devI, method = "L-BFGS-B", lower = array(0,20), upper = array(0.5,20), control = list(trace = T))
  optimiz <- function(ini, iter)
    for (i in 1:iter) {
      if (i == 1) {
        ini <- ini
      } else {
        ini <- opt$par
      }
    optMin <- auglag(par = ini, fn = g, gr = devI, hin = constrIeq, hin.jac = constrIeq.jac, heq = constrEq, heq.jac = constrEq.jac, control.outer = list(trace = T, eqs = 1E-100, mu0 = 1E-2, method = "BFGS"), control.optim = list(fnscale = 1))
#    print(opt[1])
#    print(optMin)
#    optMax <- auglag(par = ini, fn = g, gr = devI, hin = constrIeq, hin.jac = constrIeq.jac, heq = constrEq, heq.jac = constrEq.jac, control.outer = list(trace = T, eqs = 1E-100, method = "BFGS"), control.optim = list(fnscale = -1))
#    print(optMax)
    }
#  print(opt[1])
#  print(opt[2])
#  print(opt[3])
#  optimiz(ini, 1)
#  print(groDist["GroEL_Close",])
#  print(g(groDist["GroEL_Close",]))

  print(devI(groDist["GroEL_Close",]))

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

#  print(numDev(groDist["GroEL_Open",], 0.0001))
}

normalPlot <- function(ddG,name) {
  x <- c(ddG$ecoli[,1],ddG$assist[,1])
  y <- c(ddG$ecoli[,2],ddG$assist[,2])
  col <- rep("red",nrow(ddG$ecoli))
  col <- c(col,rep("blue",nrow(ddG$assist)))
  cairo_pdf(paste(name,".pdf",sep=""), width=7, height=5)
  par(family="LMRoman10")
  plot(x,y,type="p",col=col,xlab="ddG",ylab="length",xlim=c(quantile(x[which(y<1000)],probs=c(0.05,0.95))["5%"],quantile(x[which(y<1000)],probs=c(0.05,0.95))["95%"]),ylim=c(0,1000))
  dev.off()
}

plotPS <- function(x, y,xpoints=NULL,ypoints=NULL, xlab, ylab,  plotName1, plotName2) {

  xlim <- c(quantile(x[which(y<1000)],probs=c(0.025,0.975))[1],quantile(x[which(y<1000)],probs=c(0.025,0.975))[2])
  ylim <- c(0,1000)
#  xlim <- c(quantile(x,probs=c(0.025,0.975))[1],quantile(x,probs=c(0.025,0.975))[2])
#  ylim <- c(quantile(y,probs=c(0.025,0.975))[1],quantile(y,probs=c(0.025,0.975))[2])
#  xlim <- c(-100,0)
#  ylim <- c(0,8000)

  mest <- bkde2D(x=cbind(x,y), bandwidth=c(2,50), gridsize=c(1000,1000), range.x=list(xlim, ylim))
  print(paste(mest$x1[floor(which.max(mest$fhat) %% length(mest$x1))], mest$x2[floor(which.max(mest$fhat) / length(mest$x2))]))
  nlevels <- 30

  colGrad <- colorRampPalette(c(
  rgb(0.0, 0.000, 0.0),
  rgb(0.142, 0.000, 0.850),
  rgb(0.097, 0.112, 0.970),
  rgb(0.160, 0.342, 1.000),
  rgb(0.240, 0.531, 1.000),
  rgb(0.340, 0.692, 1.000),
  rgb(0.460, 0.829, 1.000),
  rgb(0.600, 0.920, 1.000),
  rgb(0.740, 0.978, 1.000),
  rgb(0.920, 1.000, 1.000),
  rgb(1.000, 1.000, 0.920),
  rgb(1.000, 0.948, 0.740),
  rgb(1.000, 0.840, 0.600),
  rgb(1.000, 0.676, 0.460),
  rgb(1.000, 0.472, 0.340),
  rgb(1.000, 0.240, 0.240),
  rgb(0.970, 0.155, 0.210),
  rgb(0.850, 0.085, 0.187),
  rgb(0.650, 0.000, 0.130)));

  colors <- colGrad(nlevels * 10)
  
  png(paste(plotName1, ".png", sep=""), width=1750, height=1750, res=250)
  par(family="LMRoman10", fg="dark gray")
  image(mest$x1, mest$x2, mest$fhat, col=colors, bg="black", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
  points(x=xpoints,y=ypoints,type="p",col="green", pch=19)
  graphics.off()

  png(paste(plotName2,".png",sep=""),width=1750, height=1750, res=250)
  par(family="LMRoman10", fg="dark gray")
  plot(density(x), xlab=xlab, col="red", type="l", lwd=4, main="Density Plot", xlim=c(-35,5),ylim=c(0,0.15))
#  points(x=density(xpoints),type="l",lwd=4,col="green") 
  for (i in 1:length(xpoints)) {
    abline(v = xpoints[i])
  }
  legend("topleft", col=c("red", "green"), legend=c("E.Coli Density", "Assisted Points"), pch=15)
  graphics.off()

}

ellipsoidPlot <- function(ddG,name) {
  x <- c(ddG$ecoli[,1],ddG$assist[,1])
  y <- c(ddG$ecoli[,2],ddG$assist[,2])
  col <- rep("red",nrow(ddG$ecoli))
  col <- c(col,rep("blue",nrow(ddG$assist)))

  r1 <- (quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),1],prob=c(0.025,0.5,0.975))[3] - quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),1],prob=c(0.025,0.5,0.975))[1]) / 2
  r2 <- (quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),2],prob=c(0.025,0.5,0.975))[3] - quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),2],prob=c(0.025,0.5,0.975))[1]) / 2
  m1 <- (quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),1],prob=c(0.025,0.5,0.975))[3] + quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),1],prob=c(0.025,0.5,0.975))[1]) / 2
  m2 <- (quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),2],prob=c(0.025,0.5,0.975))[3] + quantile(ddG$ecoli[which(ddG$ecoli[,2]<1000),2],prob=c(0.025,0.5,0.975))[1]) / 2

  a <- seq(-21.7,7.7,0.01)
  b1 <- sqrt((1 - (((a-m1)**2) / (r1**2))) * (r2**2)) + m2
  b2 <- -sqrt((1 - (((a-m1)**2) / (r1**2))) * (r2**2)) + m2

  a <- c(a,a)
  b <- c(b1,b2)
  cairo_pdf(paste(name,".pdf",sep=""), width=7, height=5)
  par(family="LMRoman10")
  plot(x,y,type="p",col=col,xlab="ddG",ylab="length",xlim=c(quantile(x[which(y<1000)],probs=c(0.025,0.975))[1],quantile(x[which(y<1000)],probs=c(0.025,0.975))[2]),ylim=c(0,1000))
  points(x=a,y=b,type="p",col="green")
  dev.off()
}

normPlot <- function(data1, data2, plotName) {
  png(paste(plotName, ".png", sep=""), width=1750, height=1750, res=250)
  par(family="LMRoman10", fg="dark gray")
  plot(x=seq(min(data1[,1]),max(data1[,1]),0.01),y=dnorm(seq(min(data1[,1]),max(data1[,1]),0.01), mean(data1[,1]), sqrt(var(data1[,1]))),type="l", lwd=4, col="red", xlab="ddG",ylab="Density", ylim=c(0,0.08),xlim=c(-40,40))
  points(x=seq(min(data1[,1]),max(data1[,1]),0.01),y=dnorm(seq(min(data1[,1]),max(data1[,1]),0.01), mean(data2 [,1]), sqrt(var(data2[,1]))), type="l", lwd=4, col="green")
  legend("topleft", col=c("red", "green"), legend=c("E.Coli Distribution", "Assisted Distribution"), pch=15)
  graphics.off()
}

phobicity <- function(username, dataset1="ecoli", dataset2="assist") {
  pidsecoli <- fetchPDBIDs(dataset1, username)
  cat("Fetching Data...")
  
  cutoff <- -1.0
  cunfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  cunfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)

  pho1 <- matrix(0,nrow(cunfold),2)
  rownames(pho1) <- rownames(cunfold)
  colnames(pho1) <- c("phobicity","length")

  pho2 <- matrix(0,nrow(cunfoldassist),2)
  rownames(pho2) <- rownames(cunfoldassist)
  colnames(pho2) <- c("phobicity","length")

  cat(" Done!\n")
  cat("Processing...")
  
  for (i in 1:nrow(cunfold)) {
    if(rownames(cunfold)[i] %in% rownames(cunfoldassist)) {
      pho1[i,1] <- 0
      cat(paste("\rSkipping...","         ", i,"/",nrow(cunfold)))
    }
    else {
      cat(paste("\rProcessing DataSet1...", i,"/",nrow(cunfold)))
      foo <- sum(cunfold[i,c("MET","ALA","VAL","LEU","ILE","PRO")]) / sum(cunfold[i,])
      pho1[i,1] <- log(foo)
      pho1[i,2] <- sum(cunfold[i,])
    }
  }
  cat(" \n")
  pho1 <- pho1[-which(pho1[,1] == 0),]

  cat("Processing DataSet2...")
  for (i in 1:nrow(cunfoldassist)) {
    cat(paste("\rProcessing DataSet2...", i,"/",nrow(cunfoldassist)))
    foo <- sum(cunfoldassist[i,c("MET","ALA","VAL","LEU","ILE","PRO")]) / sum(cunfoldassist[i,])
    pho2[i,1] <- log(foo)
    pho2[i,2] <- sum(cunfoldassist[i,])
  }
  cat(" \n")
  cat("Completed!\n")

  data <- list(pho1,pho2)
  names(data) <- c("ecoli","assist")
  
  return(data)
}

netCharge <- function(username, dataset1="ecoli", dataset2="assist") {
  pidsecoli <- fetchPDBIDs(dataset1, username)
  cat("Fetching Data...")
  
  cutoff <- -1.0
  cunfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  cunfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)

  cha1 <- matrix(0,nrow(cunfold),2)
  rownames(cha1) <- rownames(cunfold)
  colnames(cha1) <- c("netcharge","length")

  cha2 <- matrix(0,nrow(cunfoldassist),2)
  rownames(cha2) <- rownames(cunfoldassist)
  colnames(cha2) <- c("netcharge","length")

  cat(" Done!\n")
  cat("Processing...")
  
  for (i in 1:nrow(cunfold)) {
    if(rownames(cunfold)[i] %in% rownames(cunfoldassist)) {
      cha1[i,1] <- 999
      cat(paste("\rSkipping...","         ", i,"/",nrow(cunfold)))
    }
    else {
      cat(paste("\rProcessing DataSet1...", i,"/",nrow(cunfold)))
      cha1[i,1] <- sum(cunfold[i,c("LYS","ARG")]) - sum(cunfold[i, c("ASP","GLU")])
      cha1[i,2] <- sum(cunfold[i,])
    }
  }
  cat(" \n")
  cha1 <- cha1[-which(cha1[,1] == 999),]

  cat("Processing DataSet2...")
  for (i in 1:nrow(cunfoldassist)) {
    cat(paste("\rProcessing DataSet2...", i,"/",nrow(cunfoldassist)))
    cha2[i,1] <- sum(cunfoldassist[i,c("LYS","ARG")]) - sum(cunfoldassist[i, c("ASP","GLU")])
    cha2[i,2] <- sum(cunfoldassist[i,])
  }
  cat(" \n")
  cat("Completed!\n")

  surfCharge <- list(cha1,cha2)
  names(surfCharge) <- c("ecoli","assist")
  
  return(surfCharge)
}

cysFrac <- function(username, dataset1="ecoli", dataset2="assist") {
  pidsecoli <- fetchPDBIDs(dataset1, username)
  cat("Fetching Data...")
  
  cutoff <- -1.0
  cunfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  cunfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)

  cys1 <- matrix(0,nrow(cunfold),2)
  rownames(cys1) <- rownames(cunfold)
  colnames(cys1) <- c("CYS fraction","length")

  cys2 <- matrix(0,nrow(cunfoldassist),2)
  rownames(cys2) <- rownames(cunfoldassist)
  colnames(cys2) <- c("CYS fraction","length")

  cat(" Done!\n")
  cat("Processing...")
  
  for (i in 1:nrow(cunfold)) {
    if(rownames(cunfold)[i] %in% rownames(cunfoldassist)) {
      cys1[i,1] <- 999
      cat(paste("\rSkipping...","         ", i,"/",nrow(cunfold)))
    }
    else {
      cat(paste("\rProcessing DataSet1...", i,"/",nrow(cunfold)))
      cys1[i,1] <- cunfold[i,"CYS"] / sum(cunfold[i,])
      cys1[i,2] <- sum(cunfold[i,])
    }
  }
  cat(" \n")
  cys1 <- cys1[-which(cys1[,1] == 999),]

  cat("Processing DataSet2...")
  for (i in 1:nrow(cunfoldassist)) {
    cat(paste("\rProcessing DataSet2...", i,"/",nrow(cunfoldassist)))
    cys2[i,1] <- cunfoldassist[i,"CYS"] / sum(cunfoldassist[i,])
    cys2[i,2] <- sum(cunfoldassist[i,])
  }
  cat(" \n")
  cat("Completed!\n")

  cysFraction <- list(cys1,cys2)
  names(cysFraction) <- c("ecoli","assist")
  
  return(cysFraction)
}


#Below are used for actuall program running
#ddG <- energyCycle("ecoli","wenjunh",split=TRUE)
#load("pidsecoli.txt")
#load("pidsassist.txt")
#ddG2 <- proteinEnergyCycle("wenjunh")#, pidsecoli=pidsecoli, pidsassist=pidsassist)
#ddG_assist <- bootstrapEnergyCycle("assist",bootstrap=100, username="wenjunh", contacts=fetchContacts("ecoli_surface_contacts.csv", "wenjunh"))
#save(ddG_assist, file="ddG_boots_assist.Rdata")
#ddG_nonassist <- bootstrapEnergyCycle("nonassist",bootstrap=100, username="whitead", contacts=fetchContacts("ecoli_surface_contacts.csv", "wenjunh"))
#save(ddG_nonassist, file="ddG_boots_nonassist.Rdata")
#ddG
#q(save="yes")
#normalPlot(ddG,"ddGquan")
#plotPS(x=ddG$ecoli[,1],y=ddG$ecoli[,2],xpoints=ddG$assist[,1],ypoints=ddG$assist[,2], xlab="ddG",ylab="length",plotName1="trial1",plotName2="trial2")
#ellipsoidPlot(ddG,"ellipsoid")
#plotPS(x=data$ecoli[,1],y=data$ecoli[,2],xpoints=data$assist[,1],ypoints=data$assist[,2], xlab="phobicity",ylab="length",plotName1="trial1", plotName2="trial2")
#data <- phobicity("wenjunh")
#surfCharge <- netCharge("wenjunh")
#cys <- cysFrac("wenjunh")
#plotPS(x=surfCharge$ecoli[,1],y=surfCharge$ecoli[,2],xpoints=surfCharge$assist[,1],ypoints=surfCharge$assist[,2], xlab="netCharge",ylab="length",plotName1="trial1", plotName2="trial2")
#normPlot(ddG$ecoli,ddG$assist,"trial")
#plotPS(x=surfCharge$ecoli[,1],y=data$ecoli[,1],xpoints=surfCharge$assist[,1],ypoints=data$assist[,1], xlab="netCharge",ylab="Hydrophobicity",plotName1="trial1", plotName2="trial2")
#energyCycle("assist", username="wenjunh", contacts=fetchContacts("ecoli_surface_contacts.csv", "wenjunh"))
#energyBootstrap(1000, "ecoli", username="wenjunh")
#minimizeEnergy("ecoli", "wenjunh")


#Used to output value for python code
contacts = fetchContacts("ecoli_surface_contacts.csv", "wenjunh")
countMatrix <- sampleContacts(contacts)

gfold <- hydrationModel(countMatrix,groDist["GroEL_Open",])
gunfold <- hydrationModel(countMatrix,groDist["GroEL_Close",])
dG <- gfold-gunfold

#hydration <- array(0,ncol(countMatrix))
#  contact <- array(0,ncol(countMatrix))
#  for (l in 1:ncol(countMatrix)) {
#    hydration[l] <- countMatrix["WATER",l] / sum(countMatrix[c(1:20,which(rownames(countMatrix) == "WATER")),l]) * (1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]) + countMatrix["FREE",l] / countMatrix["TOTAL",l]
#    contact[l] <- 1. - countMatrix["FREE",l] / countMatrix["TOTAL",l]
#  }
#  names(hydration) <- colnames(countMatrix)
#  names(contact) <- colnames(countMatrix)
#  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )
#write.table(countMatrix, file = "countMatrix.txt", row.names = FALSE, col.names = FALSE)
#write.table(contact, file = "contact.txt", row.names = FALSE, col.names = FALSE)
#write.table(hydration, file = "hydration.txt", row.names = FALSE, col.names = FALSE)

#dataset <- "ecoli"
#username <- "wenjunh"

#cutoff <- 0.3
#countsurf <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)
#pidsfold <- fetchPDBIDs(dataset, username)

#cutoff <- -1.0
#countunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)
#pidsunfold <- fetchPDBIDs(dataset, username)

#countburied <- countunfold - countsurf
#surfrac <- sum(countsurf) / (sum(countsurf) + sum(countburied))

#counts <- list(countsurf, countburied, countunfold)
#pfracs <- counts
#names(pfracs) <- c("psurf","pburied","punfold")
    
#for (i in 1:3) {
#  for (j in 1:nrow(pfracs[[i]])) {
#    pfracs[[i]][j,] <- pfracs[[i]][j,] / sum(pfracs[[i]][j,])
#  }
#}

#groDist["GroEL_Open", ] <- groDist["GroEL_Open",match(colnames(countMatrix),nam#es(groDist["GroEL_Open",]))]
#groDist["GroEL_Close", ] <- groDist["GroEL_Close",match(colnames(countMatrix),names(groDist["GroEL_Close",]))]
#colnames(groDist) <- colnames(countMatrix)

#write(surfrac, file = "surfrac.txt")
#write.table(pfracs$psurf, file = "psurf.txt", row.names = FALSE, col.names = FALSE)
#write.table(pfracs$pburied, file = "pburied.txt", row.names = FALSE, col.names = FALSE)
#write.table(pfracs$punfold, file = "punfold.txt", row.names = FALSE, col.names = FALSE)
#write.table(groDist["GroEL_Open", ], file = "GroOp.txt", row.names = FALSE, col.names = FALSE)
#write.table(groDist["GroEL_Close", ], file = "GroCl.txt", row.names = FALSE, col.names = FALSE)
