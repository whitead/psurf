library(alabama)
library(KernSmooth)
source("SQLShareLib.R")

##Obtain the raw counts of GroEL inside surface residues, both open and close form
sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))


##This is the main code in this script for free energy model at 903 individual protein level,  It is linked to the "proteinFreeEnergyModel" script below 
proteinEnergyCycle <- function(username, dataset1="ecoli40", dataset2="assist", contacts=fetchContacts(paste("ecoli", "_total_contacts.csv",sep=""), username), countMatrix = FALSE, pidsecoli=NULL, pidsassist=NULL, sample=FALSE) {
  if(countMatrix == FALSE) {
    countMatrix <- sampleContacts(contacts)
  }
  #Draw contact matrix from the SQLShare if no matrix is specified
  
  #obtain surface residues for both E.Coli fold and E.Coli unfold
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
  
  cutoff <- 0.3  #Set the surface cutoff
  csurf <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  csurf <- csurf[pidsecoli,]
  csurfassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)
  csurfassist <- csurfassist[pidsassist,]
  
  cutoff <- -1.0  #Set the unfolded cutoff
  cunfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  cunfold <- cunfold[pidsecoli,]
  cunfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)
  cunfoldassist <- cunfoldassist[pidsassist,]
  
  cburied <- cunfold - csurf
  cburiedassist <- cunfoldassist - csurfassist
  cat(" Done!\n")

  #Generate the raw counts lists
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

  #Generate the probability fraction list
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

  #Generate surface fraction for both datasets
  surFrac <- sum(counts$csurf) / (sum(counts$csurf) + sum(counts$cburied))
  surFracAssist <- sum(countsassist$csurfassist) / (sum(countsassist$csurfassist) + sum(countsassist$cburiedassist));

  #Generate ddG matrix
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
      #Use the "proteinFreeEnergyModel" script to calculate ddG value for specific protein using either GroEL_open distribution or GroEL_close distribution.
      ener <- 0.
      ener <- proteinFreeEnergyModel(pidsecoli[i], countMatrix, groDist["GroEL_Close", ], pfracs, counts, surFrac, hydration=TRUE)
#      ener <- proteinFreeEnergyModel(pidsecoli[i], countMatrix, groDist["GroEL_Open", ], pfracs, counts, surFrac, hydration=TRUE)
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
    ener <- proteinFreeEnergyModel(pidsassist[i], countMatrix, groDist["GroEL_Close", ], pfracsassist, countsassist, surFracAssist, hydration=TRUE)
#    ener <- proteinFreeEnergyModel(pidsassist[i], countMatrix, groDist["GroEL_Open", ], pfracsassist, countsassist, surFracAssist, hydration=TRUE)
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



#get the interaction energies between residue types
getInteractionEnergys <- function(dataset, username=NULL) {

  dataset <- paste(paste(dataset,"surface","contacts", sep="_"), "csv" ,sep=".")
 
  #get the data
  if(is.null(username)) {
    countMatrix <- fetchContacts(dataset)
  }
  else{
    countMatrix <- fetchContacts(dataset, username)
  }

  #re-order it
  anames <- colnames(countMatrix[-c(1,2)])
  anames.ord <- order(anames)

  countMatrix <- countMatrix[c(1,2,anames.ord + 2)]

  #sum it
  contactMatrix <- sampleContacts(countMatrix, random=FALSE)

  contactMatrix <- contactMatrix[order(rownames(contactMatrix)[1:20]),]

  #normalize it
  mat <- t(apply(contactMatrix, MARGIN=1, FUN=function(x){ x / sum(x)}))

  print(mat)
  print(sum(mat[1,]))

  mat <- -log(mat)
  return(mat)
}

#The mean enegycycle code, worked for ecoli, assisted proteins, also used by the stochastic annealing code (Python).
proteinFreeEnergyModel <- function(PDBID,countMatrix,yDist,pfracs,counts,surFrac,hydration=FALSE) {
  if (hydration == TRUE) {
    #process the contact matrix, making countMatrix, hydration, and contact.
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
    #Calculation of ddG for individual protein
    sum <- 0.
    for (i in 1:ncol(pfracs$psurf)) {
      pxy <- 0.
      for (j in 1:length(yDist)) {
        t1 <- countMatrix[i,j] * yDist[j] * contact[j]
        t2 <- hydration[i] / sum(hydration) * hydration[j] * yDist[j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(pfracs$psurf[PDBID,i] * sigmaY)
    }
    Xsurf <- sum * surFrac

    sum <- 0.
    for (i in 1:ncol(pfracs$pburied)) {
      pxy <- 0.
      for (j in 1:ncol(pfracs$pburied)) {
        t1 <- countMatrix[i,j] * pfracs$pburied[PDBID,j] * contact[j]
        t2 <- hydration[i] / sum(hydration) * hydration[j] * pfracs$pburied[PDBID,j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(pfracs$pburied[PDBID,i] * sigmaY)
    }
    Xburied <- sum * (1-surFrac)
    ener[1] <- Xsurf + Xburied

    sum <- 0.
    for (i in 1:ncol(pfracs$punfold)) {
      pxy <- 0.
      for (j in 1:length(yDist)) {
        t1 <- countMatrix[i,j] * yDist[j] * contact[j]
        t2 <- hydration[i] / sum(hydration) * hydration[j] * yDist[j]
        pxy <- pxy + as.double(t1) + as.double(t2)
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(pfracs$punfold[PDBID,i] * sigmaY)
    }
    Xunfold <- sum * surFrac
    ener[2] <- Xunfold
   
  } else {

    countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )
    yDist <- yDist[match(colnames(countMatrix),names(yDist))]
    surFrac <- sum(counts$surf) / (sum(counts$surf) + sum(counts$buried))
  
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
    Xsurf <- sum * surFrac

    sum <- 0.
    for (i in 1:ncol(counts$cburied)) {
      pxy <- 0.
      for (j in 1:ncol(pfracs$pburied)) {
        pxy <- pxy + countMatrix[i,j] * pfracs$pburied[PDBID,j]
      }
      sigmaY <- log(pxy)
      sum <- sum + as.double(counts$cburied[PDBID,i] * sigmaY)
    }
    Xburied <- sum * (1-surFrac)
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
    Xunfold <- sum * surFrac
    ener[2] <- Xunfold
    
  }
			
  return(ener[1] - ener[2])
}

##This part of the code is the calculation for ddG of all the E.Coli protein, available for incorprating lambdaf and lambdau
energyCycle <- function(dataset, username=myUsername, countMatrix=FALSE, lambdaf=1, lambdau=1, split=FALSE) {
  if (countMatrix == FALSE) {
    contacts <- fetchContacts(paste(dataset, "_total_contacts.csv",sep=""), username)
    countMatrix <- sampleContacts(contacts)
  }

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
    names(counts) <- c("csurf","cburied","cunfold")
    for (i in 1:3) {
      counts[[i]] <- counts[[i]][match(colnames(countMatrix),names(counts[[i]]))]
    }
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
  ddG[1] <- freeEnergyModel(countMatrix, groDist["GroEL_Close", ], pfracs, counts=counts,lambdaf, lambdau, split, surfrac)
  #ddG[2] <- freeEnergyModel(countMatrix, groDist["GroEL_Open", ], pfracs, counts=counts,lambdaf, lambdau, split, surfrac)

  return(ddG)
}

#This function is written on Aug 18th, 2011, based on hydration factors
freeEnergyModel <- function(countMatrix,yDist,pfracs,counts=NA,lambdaf,lambdau,split,surfrac) {
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
      X[i] <- mean(pfracs$psurf[,i]) * sigmaY
    }
    sigmaXln[1] <- sum(X) * surfrac

    X <- array(0, ncol(pfracs$pburied))
    for (i in 1:ncol(pfracs$pburied)) {
      Y <- rep(0, ncol(pfracs$pburied))
      for (j in 1:ncol(pfracs$pburied)) {
        t1 <- countMatrix[colnames(pfracs$pburied)[i], colnames(pfracs$pburied)[j]] * pfracs$pburied[i,j] * contact[colnames(pfracs$pburied)[j]]
        t2 <- hydration[colnames(pfracs$psurf)[i]] / sum(hydration) * hydration[names(yDist)[j]] * yDist[j]
        Y[j] <- as.double(t1) + as.double(t2)
      }
      sigmaY <- log(sum(Y))
      X[i] <- mean(pfracs$pburied[,i]) * sigmaY
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
      X[i] <- mean(pfracs$punfold[,i]) * sigmaY
    }
    sigmaXln[3] <- sum(X) * surfrac
    
    deltaG = ((sigmaXln[1]+sigmaXln[2]) * lambdaf - sigmaXln[3] * lambdau)
  }
  return (deltaG)
}


#Below are used for actuall program running
#ddG <- energyCycle("ecoli","wenjunh",split=TRUE)    

##Main excecuting code of this script, energycycle also used by Python code 
#load("pidsecoli.txt")
#load("pidsassist.txt")
#ddG1 <- proteinEnergyCycle("wenjunh")#, pidsecoli=pidsecoli, pidsassist=pidsassist)

#q(save="yes")

##Executing codes for the energy model that NOT USED any more
#ddG2 <- energyCycle("ecoli", username="wenjunh", split=TRUE)#,countMatrix=read.table('countMatrix.txt'))


##Used to output value for python code
#contacts = fetchContacts("ecoli_surface_contacts.csv", "wenjunh")
#countMatrix <- sampleContacts(contacts)
#
#
#
#hydration <- array(0,ncol(countMatrix))
#contact <- array(0,ncol(countMatrix))
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
#
#dataset <- "ecoli"
#username <- "wenjunh"
#
#cutoff <- 0.3
#countsurf <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)
#pidsfold <- fetchPDBIDs(dataset, username)
#
#cutoff <- -1.0
#countunfold <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)
#pidsunfold <- fetchPDBIDs(dataset, username)
#
#countburied <- countunfold - countsurf
#surfrac <- sum(countsurf) / (sum(countsurf) + sum(countburied))
#
#counts <- list(countsurf, countburied, countunfold)
#for (i in 1:3) {
#  counts[[i]] <- counts[[i]][match(colnames(countMatrix),names(counts[[i]]))]
#    }
#pfracs <- counts
#names(counts) <- c("csurf","cburied","cunfold")
#names(pfracs) <- c("psurf","pburied","punfold")
#    
#for (i in 1:3) {
#  for (j in 1:nrow(pfracs[[i]])) {
#    pfracs[[i]][j,] <- pfracs[[i]][j,] / sum(pfracs[[i]][j,])
#  }
#}
#print(pfracs$psurf)
#
#yDist <- groDist["GroEL_Close",]
#yDist <- yDist[match(colnames(countMatrix),names(yDist))]
#print(yDist)
#PDBID <- rownames(pfracs[[1]])
#
#proteinFreeEnergyModel(PDBID,countMatrix,yDist,pfracs,counts)
#groDist_Open <- groDist["GroEL_Open",]
#groDist_Open <- groDist_Open[match(colnames(countMatrix),names(groDist_Open))]
#groDist_Closed <- groDist["GroEL_Close",]
#groDist_Closed <- groDist_Closed[match(colnames(countMatrix),names(groDist_Closed))]
#colnames(groDist) <- colnames(countMatrix)

#write(surfrac, file = "surfrac.txt")
#write.table(pfracs$psurf, file = "psurf.txt", row.names = FALSE, col.names = FALSE)
#write.table(pfracs$pburied, file = "pburied.txt", row.names = FALSE, col.names = FALSE)
#write.table(pfracs$punfold, file = "punfold.txt", row.names = FALSE, col.names = FALSE)
#write.table(groDist_Open, file = "GroOp.txt", row.names = FALSE, col.names = FALSE)
#write.table(groDist_Closed, file = "GroCl.txt", row.names = FALSE, col.names = FALSE)

#end of script


##Functions that are no longer being used





