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


##This is the main code in this script for free energy model at 903 individual protein level, Modified Dec 19th, 2011
proteinEnergyCycle <- function(username, dataset1="ecoli", dataset2="assist") {
  surfRough <- 3
  charaLength <- 40000    #These numbers are subject to change
  thetaF <- 1             #Need a method/detail to get this number (for ennergy calcualtion)    
  confineExp <- 3.25
  
  gyration <- getGyrationRadius(dataset1, username)
  print(gyration[300,])

  Tao <- getContactFractionTime(gyration, surfRough, charaLength)

  lambda <- getAreaofContact(gyration, dataset1, username)

  Xmatrix <- getInteractionEnergys(dataset1, username)

  
  cat("Fetching Data...")
  
  cutoff <- 0.3  #Set the surface cutoff
  psurf <- fetchAllSurfResidues(dataset1, cutoff, normalize=TRUE, username)
  psurfassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=TRUE, username)
  
  cutoff <- -1.0  #Set the unfolded cutoff
  punfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=TRUE, username)
  punfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=TRUE, username)

  cat(" Done!\n")

  #Generate ddG matrix
  cat("Processing DataSet...")

  energy <- matrix(0,nrow(gyration),2)
  rownames(energy) <- rownames(psurf)
  colnames(energy) <- c("Efold","Eunfold")

  for(i in 1:nrow(psurf)) {
    PDBID <- rownames(psurf)[i]
    energy[i,"Efold"] <- getEnergyforFoldProtein(PDBID,psurf,groDist["GroEL_Close",],Xmatrix,thetaF)
    energy[i,"Eunfold"] <- getEnergyforUnfoldProtein(PDBID,punfold,groDist["GroEL_Close",],Xmatrix,thetaF,gyration)
  }

  #calculate the final DDG matrix
  ddG <- matrix(0,nrow(energy),1)
  rownames(ddG) <- rownames(energy)
  colnames <- "ddG"

  for (i in 1:nrow(ddG)) {
    ddG[i,1] <-  Tao[i,2] * lambda[i,2] * energy[i,2] - Tao[i,1] * lambda[i,1] * energy[i,1] - (gyration[i,1] / charaLength)^confineExp
  }


  return(ddG)
}
 

#get the general gyration radius (Rg = N^nv*l) (Rf = sqrt(lambda_x^2+lambda_y^2+lambda_z^2))
getGyrationRadius <- function(dataset, username) {

  #Get Rg
  Nv <- 0.6
  l <- 2.5    #based on estimation, subject to change in the future
  N <- fetchResNum(dataset, username)
  
  Rg <- matrix(0,nrow(N),1)
  rownames(Rg) <- rownames(N)
  
  for (i in 1:nrow(N)) {
    Rg[i,1] <- (as.numeric(N[i,1])^Nv)*l
  }
 
  #Get Rf
  lambda <- fetchGyration(dataset,username)

  Rf <- matrix(0,nrow(lambda),1)
  rownames(Rf) <-  rownames(lambda)
  
  for (i in 1:nrow(lambda)) {
    Rf[i,1] <- sqrt(as.numeric(lambda[i,1])^2+as.numeric(lambda[i,2])^2+as.numeric(lambda[i,3])^2)
  }

  #Put Rf, Rg into one matrix
  data <- empty.df(c("Rg","Rf"), rownames(lambda))
  for (i in 1:nrow(lambda)) {
    data[i,1] <- N[i,1]
    data[i,2] <- Rf[i,1]
  }

  return(data)
}

#Get the Tao value for both fold and unfolded states
getContactFractionTime <- function(gyration, surfRough, charaLength){
  Tao <- matrix(0,nrow(gyration),ncol(gyration))
  rownames(Tao) <- rownames(gyration)
  colnames(Tao) <- c("Tf","Tu")

  for (i in 1:nrow(Tao)) {
      Tao[i,1] <- 1 - (1 - surfRough / (charaLength - as.numeric(gyration[i,1])))^3
      Tao[i,2] <- 1 - (1 - surfRough / (charaLength - as.numeric(gyration[i,2])))^3
    }

  return(Tao)
}

#get the Area of Contact for both fold and unfolded state
getAreaofContact <- function(gyration,dataset,username) {
  #acquire gyration info
  lambda <- fetchGyration(dataset,username)

  a <- matrix(0,nrow(gyration),2)
  rownames(a) <- rownames(gyration)
  colnames(a) <- c("af","au")

  for (i in 1:nrow(a)) {
    a[i,1] <- (0.5*(as.numeric(lambda[i,"lambda_y"])^2+as.numeric(lambda[i,"lambda_z"])^2) - as.numeric(lambda[i,"lambda_x"])^2) / as.numeric(gyration[i,"Rf"])
#    a[i,1] <- a[i,1] * (4*pi*as.numeric(gyration[i,"Rf"])^2)
    a[i,2] <- 2 * pi * as.numeric(gyration[i,"Rg"])^2
  }
  
  return(a)
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

  mat <- -log(mat)
  return(mat)
}

#Function to calculate the energy for folded protein
getEnergyforFoldProtein <- function(PDBID,psurf,groDist,Xmatrix,thetaF) {
  sum <- 0
  for (i in 1:20) {
    for (j in 1:20) {
      sum <- sum + psurf[PDBID,i] * Xmatrix[colnames(psurf)[i],colnames(groDist)[j]] * groDist[j]
    }
  }
  sum <- sum * thetaF
  return(sum)
}

#Function to calculate the energy for unfolded protein
getEnergyforUnfoldProtein <- function(PDBID,punfold,groDist,Xmatrix,thetaF,gyration) {
  sum <- 0
  for (i in 1:20) {
    for (j in 1:20) {
      sum <- sum + punfold[PDBID,i] * Xmatrix[colnames(punfold)[i],colnames(groDist)[j]] * groDist[j]
    }
  }
  sum <- sum * as.numeric(gyration[PDBID,"Rf"])^3 / as.numeric(gyration[PDBID,"Rg"])^3 * thetaF
  return(sum)
}


#Below are used for actuall program running
#ddG <- energyCycle("ecoli40","wenjunh",split=TRUE)    

##Main excecuting code of this script, energycycle also used by Python code 
#load("pidsecoli.txt")
#load("pidsassist.txt")
ddG1 <- proteinEnergyCycle("wenjunh")#, pidsecoli=pidsecoli, pidsassist=pidsassist)

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





