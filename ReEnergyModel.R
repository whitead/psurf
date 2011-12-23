library(alabama)
library(KernSmooth)
source("SQLShareLib.R")

surfCutoff <- 0.3
kuhnLength <- 1.83
surfaceRoughness <- 3.5
sawExponent <- 3/5.
confinementExponent <- 3.25
ellipExponent <- 1.6
charaLength <- 92.6    #These numbers are subject to change


##Obtain the raw counts of GroEL inside surface residues, both open and close form
sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))


##This is the main code in this script for free energy model at 903 individual protein level, Modified Dec 19th, 2011
proteinEnergyCycle <- function(username, groDistForm, derivative=FALSE, dataset1="ecoli", dataset2="assist") {
  
  gyration <- getGyrationRadii(dataset1, username)

  timeFraction <- getContactFractionTime(gyration, surfaceRoughness, charaLength)
  
  contactArea <- getContactArea(gyration, dataset1, username)
  
  interactionMatrix <- getInteractionEnergy(dataset1, username)

  surfDensities <- getSurfaceDensities(dataset1, username)
  cat("Fetching Data...")
  
  cutoff <- surfCutoff  #Set the surface cutoff
  psurf <- fetchAllSurfResidues(dataset1, cutoff, normalize=TRUE, username)
  psurfassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=TRUE, username)
  
  cutoff <- -1.0  #Set the unfolded cutoff
  punfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=TRUE, username)
  punfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=TRUE, username)

  cat(" Done!\n")

  if (derivative == FALSE) {
    #Generate ddG matrix
    cat("Processing DataSet...\n")
    
    energy <- empty.df(rnames=rownames(gyration), cnames=c("Efold","Eunfold"))
    
    for(i in 1:nrow(psurf)) {
      cat(paste("\r", i,"/", nrow(psurf)))
      
      PDBID <- rownames(psurf)[i]
      
      energy[i,"Efold"] <-
        getContactEnergy(psurf[PDBID,],
                         groDist[groDistForm,],
                         interactionMatrix,
                         surfDensities[PDBID])
    
      energy[i,"Eunfold"] <-
        getContactEnergy(punfold[PDBID, ],
                         groDist[groDistForm,],
                         interactionMatrix,
                         surfDensities[PDBID] * (gyration[PDBID, "Rf"] / gyration[PDBID, "Rg"]) ^ 3)
    }
    
    #calculate the final DDG matrix
    ddG <- rep(0, nrow(energy))
    names(ddG) <- rownames(energy)
    
    for (i in 1:length(ddG)) {
      ddG[i] <-  -timeFraction[i,"Tu"] * contactArea[i,"Au"] * energy[i,"Eunfold"] + timeFraction[i,"Tf"] * contactArea[i,"Af"] * energy[i,"Efold"] + (gyration[i,"Rg"] / charaLength)^confinementExponent
    }
    return(ddG)
  } else {
    #method to get the derivative of the ddG and make plot of that
    cat("Calculating Derivative...")

    proteinResDev <- matrix(0,nrow(psurf),ncol(psurf)) #used matrix here because we will do matrix multiplication on that later
    rownames(proteinResDev) <- rownames(psurf)
    colnames(proteinResDev) <- colnames(psurf)

    for (i in 1:nrow(proteinResDev)) {
      #calculate the prefactor of the energy matrix
      temp <- -timeFraction[i,"Tu"] * contactArea[i,"Au"] * punfold[i,] + timeFraction[,"Tf"] * contactArea[i,"Af"] * psurf[i,]  
      proteinResDev[i,] <- as.matrix(temp) %*% interactionMatrix
      }

    resDev <- apply(proteinResDev, MARGIN=2, FUN=mean)

    cat("Done! \n")
    
    #make plot
    cairo_pdf('Ecoli_ResDiffDevFull.pdf', width=8, height=5)
    par(family='LMSans10', cex.axis=0.65, ps=11)
    barplot(resDev, col='gray')#, ylim=c(-1.0,2.0))
    graphics.off()
  }
} 

#get the general gyration radius (Rg = N^nv*l) (Rf = sqrt(lambda_x^2+lambda_y^2+lambda_z^2))
getGyrationRadii <- function(dataset, username) {

  #get data from SQLShare
  resNumbers <- fetchResNum(dataset, username)
  #get principal components
  lambda <- fetchGyration(dataset,username)
  #make data frame with two columns
  gyrationRadii <- empty.df(cnames=c("Rf", "Rg", "lambda.x", "lambda.y", "lambda.z"), rnames=rownames(lambda))
  
  for (i in 1:nrow(gyrationRadii)) {
    gyrationRadii[i, "Rg"] <- (resNumbers[i]^sawExponent)*kuhnLength #equation for Rg
    #convert the principal components into numbers
    gyrationRadii[i, 3:5] <- sapply(lambda[i,], sqrt)
    #calculate the empirical radius of gyration from the PCs.
    gyrationRadii[i, "Rf"] <- sqrt(sum(gyrationRadii[i, 3:5]^2))

  }

  return(gyrationRadii)
}

#Get the Tao value for both fold and unfolded states
getContactFractionTime <- function(gyration, surfRough, charaLength){
  Tao <- empty.df(rnames=rownames(gyration), cnames=c("Tf", "Tu"))
  
  for (i in 1:nrow(Tao)) {
      Tao[i,"Tf"] <- 1 - (1 - surfRough / (charaLength - as.numeric(gyration[i,"Rg"])))^3
      Tao[i,"Tu"] <- 1 - (1 - surfRough / (charaLength - as.numeric(gyration[i,"Rf"])))^3
    }

  return(Tao)
}

#get the Area of Contact for both fold and unfolded state
getContactArea <- function(gyration,dataset,username) {

  #get empirical areas
  surfaceAreas <- fetchChargeAndSA(dataset, username)$surface_area
  areas <- empty.df(rnames=rownames(gyration), cnames=c("Af", "Au"))
  
  for (i in 1:nrow(areas)) {
    #calculate the shape parameter, similar to asphericity.x
    a <- (0.5 * (gyration[i, "lambda.x"] ^ 2 + gyration[i, "lambda.y"] ^ 2)
                               - gyration[i, "lambda.z"] ^ 2) / gyration[i, "Rf"]
    areas[i,"Af"] <- a * surfaceAreas[i]
    ellip <- 4 * pi * (((gyration[i, "lambda.x"]*gyration[i, "lambda.y"])^ellipExponent + (gyration[i, "lambda.x"]*gyration[i, "lambda.z"])^ellipExponent + (gyration[i, "lambda.y"]*gyration[i, "lambda.z"])^ellipExponent) / 3)^(1 / ellipExponent)
    ratio <- surfaceAreas[i] / ellip
    areas[i, "Au"] <- ratio * 0.5 * 4 * pi * gyration[i, "Rg"] ^ 2 #see equation in paper  #added the ratio term to adjust the magnitude of Au
  }

  
  return(areas)
}
  

#calucate the surface density for all proteins in a given dataset
getSurfaceDensities <- function(dataset, username=NULL) {

  surfResidues <- fetchAllSurfResidues(dataset, cutoff, username=username, normalize=FALSE)
  if(is.null(username)) {
    chargeAndSA <- fetchChargeAndSA(dataset, username)
  } else {
    chargeAndSA <- fetchChargeAndSA(dataset, username)
  }

  densities <- rep(0, nrow(surfResidues))
  names(densities) <- rownames(surfResidues)
  for(i in 1:length(densities)) {
    densities[i] <- sum(surfResidues[i,]) / chargeAndSA$surface_area[i]
  }

  return(densities)
}

#get the interaction energies between residue types
getInteractionEnergy <- function(dataset, username=NULL, glycine=FALSE) {

  dataset <- paste(paste(dataset,"total","contacts", sep="_"), "csv" ,sep=".")
 
  #get the data
  if(is.null(username)) {
    countMatrix <- fetchContacts(dataset)
  }
  else{
    countMatrix <- fetchContacts(dataset, username)
  }

  #re-order it
  anames <- colnames(countMatrix[-c(1,2)])
  aanum <- length(anames)
  anames.ord <- order(anames)

  countMatrix <- countMatrix[c(1,2,anames.ord + 2)]

  #sum it
  contactMatrix <- sampleContacts(countMatrix, random=FALSE)
  contactMatrix[1:aanum,] <- contactMatrix[order(rownames(contactMatrix)[1:aanum]),]
  rownames(contactMatrix) <- c(sort(rownames(contactMatrix)[1:aanum]), rownames(contactMatrix)[-(1:aanum)])
  contactMatrix[1:aanum, 1:aanum] <- contactMatrix[1:aanum,1:aanum] + t(contactMatrix[1:aanum, 1:aanum])



  
  #normalize it to the effects remove amounts of each amino acid
  #add the effect of the free residues
  for(i in 1:aanum) {
    contactMatrix[1:aanum, i] <- contactMatrix[1:aanum, i] / sum(contactMatrix[1:aanum, i])
    contactMatrix[1:aanum, i] <- contactMatrix[1:aanum, i] * (1 - contactMatrix["FREE", i] / contactMatrix["TOTAL", i])
    contactMatrix["FREE", i] <- contactMatrix["FREE", i] / contactMatrix["TOTAL", i]
  }

  

  #normalize it so all events sum to 1
  normRows <- c(1:aanum, which(rownames(contactMatrix) == "FREE"))
  contactMatrix[normRows, 1:aanum] <- contactMatrix[normRows, 1:aanum] / sum(contactMatrix[normRows, 1:aanum])
  
  #make it relative to being a free residue
  mat <- matrix(rep(0, aanum**2), nrow=aanum)
  rownames(mat) <- sort(anames)
  colnames(mat) <- sort(anames)
  for(i in 1:aanum) {
    for(j in 1:aanum) {
      mat[i,j] <- contactMatrix[i,j] * contactMatrix[j,i] / (contactMatrix["FREE",i] * contactMatrix["FREE", j])
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

#Function to calculate the energy for folded protein
getContactEnergy <- function(proteinDist,groDist,interactionMatrix,surfDensity) {
  sum <- 0
  for (i in 1:20) {
    for (j in 1:20) {
      sum <- sum + proteinDist[i] * interactionMatrix[colnames(proteinDist)[i],colnames(groDist)[j]] * groDist[j]
    }
  }
  sum <- sum * surfDensity
  return(sum)
}

#Method to calculate the first derivative of Residue and plot the bar graph
getSurfResFirstDev <- function(username, dataset) {
  interactionMatrix <- getInteractionEnergy(dataset, username)

  cat("Processing Data...")
  cutoff <- surfCutoff
  psurf <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)

## Code for Andrew's special demand 
  cutoff <- -1.0
  punfold <- fetchAllSurfResidues(dataset, cutoff, normalize=TRUE, username)

  meanFold <- apply(psurf, MARGIN=2, FUN=mean)
  meanUnfold <- apply(punfold, MARGIN=2, FUN=mean)

  meanDiff <- meanFold - meanUnfold
  cat("Done!\n")
#######
  
  #plot 1
  cairo_pdf('EColi_ResDiffDev.pdf',width=8, height=5)
  par(family='LMSans10', cex.axis=0.65, ps=11)
  barplot(meanDiff, col='gray', ylim=c(-1.0,2.0))
  graphics.off()
  
  ecoliDist <- (meanFold - meanUnfold) %*% interactionMatrix

  cairo_pdf('EColi_SpecialResFirstDev.pdf',width=8, height=5)
  par(family='LMSans10', cex.axis=0.65, ps=11)
  barplot(ecoliDist, col='gray', ylim=c(-1.0,2.0))
  graphics.off()

  cat(paste(dataset,"_ResFirstDev.pdf has been added to your current directory", sep=""))
  cat(" \n")

  #plot 2 : GroEL Open residue times interaction Matrix
  groOpen <- as.numeric(groDist["GroEL_Open",]) %*% interactionMatrix

  cairo_pdf('groOpen_ResFirstDev.pdf',width=8, height=5)
  par(family='LMSans10', cex.axis=0.65, ps=11)
  barplot(groOpen, col='gray', ylim=c(-1.0,2.0))
  dev.off()

  cat("groOpen_ResFirstDev.pdf has been added to your current directory\n")

  #plot 3 : GroEL Close residue times interaction Matrix
  groClose <- as.numeric(groDist["GroEL_Close",]) %*% interactionMatrix

  cairo_pdf('groClose_ResFirstDev.pdf',width=8, height=5)
  par(family='LMSans10', cex.axis=0.65, ps=11)
  barplot(groClose, col='gray', ylim=c(-1.0,2.0))
  dev.off()

  cat("groClose_ResFirstDev.pdf has been added to your current directory\n")
}



#Below are used for actuall program running
ddG1 <- proteinEnergyCycle("wenjunh", "GroEL_Close", derivative=FALSE)
ddG2 <- proteinEnergyCycle("wenjunh", "GroEL_Open")

#getSurfResFirstDev("wenjunh","ecoli")

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





