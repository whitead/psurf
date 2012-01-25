source("SQLShareLib.R")

surfCutoff <- 0.3
kuhnLength <- 1.83
surfaceRoughness <- 3.5
sawExponent <- 3/5.
closeConfinementExponent <- 3.25
openConfinementExponent <- 5/3. # this number comes frmo De gennes on page 49 , not sure why but I used to have 1.5
closeCharaLength <- 47.9    #These numbers are subject to change
openCharaLength <- 36.3

##Obtain the raw counts of GroEL inside surface residues, both open and close form
sql <- paste("select * FROM [whitead@washington.edu].[GroEL_counts.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))


##This is the main code in this script for free energy model at 903 individual protein level, Modified Dec 19th, 2011
proteinEnergyCycle <- function(username, groForm, groDistribution=NULL, derivative=FALSE, dataset="ecoli_nogaps", contactDataset="ecoli", contactUsername="wenjunh") {


  #set up the exponents and distribution
  groDistForm <- "GroEL_Close"
  confinementExponent <- closeConfinementExponent
  charaLength <- closeCharaLength
  
  if(groForm == "Open") {
    groDistForm <- "GroEL_Open"
    confinementExponent <- openConfinementExponent
    charaLength <- openCharaLength
  }

  if(is.null(groDistribution)) {
    groDistribution <- groDist[groDistForm,]
  }

  #get the variables needed for the equations
  gyration <- getGyrationRadii(dataset, username)

  timeFraction <- getContactFractionTime(gyration, surfaceRoughness, charaLength)
  
  shapeFactors <- getShapeFactor(gyration)
  
  interactionMatrix <- getInteractionEnergy(contactDataset, contactUsername)


  #get the residue fractions
  cat("Fetching Data...")
  
  cutoff <- surfCutoff  #Set the surface cutoff
  psurf <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)


  resNumbersFold <- apply(psurf, MARGIN=1, sum)
  psurf <- t(apply(psurf, MARGIN=1, function(x) {x / sum(x)}))
  
  cutoff <- -1.0  #Set the unfolded cutoff
  punfold <- fetchAllSurfResidues(dataset, cutoff, normalize=FALSE, username)

  resNumbersUnfold <- apply(punfold, MARGIN=1, sum)
  punfold <- t(apply(punfold, 1, function(x) { x / sum(x)}))
 

  cat(" Done!\n")

  
  cat("Processing Dataset...")

  #calculate the energies from the distributions
  energy <- empty.df(rnames=rownames(gyration), cnames=c("Efold","Eunfold"))
  
  for(i in 1:nrow(psurf)) {
    cat(paste("\rProcessing Dataset...", i,"/", nrow(psurf)))
    
    PDBID <- rownames(gyration)[i]
    
    energy[i,"Efold"] <-
      getContactEnergy(psurf[PDBID,],
                       groDistribution,
                       interactionMatrix,
                       resNumbersFold[PDBID])
    energy[i,"Eunfold"] <-
      getContactEnergy(punfold[PDBID, ],
                       groDistribution,
                       interactionMatrix,
                       resNumbersUnfold[PDBID] * gyration[i,"Rf"] / gyration[i, "Rg"])
  }
  cat("\n")

  
  #calculate the final DDG vector
  ddG <- rep(0, nrow(energy))
  names(ddG) <- rownames(energy)
  
  for (i in 1:length(ddG)) {
    ddG[i] <-  timeFraction[i,"Tf"] * shapeFactors[i,"Af"] * energy[i,"Efold"] - timeFraction[i,"Tu"] * shapeFactors[i,"Au"] * energy[i,"Eunfold"] - (gyration[i,"Rg"] / charaLength)^confinementExponent
  }
 #put the various values into a dataframe
  allData <- data.frame(ddG=ddG, Tu=timeFraction[,"Tu"], resNumber=gyration[,"resNumber"], Au=shapeFactors[,"Au"], Eu=energy[,"Eunfold"], Tf=timeFraction[,"Tf"], Af=shapeFactors[,"Af"], Ef=energy[,"Efold"], Rg=gyration[,"Rg"], Rf=gyration[,"Rf"])
  return(allData)
} 

#get the general gyration radius (Rg = N^nv*l) (Rf = sqrt(lambda_x^2+lambda_y^2+lambda_z^2))
getGyrationRadii <- function(dataset, username) {

  #get data from SQLShare
  resNumbers <- fetchResNum(dataset, username)
  #get principal components
  lambda <- fetchGyration(dataset,username)
  #make data frame with two columns
  gyrationRadii <- empty.df(cnames=c("resNumber", "Rf", "Rg", "lambda.x", "lambda.y", "lambda.z"), rnames=rownames(lambda))
  
  for (i in 1:nrow(gyrationRadii)) {
    gyrationRadii[i, "resNumber"] <- resNumbers[i]
    gyrationRadii[i, "Rg"] <- (resNumbers[i]^sawExponent)*kuhnLength #equation for Rg
    #convert the principal components into numbers
    gyrationRadii[i, 4:6] <- sapply(lambda[i,], sqrt)
    #calculate the empirical radius of gyration from the PCs.
    gyrationRadii[i, "Rf"] <- sqrt(sum(gyrationRadii[i, 4:6]^2))

  }

  
  return(gyrationRadii)
}

#Get the Tao value for both fold and unfolded states
getContactFractionTime <- function(gyration, surfRough, charaLength){
  Tao <- empty.df(rnames=rownames(gyration), cnames=c("Tf", "Tu"))
  
  for (i in 1:nrow(Tao)) {


    if(gyration[i, "Rf"] + surfRough >= charaLength) {
      Tao[i, "Tf"] <- 1
    } else {
      Tao[i,"Tf"] <- 1 - (1 - surfRough / (charaLength - (gyration[i,"Rf"])))^3
    }


    if(gyration[i, "Rg"] + surfRough >= charaLength) {
      Tao[i, "Tu"] <- 1
    }  else {
      Tao[i,"Tu"] <- 1 - (1 - surfRough / (charaLength - (gyration[i,"Rg"])))^3
    }
  }

  
  return(Tao)
}

#get the Area of Contact for both fold and unfolded state
getShapeFactor <- function(gyration) {

  areas <- empty.df(rnames=rownames(gyration), cnames=c("Af", "Au"))
  
  for (i in 1:nrow(areas)) {
    #calculate the shape parameter, similar to asphericity.x
    a <- (0.5 * (gyration[i, "lambda.x"] ^ 2 + gyration[i, "lambda.y"] ^ 2) - gyration[i, "lambda.z"] ^ 2) / gyration[i, "Rf"] ^ 2
    areas[i,"Af"] <- a
    areas[i, "Au"] <- 0.5
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
getContactEnergy <- function(proteinDist,groDist,interactionMatrix,resNumber) {
  sum <- (proteinDist %*% interactionMatrix) %*% t(groDist)
  sum <- sum * resNumber
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

#Below are used for actual program
ddG1 <- proteinEnergyCycle("whitead", "Close", derivative=FALSE)
ddG2 <- proteinEnergyCycle("whitead", "Open", derivative=FALSE)


ddG1.assist <- proteinEnergyCycle("whitead", "Close", derivative=FALSE, dataset="assist_nogaps")
ddG2.assist <- proteinEnergyCycle("whitead", "Open", derivative=FALSE, dataset="assist_nogaps")

#Remove all proteins that should not fit within the cavity
zddG1.trunc <- ddG1[ddG1$Rf < closeCharaLength,]
ddG2.trunc <- ddG2[ddG2$Rf < closeCharaLength,]

#make plot of two plots

#setup color gradient, just in case
colorNumber <- 100
cuts <- cut(ddG1.trunc$resNumber, colorNumber, labels=F)
colorGrad <- colorRampPalette(c("black", "blue", "purple", "red"), space="Lab")(colorNumber)[cuts]

cairo_pdf("entropy_enthalpy_groel.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
ddG1.entropy <- sapply(ddG1.trunc$Rg, FUN=function(x) { (x / closeCharaLength)^closeConfinementExponent})
plot(0,0, xlab=expression(paste(-T * Delta * Delta * S / kT)), ylab=expression(paste(Delta * Delta * U / kT)), xlim=c(-50,0), ylim=c(-50, 0), col="white")
abline(v=0, lty=1, col="light gray", lwd=0.5)
abline(h=0, lty=1, col="light gray", lwd=0.5)
lines(seq(-150,50), seq(-150, 50), col="light gray", lty=2, xlim=c(-50,0), ylim=c(-50, 0))
points(-ddG1.entropy, ddG1.trunc$ddG + ddG1.entropy, xlim=c(-50,0), ylim=c(-50, 0), cex=0.75, lwd=0.5, col="gray25")
graphics.off()

cairo_pdf("open_close_entropy_enthalpy.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
ddG1.entropy <- sapply(ddG1.trunc$Rg, FUN=function(x) { (x / closeCharaLength)^closeConfinementExponent})
ddG2.entropy <- sapply(ddG2.trunc$Rg, FUN=function(x) { (x / openCharaLength)^openConfinementExponent})
plot(0, 0, xlab=expression(paste(-T * Delta * Delta * Delta * S / kT)), ylab=expression(paste(Delta * Delta * Delta * U / kT)), xlim=c(-20,5), ylim=c(-50, 40), col="white")
abline(v=0, lty=1, col="light gray", lwd=0.5, xlim=c(-20, 5))
abline(h=0, lty=1, col="light gray", lwd=0.5, xlim=c(-20, 5))
points(ddG2.entropy - ddG1.entropy, (ddG1.trunc$ddG - ddG2.trunc$ddG) + (ddG1.entropy - ddG2.entropy), xlim=c(-20,5), ylim=c(-50, 40), cex=0.75, lwd=0.5, col="gray25")
graphics.off()

cairo_pdf("open_close_1.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
plot(0,0, xlab=expression(paste("Close ", Delta * Delta * A, " [kT]")), ylab=expression(paste("Open ", Delta * Delta * A, " [kT]")), xlim=c(-100,25), ylim=c(-50, 50), pch=0, type="p", col="white")
lines(seq(-150,50), seq(-150, 50), col="light gray", lty=2, xlim=c(-100,25), ylim=c(-100, 25))
abline(v=0, lty=1, col="green", lwd=0.5)
abline(h=0, lty=1, col="blue", lwd=0.5)
points(ddG1.trunc$ddG, ddG2.trunc$ddG, cex=0.75, lwd=0.5, col="gray25", xlim=c(-100, 25), ylim=c(-100, 25))
points(ddG1.assist$ddG, ddG2.assist$ddG, xlim=c(-100, 25), ylim=c(-100, 25), pch=19, cex=0.75, lwd=0.5, col="red")
graphics.off()

print(as.double(sum(ddG1.trunc$ddG < ddG2.trunc$ddG)) / nrow(ddG1.trunc))

cairo_pdf("open_close_2.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
hist(ddG1.trunc$ddG - ddG2.trunc$ddG, xlab=expression(paste(Delta * Delta * Delta * A, " [KT]")), main="", border="black", col="dark gray")
graphics.off()

print(paste("Median close:", median(ddG1.trunc$ddG)))
print(paste("Median diff:", median(ddG1.trunc$ddG - ddG2.trunc$ddG)))

cairo_pdf("open_close_3.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
plot(ddG1.trunc$Rf, ddG2.entropy - ddG1.entropy, xlab=expression(paste(R[g], "[Å]")), ylab=expression(paste(Delta * Delta * Delta * A, " [kT]")), cex=0.75, lwd=0.5, col="red", ylim=range(ddG1.trunc$ddG - ddG2.trunc$ddG))

points(ddG1.trunc$Rf, ddG1.trunc$ddG - ddG2.trunc$ddG, lwd=0.5, col="gray25",cex=0.75)
graphics.off()

cairo_pdf("open_close_4.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
plot(ddG1.trunc$Rf, ddG1.trunc$ddG, xlab=expression(paste(R[g], "[Å]")), ylab=expression(paste("Close ", Delta * Delta * A, " [kT]")), cex=0.75, lwd=0.5, col=colorGrad, ylim=range(ddG1$ddG))


cairo_pdf("open_close_5.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
plot(ddG2$Rf, ddG2$ddG, xlab=expression(paste(R[g], "[Å]")), ylab=expression(paste("Open ", Delta * Delta * A, " [kT]")), cex=0.75, lwd=0.5, col="dark gray", ylim=range(ddG2$ddG))
graphics.off()

#Same as open_close_1 with color gradient for number of residues
cairo_pdf("open_close_6.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.65, fg="dark gray")
plot(0,0, xlab=expression(paste("Close ", Delta * Delta * A, " [kT]")), ylab=expression(paste("Open ", Delta * Delta * A, " [kT]")), xlim=c(-100,25), ylim=c(-50, 50), pch=0, type="p", col="white")
lines(seq(-150,50), seq(-150, 50), col="light gray", lty=2, xlim=c(-100,25), ylim=c(-100, 25))
abline(v=0, lty=1, col="blue", lwd=0.5)
abline(h=0, lty=1, col="green", lwd=0.5)
points(ddG1.trunc$ddG, ddG2.trunc$ddG, cex=0.75, lwd=0.5, xlim=c(-100, 25), ylim=c(-100, 25), col=colorGrad)
points(ddG1.assist$ddG, ddG2.assist$ddG, xlim=c(-100, 25), ylim=c(-100, 25), pch=19, cex=0.75, lwd=0.5, col="orange")
graphics.off()



#Now, find the optimal GroEL distribution
propDists <- matrix(rep(0,20**2), nrow=20)
propDists[sapply(1:20, FUN=function(x) {x + (x - 1) * 20})] <- 1
propDists <- data.frame(propDists)
colnames(propDists) <- colnames(groDist)

propddG <- rep(0,20)
names(propddG) <- colnames(propDists)

for(i in 1:20) {
 
  ddG <- proteinEnergyCycle("whitead", "Close", groDistribution=propDists[i,])
  propddG[i] <- median(ddG$ddG)
  
}

print(propddG)

cairo_pdf("optim_groeol.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.55, fg="dark gray")
barplot(propddG, names.arg=aalist.sh, border="black", xlab="Residue", ylab=expression(paste(Delta * Delta * A, "[kT]")))
graphics.off()


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
