library(alabama)
library(KernSmooth)
source("SQLShareLib_Wenjun.R")

sql <- paste("select * FROM [wenjunh@washington.edu].[groel_insurfres_count.csv]")
rawData <- fetchdata(sql)
groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist["GroEL_Close", ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroEL_Close", ]))
groDist["GroEL_Open", ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroEL_Open", ]))



energyCycle <- function(username, dataset1="ecoli",dataset2="assist", contacts=fetchContacts(paste(dataset1, "_surface_contacts.csv",sep=""))) {
  countMatrix <- sampleContacts(contacts)

  #  obtain surface residues for both E.Coli fold and E.Coli unfold
  pidsecoli <- fetchPDBIDs(dataset1, username)
  pidsassist <- fetchPDBIDs(dataset2, username)
  cat("Fetching Data...")
  
  cutoff <- 0.3
  csurf <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  csurfassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)
  
  cutoff <- -1.0
  cunfold <- fetchAllSurfResidues(dataset1, cutoff, normalize=FALSE, username)
  cunfoldassist <- fetchAllSurfResidues(dataset2, cutoff, normalize=FALSE, username)
  
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
  
#  pidsecoli <- c("1AF6","1AZO", "1BYI", "1D9E", "1EUW", "1FL2", "1I78", "1JW9", "1K6D", "1K6W","1K92", "1KMJ", "1KQF", "1LV7", "1M40", "1NI5", "1NUL", "1OTS", "1PW4", "1Q16","1Q7E", "1QD5", "1QJP", "1RC6", "1SCZ", "1T16", "1TLY", "1U5W", "1U7G", "1UJC","1UWF", "1VI3", "1WXI", "1YOE", "1YW6", "2ATE", "2AU7", "2B82", "2BZ1", "2CFQ","2D4U", "2E85", "2F1V", "2GFP", "2GRX", "2GUF", "2HDI", "2J1N", "2JF2", "2PAN","2PGX", "2QI9", "2QOM", "2VGD", "2VQI", "2VV5", "2VWS", "2WCD", "2WJ9", "2WJR","2X9K", "2XE3", "2XOV", "2Z98", "2ZCU", "2ZFG", "2ZHH", "3AEH", "3BBY", "3FV5","3GEA", "3GP6", "3H90", "3HFX", "3HO9", "3I87", "3IIQ", "3IP0", "3JQO", "3KCU","3L1L", "3LGI", "3NKA", "3O7Q", "3OHN", "3PIK", "3QE7", "3RFZ")
#  pidsecoli <- c("1AB4", "1AT1", "1B5T", "1BDF", "1BS0", "1DC3", "1DKG", "1E9I", "1EVL", "1GG1","1MXB", "1XEY", "2C4N", "2EHJ", "2GQR", "3LTI", "3PCO", "3Q9L")
#  pidsecoli <- c("1BDF","1DKG", "1EVL", "1XEY", "2C4N", "2EHJ", "2GQR", "3LTI", "3PCO", "3Q9L")
#  pidsecoli <- c("2PAN","2WCD")
  
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
      ener <- c(0,0)
 #     ener[1] <- newFreeEnergyModel(pidsecoli[i], countMatrix, groDist["GroEL_Close", ], pfracs, counts)
      ener[1] <- newFreeEnergyModel(pidsecoli[i], countMatrix, groDist["GroEL_Open", ], pfracs, counts)
      ddGecoli[i,1] <- -ener[1]# - ener[2]
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
    ener <- c(0,0)
#    ener[1] <- newFreeEnergyModel(pidsassist[i], countMatrix, groDist["GroEL_Close", ], pfracsassist, countsassist)
    ener[1] <- newFreeEnergyModel(pidsassist[i], countMatrix, groDist["GroEL_Open", ], pfracsassist, countsassist)
    ddGassist[i,1] <- -ener[1]# - ener[2]
    cat(paste("\rProcessing DataSet2...", i,"/",length(pidsassist)))
    ddGassist[i,2] <- sum(countsassist$cunfoldassist[pidsassist[i],])
  }
  cat(" \n")
  cat("Completed! \n")
  ddG <- list(ddGecoli, ddGassist)
  names(ddG) <- c("ecoli","assist")

  return(ddG)
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

newFreeEnergyModel <- function(PDBID,countMatrix,yDist,pfracs,counts) {
  countMatrix <- apply(countMatrix[c(1:20, which(rownames(countMatrix) == "WATER")),], 2, function(row) {row / sum(row)} )
  yDist <- yDist[match(colnames(countMatrix),names(yDist))]
#  print(yDist)
#  print(paste(pfracs$psurf[PDBID,],pfracs$pburied[PDBID,],pfracs$punfold[PDBID,]))
  
  ener <- c(0,0)
  names(ener) <- list("Gfold","Gunfold")
  sigmaXln <- array(0,2)
  sigmaX <- array(0,2)

  sum <- 0.
  for (i in 1:ncol(counts$csurf)) {
    pxy <- 0.
    for (j in 1:length(yDist)) {
      pxy <- pxy + countMatrix[i,j] * yDist[j]#pfracs$psurf[PDBID,j]#yDist[j]
    }
    sigmaY <- log(pxy)
    sum <- sum + as.double(counts$csurf[PDBID,i] * sigmaY)
  }
  Xsurf <- sum
#  print(Xsurf)

  sum <- 0.
  for (i in 1:ncol(counts$cburied)) {
    pxy <- 0.
    for (j in 1:length(pfracs$pburied)) {
      pxy <- pxy + countMatrix[i,j] * pfracs$pburied[PDBID,j]
    }
    sigmaY <- log(pxy)
    sum <- sum + as.double(counts$cburied[PDBID,i] * sigmaY)
  }
  Xburied <- sum
#  print(Xburied)
  ener[1] <- Xsurf + Xburied

  sum <- 0.
  for (i in 1:ncol(counts$cunfold)) {
    pxy <- 0.
    for (j in 1:length(pfracs$punfold)) {
      pxy <- pxy + countMatrix[i,j] * pfracs$punfold[PDBID,j]
    }
    sigmaY <- log(pxy)
    sum <- sum + as.double(counts$cunfold[PDBID,i] * sigmaY)
  }
  Xunfold <- sum
#  print(Xunfold)
  ener[2] <- Xunfold

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
    print(optMin)
#    optMax <- auglag(par = ini, fn = g, gr = devI, hin = constrIeq, hin.jac = constrIeq.jac, heq = constrEq, heq.jac = constrEq.jac, control.outer = list(trace = T, eqs = 1E-100, method = "BFGS"), control.optim = list(fnscale = -1))
#    print(optMax)
    }
#  print(opt[1])
#  print(opt[2])
#  print(opt[3])
  optimiz(ini, 1)
  print(groDist["GroEL_Close",])
  print(g(groDist["GroEL_Close",]))

#  print(devI(groDist["GroEL_Open",]))

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

#  xlim <- c(quantile(x[which(y<1000)],probs=c(0.025,0.975))[1],quantile(x[which(y<1000)],probs=c(0.025,0.975))[2])
#  ylim <- c(0,1000)
  xlim <- c(quantile(x,probs=c(0.025,0.975))[1],quantile(x,probs=c(0.025,0.975))[2])
  ylim <- c(quantile(y,probs=c(0.025,0.975))[1],quantile(y,probs=c(0.025,0.975))[2])
#  xlim <- c(-100,0)
#  ylim <- c(0,8000)

  mest <- bkde2D(x=cbind(x,y), bandwidth=c(5,0.05), gridsize=c(1000,1000), range.x=list(xlim, ylim))
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
  plot(density(x), xlab=xlab, col="red", type="l", lwd=4, main="Density Plot", xlim=c(-60,20),ylim=c(0,0.04))
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
  plot(x=seq(min(data1[,1]),max(data1[,1]),0.01),y=dnorm(seq(min(data1[,1]),max(data1[,1]),0.01), mean(data1[,1]), sqrt(var(data1[,1]))),type="l", lwd=4, col="red", xlab="Hydrophobicity",ylab="Density", ylim=c(0,0.04),xlim=c(-50,50))
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

ddG2 <- energyCycle(username="wenjunh")
#normalPlot(ddG,"ddGquan")
#plotPS(x=ddG$ecoli[,1],y=ddG$ecoli[,2],xpoints=ddG$assist[,1],ypoints=ddG$assist[,2], xlab="ddG",ylab="length",plotName1="trial1",plotName2="trial2")
#ellipsoidPlot(ddG,"ellipsoid")
#plotPS(x=data$ecoli[,1],y=data$ecoli[,2],xpoints=data$assist[,1],ypoints=data$assist[,2], xlab="phobicity",ylab="length",plotName1="trial1", plotName2="trial2")
#data <- phobicity("wenjunh")
#surfCharge <- netCharge("wenjunh") 
#plotPS(x=surfCharge$ecoli[,1],y=surfCharge$ecoli[,2],xpoints=surfCharge$assist[,1],ypoints=surfCharge$assist[,2], xlab="netCharge",ylab="length",plotName1="trial1", plotName2="trial2")
#normPlot(ddG$ecoli,ddG$assist,"trial")
#plotPS(x=surfCharge$ecoli[,1],y=data$ecoli[,1],xpoints=surfCharge$assist[,1],ypoints=data$assist[,1], xlab="netCharge",ylab="Hydrophobicity",plotName1="trial1", plotName2="trial2")
#energyCycle("assist", username="wenjunh", contacts=fetchContacts("ecoli_surface_contacts.csv", "wenjunh"))
#energyBootstrap(1000, "ecoli", username="wenjunh")
#minimizeEnergy("ecoli", "wenjunh")
