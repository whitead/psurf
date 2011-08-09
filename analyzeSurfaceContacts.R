#! /usr/bin/env Rscript

sourcex("SQLShareLib.R")

job <- commandArgs()[6]
data <- fetchFrame(job)
anames <- colnames(data[-c(1,2)])
rnum <- length(anames)

matPlot <- function(x) {

  rownames(x) <- aalist.one
  colnames(x) <- c(aalist.one, "*")

  return(as.matrix(x[labels.ord, labels.ord.f]))

}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

calculateBackground <- function(counts) {
  bpairs <- empty.df(c(anames, "FREE"), anames, default=0)
  for(i in 1:(rnum)) {
    for(j in 1:rnum) {
      bpairs[i,j] <- as.double(sum(counts[j,])) / sum(counts[anames,]) * (counts["TOTAL", i] - counts["FREE", i]) / counts["TOTAL",i]
    }
    bpairs[i,j + 1] <- as.double(counts["FREE", i]) / counts["TOTAL",i]
  }

  return(bpairs)
}

calculatePairs <- function(counts) {
  ppairs <- empty.df(c(anames, "FREE"), anames, default=0)

  for(i in 1:rnum) {
    for(j in 1:rnum) {
      ppairs[i,j] <- as.double(counts[i,j]) / sum(counts[i,]) * (counts["TOTAL", i] - counts["FREE", i]) / counts["TOTAL",i]
    }
    ppairs[i,rnum + 1] <- counts["FREE", i] / counts["TOTAL", i]
  }

  return(ppairs)
}

empty.df<- function(cnames, rnames, default=NA){
  
  df<-data.frame(matrix(rep(default,length(cnames)*length(rnames)), nrow=length(rnames)))
  colnames(df)<-cnames
  rownames(df) <- rnames
  return(df)
}

sampleCounts <- function(turnOffC=FALSE) {

  ids <- unique(data[,"pdb_id"])
  indices <- sample(length(ids), replace=TRUE)

  counts <- empty.df(anames, data[data[,"pdb_id"] == ids[1],"res_type"], default=0)

  cat("\n")
  for(i in 1:length(indices)) {
    temp <- data[data[,"pdb_id"] == ids[indices[i]], anames]
    for(j in 1:length(temp)) {
      counts[j] <- counts[j] + temp[j]
    }
    cat(paste("\r", i,"/", length(indices)))
  }
  cat("\n")

  if(turnOffC) {
    counts["CYS", "CYS"] <- 0 
  }
  
  return(counts)

}

save.image()

bootstrap <- 250
interaction <- data.frame(matrix(rep(0,bootstrap * length(aalist.sh)), nrow=bootstrap))
colnames(interaction) <- aalist.sh

water <- data.frame(matrix(rep(0,bootstrap * length(aalist.sh)), nrow=bootstrap))
colnames(water) <- aalist.sh

for(i in 1:bootstrap) {
  counts <- sampleCounts()
  interaction[i,] <- 1 - counts["FREE",] / counts["TOTAL",]
  water[i,] <- counts["WATER",] / apply(counts[c(anames, "WATER"),,], 2, sum)
  cat("\r",i,"/",bootstrap)
}
cat("\n")
save.image()

wyy <- apply(water, 2, function(x) {quantile(x,c(0.025,0.5,0.975), na.rm=T)})
iyy <- apply(interaction, 2, function(x) {quantile(x,c(0.025,0.5,0.975), na.rm=F)})

gi <- which(colnames(wyy) == "G")
wyy <- wyy[,-gi]
iyy <- iyy[,-gi]

tiff("contact_plots.tiff", width=3.3*500,height=4.5*500, pointsize=8, res=500)
par(family="LMSans10", mfrow=c(2,1), cex.axis=0.65, mar=c(2,4,1,0.5))

xp <- barplot(iyy[2,], ylim=c(0,max(iyy[3,] + 0.03)), ylab="Proportion Interacting", col="light gray")
error.bar(xp, iyy[2,], lower=(iyy[2,] - iyy[1,]), upper=(iyy[3,] - iyy[2,]), length=0.03)

par(mar=c(3,4,2,0.5))

xp <- barplot(wyy[2,], ylim=c(0,max(wyy[3,] + 0.03)), ylab="Water per Contact", col="light gray")
error.bar(xp, wyy[2,], lower=(wyy[2,] - wyy[1,]), upper=(wyy[3,] - wyy[2,]), length=0.03)
graphics.off()

