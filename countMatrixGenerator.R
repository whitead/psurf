#! /usr/bin/env Rscript

#job <- commandArgs()[6]

empty.df<- function(cnames, rnames, default=NA){
  
  df<-data.frame(matrix(rep(default,length(cnames)*length(rnames)), nrow=length(rnames)))
  colnames(df)<-cnames
  rownames(df) <- rnames
  return(df)
}

sampleCounts <- function(data,anames,turnOffC=FALSE) {

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

matrixGenerator <- function(lst) {
  data <- read.table(lst, header=T)
  anames <- colnames(data[-c(1,2)])
  rnum <- length(anames)

  labels.ord <- c("K", "R", "H", "E", "D", "N", "Q", "G", "S", "T", "A", "I", "L", "M", "P", "V", "C", "F", "W", "Y")
  labels.one <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  labels.ord.f <- c("K", "R", "H", "E", "D", "N", "Q", "G", "S", "T", "A", "I", "L", "M", "P", "V", "C", "F", "W", "Y", "*")
  labels.one.f <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")

  bootstrap <- 1
  interaction <- data.frame(matrix(rep(0,bootstrap * length(labels.one)), nrow=bootstrap))
  colnames(interaction) <- labels.one

  water <- data.frame(matrix(rep(0,bootstrap * length(labels.one)), nrow=bootstrap))
  colnames(water) <- labels.one

  for(i in 1:bootstrap) {
    counts <- sampleCounts(data, anames)
    interaction[i,] <- 1 - counts["FREE",] / counts["TOTAL",]
    water[i,] <- counts["WATER",] / apply(counts[c(anames, "WATER"),,], 2, sum)
  }


  counts <- counts[-c(21,22,24),]

  for (i in 1:ncol(counts)) {
    counts[,i] <- counts[,i] / sum(counts[,i])
  }
  return(counts)
}

