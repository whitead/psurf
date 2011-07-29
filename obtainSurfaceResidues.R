#!/usr/bin/env Rscript

#source("/home/whitead/Documents/ProteinSurfaces/circlecorr.R")
library(RCurl)

username <- "wenjunh"
apikey <- "670e3ae6fbd5bd9efb55e47bc5b908b8"

#cutoff <- 0.0
aalist = read.table("/home/wenjunh/Documents/ProteinSurfaces/aalist.txt", header=T)$res
aalist.sh = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
labels.sh.ord <- c("K", "R", "H", "E", "D", "N", "Q", "G", "S", "T", "A", "I", "L", "M", "P", "V", "C", "F", "W", "Y")
fetchdata <- function(sql) {
  host = "https://sqlshare-rest.cloudapp.net"
  selector = "/REST.svc/v1/db/file"
  query=paste("?sql=",URLencode(sql),sep="")

  data = getURL(paste(host,selector,query,sep=""), httpheader=c(Authorization =paste("ss_apikey ", username, "@washington.edu :", apikey, sep="")), verbose = FALSE, ssl.verifypeer = FALSE, ssl.verifyhost = FALSE)

  splitrow <- function(row) strsplit(row, "\t")

  rdata = lapply(strsplit(data,"\r\n"),splitrow)

  return(rdata)
}

fetchPDBIDs <- function(dataset) {

  sql = paste("select pdb_id FROM [", username, "@washington.edu].[",dataset,"_1.csv]",sep="")
  idlist <- fetchdata(sql)
  #minus one to skip the column headers
  ids <- rep("", length(idlist[[1]]) - 1)
  for(i in 2:length(idlist[[1]])) {
    ids[i - 1] <- as.character(idlist[[1]][i])
  }
  return(ids)
}

fetchSurfResidues <- function(dataset, cutoff, pdbid) {
  sql = paste("SELECT res_type FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE pdb_id =\'", pdbid, "\' AND res_surface_area_ratio > ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")

  rlist <- fetchdata(sql)
  rs <- rep("", length(rlist[[1]]) - 1)
  for(i in 2:length(rlist[[1]])) {
    rs[i - 1] <- as.character(rlist[[1]][i])
  }
  return(rs)
}

fetchAllSurfResidues <- function(dataset, cutoff, normalize=FALSE) {
  sql = paste("SELECT pdb_id, res_type FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE res_surface_area_ratio > ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")

  rlist <- fetchdata(sql)
  #print(rlist)

  sql = paste("SELECT DISTINCT pdb_id FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE res_surface_area_ratio > ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")
  
  plist <- fetchdata(sql)
  #print(plist)
  pdbids <- rep("", length(plist[[1]]) - 1)
  for(i in 2:length(plist[[1]])) {
    pdbids[i - 1] <- plist[[1]][[i]][1]
  }
  
  data <- empty.df(aalist, pdbids, default=0)
  curid <- NULL
  i <- 2
  while(i <= length(rlist[[1]])) {

    curid <- rlist[[1]][[i]][1]
    pastid <- curid
    rs <- c()
    while(pastid == curid) {
      rs <- c(rs, rlist[[1]][[i]][2])
      if(i >= length(rlist[[1]])) {
        i <- i + 1
        break
      }
      curid <- rlist[[1]][[i + 1]][1]
      i <- i + 1
    }
    rs.table <- table(rs)
    for(j in 1:length(rs.table)) {
      if(normalize) {
        rs.table <- rs.table / sum(rs.table)
      } 
      data[curid,names(rs.table)[j]] <- rs.table[j]
    }
  }
  
  return(data)
}

empty.df<- function(cnames, rnames, default=NA){
        df<-data.frame(matrix(rep(default,length(cnames)*length(rnames)), nrow=length(rnames)))
        colnames(df)<-cnames
        rownames(df) <- rnames
        return(df)
}

