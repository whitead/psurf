#!/usr/bin/env Rscript

source("circlecorr.R")
library(RCurl)

myUsername <- "whitead"
apikey <- "58f90137316aedb538b85a54955173c0"

cutoff <- 0.3
aalist <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
aalist.sh <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
labels.sh.ord <- c("K", "R", "H", "E", "D", "N", "Q", "G", "S", "T", "A", "I", "L", "M", "P", "V", "C", "F", "W", "Y")

#Get raw data for the given query
fetchdata <- function(sql) {
  host <- "https://sqlshare-rest.cloudapp.net"
  selector <- "/REST.svc/v1/db/file"
  query<-paste("?sql=",URLencode(sql),sep="")

  data <- getURL(paste(host,selector,query,sep=""),
httpheader=c(Authorization =paste("ss_apikey ", myUsername,
"@washington.edu :", apikey, sep="")), verbose = FALSE, ssl.verifypeer
                 = FALSE, ssl.verifyhost=FALSE, .encoding="UTF-8")

  splitrow <- function(row) strsplit(row, "\t")

  rdata <- lapply(strsplit(data,"\r\n"),splitrow)

  return(rdata)
}

#Turn an SQL contacts table into a dataframe
fetchContacts <- function(tableName, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }

  sql <- paste("select * FROM [", username, "@washington.edu].[", tableName, "]", sep="")
  
  cat(paste("Fetching", tableName, "\n"))
  rawData <- fetchdata(sql)
  rnames <- (1:(length(rawData[[1]]) - 1))
  
  data <- empty.df(rawData[[1]][[1]], rnames)

  lconvert <- function(row) {sapply(row[-c(1,2)], as.integer)}

  temp <- lapply(rawData[[1]][-1], lconvert)
  data[,-c(1,2)] <- matrix(unlist(temp), ncol=ncol(data) - 2, byrow=T)

  lconvert <- function(row) {row[1:2]}
  temp <- lapply(rawData[[1]][-1], lconvert)
  
  data[,1:2] <- matrix(unlist(temp), ncol=2, byrow=T)
  
  return(data)

}

#Given a dataframe containing contact matrices for a list of pdb_ids,
#this will turn it into a single contact matrix using a random sample of the pdb_ids (with replacement)
sampleContacts <- function(countDataFrame, turnOffC=FALSE, random=TRUE) {

  ids <- unique(countDataFrame[,1])
  if(random == TRUE) {
    indices <- sample(length(ids), replace=TRUE)
  }else{
    indices <- 1:length(ids)
  }
  
  anames <- colnames(countDataFrame)[-c(1,2)]

  counts <- empty.df(anames, countDataFrame[countDataFrame[,1] == ids[1],"res_type"], default=0)

  cat("sampling...")
  for(i in 1:length(indices)) {
    temp <- countDataFrame[countDataFrame[,1] == ids[indices[i]], anames]
    for(j in 1:length(temp)) {
      counts[j] <- counts[j] + temp[j]
    }
    cat(paste("\rsampling...", i,"/", length(indices)))
  }
  cat("\n")

  if(turnOffC) {
    counts["CYS", "CYS"] <- 0 
  }
  
  return(counts)

}

#Get the PDB IDs for the given dataset
fetchPDBIDs <- function(dataset, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }

  sql = paste("select pdb_id FROM [", username, "@washington.edu].[",dataset,"_1.csv]",sep="")
  idlist <- fetchdata(sql)
  print(idlist)
  #minus one to skip the column headers
  ids <- rep("", length(idlist[[1]]) - 1)
  for(i in 2:length(idlist[[1]])) {
    ids[i - 1] <- as.character(idlist[[1]][i])
  }
  return(ids)
}

#Get the number of residues for each PDB in the given dataset
fetchLengths <- function(dataset, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }
  
  sql = paste("SELECT res_num FROM [", username, "@washington.edu].[",
    dataset,"_1.csv]",sep="")
  sqlList <- fetchdata(sql)
  lengths <- rep(NA, length(sqlList[[1]]) - 1)
  for(i in 2:length(sqlList[[1]])) {
    lengths[i - 1] <- as.integer(sqlList[[1]][[i]][1])
  }
  return(lengths)
}

#Get the charge and surface areas for each PDB in the given dataset
fetchChargeAndSA <- function(dataset, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }

  sql = paste("select charge,surface_area FROM [", username, "@washington.edu].[",dataset,"_1.csv] WHERE charge IS NOT NULL",sep="")
  sqlList <- fetchdata(sql)
  #minus one to skip the column headers
  charges <- rep(NA, length(sqlList[[1]]) - 1)
  SAs <- rep(NA, length(sqlList[[1]]) - 1)
  for(i in 2:length(sqlList[[1]])) {
    charges[i - 1] <- as.double(sqlList[[1]][[i]][1])
    SAs[i - 1] <- as.double(sqlList[[1]][[i]][2])
    
  }
  return(data.frame(charge=charges, surface_area=SAs))
  
}

#Get the list of surface residues in the given pdbid, where surface means their surface area ratio is greater than the given cutoff
fetchSurfResidues <- function(dataset, cutoff, pdbid, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }
  
  sql = paste("SELECT res_type FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE pdb_id =\'", pdbid, "\' AND res_surface_area_ratio > ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")

  rlist <- fetchdata(sql)
  rs <- rep("", length(rlist[[1]]) - 1)
  for(i in 2:length(rlist[[1]])) {
    rs[i - 1] <- as.character(rlist[[1]][i])
  }
  return(rs)
}

#Get all the surface residues in the given dataset. Normalize turns it into a list of distributions (one per PDB).
fetchAllSurfResidues <- function(dataset, cutoff, normalize=FALSE, username=FALSE, gapTol=0, secStruct=NULL) {
  if (username == FALSE) {
    username = myUsername
  }
  
  if(is.null(secStruct)) {
   sql = paste("SELECT pdb_id, res_type FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE res_surface_area_ratio > ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")
  } else{
    sql = paste("SELECT pdb_id, res_type FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE res_surface_area_ratio > ", cutoff,
      " AND res_surface_area_ratio IS NOT NULL AND structure = \'", secStruct ,"\'", sep="")
  }
  rlist <- fetchdata(sql)

  sql = paste("SELECT DISTINCT pdb_id FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE res_surface_area_ratio > ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")
  
  plist <- fetchdata(sql)
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
      data[pastid,names(rs.table)[j]] <- rs.table[j]
    }
  }
  
  return(data)
}

#Gets all the residues below the given cutoff, same as fetchAllSurfResidues
fetchAllBuriedResidues <- function(dataset, cutoff, normalize=FALSE, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }
  
  sql = paste("SELECT pdb_id, res_type FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE res_surface_area_ratio < ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")

  rlist <- fetchdata(sql)

  sql = paste("SELECT DISTINCT pdb_id FROM [", username, "@washington.edu].[",dataset, "_2.csv] WHERE res_surface_area_ratio < ", cutoff,
    " AND res_surface_area_ratio IS NOT NULL", sep="")
  
  plist <- fetchdata(sql)
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
      if(i == length(rlist[[1]])) {
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
      data[pastid,names(rs.table)[j]] <- rs.table[j]
    }
  }
  
  return(data)
}

#Fetch pairs in the given PDBid.
fetchResiduePairs <- function(dataset, cutoff, pdbid, symm=TRUE, counter=NULL, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }

  if(is.null(counter)) {
    counter <- empty.paircounter(aalist, aalist)
  }
  
  sql = paste("SELECT res_index, res_type FROM [", username, "@washington.edu].[", dataset,
    "_2.csv] WHERE pdb_id=\'", pdbid,
    "\' AND res_surface_area_ratio > ", cutoff, " AND res_surface_area_ratio IS NOT NULL", sep="")
  rlist <- fetchdata(sql)
  rs <- rep("", length(rlist[[1]]) - 1)
  indices <- rep(NA, length(rlist[[1]]) - 1)
  for(i in 2:length(rlist[[1]])) {

    indices[i - 1] <- as.integer(rlist[[1]][[i]][1])
    rs[i - 1] <- as.character(rlist[[1]][[i]][2])
    
  }
  for(i in 1:(length(indices) - 1)) {
    if(indices[i + 1] == indices[i] + 1) {
      if(sum(rownames(counter) == rs[i]) >= 1) {
        counter[rs[i],rs[i+1]] <- counter[rs[i],rs[i+1]] + 1
        if(symm == TRUE && rs[i] != rs[i + 1]) {
          counter[rs[i + 1],rs[i]] <- counter[rs[i + 1],rs[i]] + 1
        }
      }
      else {
        cat(paste("Unable to identify pair", rs[i],rs[i+1],"\n"))
      }
    }
  }
  
  return(counter)
}

#Fetch pairs for a dataset, where the pairs are coming from the dataset_pairs table in SQLShare
fetchAllResiduePairs <- function(dataset, pdbIDs=fetchPDBIDs(dataset), symm=TRUE, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }

  sql = paste("SELECT * FROM [whitead@washington.edu].[",dataset,"_pairs]" ,sep="")
  rlist <- fetchdata(sql)
  pairlist <- vector(mode="list", length=length(pdbIDs))
  index <- 0
  curID <- ""
  
  for(i in 2:length(rlist[[1]])) {
    if(rlist[[1]][[i]][1] != curID) {
      curID = rlist[[1]][[i]][1]
      index <- index + 1
      print(index)
      pairlist[[index]] <- empty.paircounter(aalist, aalist)
    }
    
    r1 <- rlist[[1]][[i]][2]
    r2 <- rlist[[1]][[i]][3]
    if(sum(aalist == r1) + sum(aalist == r2) == 2) {
      pairlist[[index]][which(r1 == aalist), which(r2 == aalist)] <-
         pairlist[[index]][which(r1 == aalist), which(r2 == aalist)] + 1
      if(symm) {
        pairlist[[index]][which(r2 == aalist), which(r1 == aalist)] <-
          pairlist[[index]][which(r2 == aalist), which(r1 == aalist)] + 1
      }
    }
  }
  
  return(pairlist)
}

#Fetch all triplets, where they are coming from the triplets table
fetchAllResidueTriplets <- function(dataset, symm=TRUE, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }

  sql = paste("SELECT * FROM [whitead@washington.edu].[",dataset,"_triplets]" ,sep="")
  rlist <- fetchdata(sql)
  triples <- rep(NA, length(rlist[[1]]) - 1)
  for(i in 2:length(rlist[[1]])) {
    triples[i - 1] <- paste(as.character(rlist[[1]][[i]][1]), as.character(rlist[[1]][[i]][2]),
            as.character(rlist[[1]][[i]][3]), sep="")
  }
  
  return(triples)
}

#Get the fraction of surface residues for the given pdbid
getSurfFrac <- function(dataset, cutoff, pdbid, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }
  residues <- fetchSurfResidues(dataset,cutoff,pdbid)
  fracs <- rep(0,length(aalist))
  names(fracs) <- aalist
  rtable <- table(residues)
  for(i in 1:length(rtable)) {
    fracs[names(rtable)[i]] = rtable[names(rtable)[i]]
  }
  fracs <- fracs / sum(fracs)
  return(fracs)
}

#Get the gyration parameter for each protein
fetchGyration <- function(dataset, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }
  sql = paste("SELECT pdb_id, lambda_x, lambda_y, lambda_z FROM [", username, "@washington.edu].[",dataset, "_1.csv]", sep="")
  lambdalist <- fetchdata(sql)
  
  pdbids <- array("", length(lambdalist[[1]]) - 1)
  for(i in 2:length(lambdalist[[1]])) {
    pdbids[i - 1] <- lambdalist[[1]][[i]][1]
  }
  
  data <- empty.df(c("lambda_x","lambda_y","lambda_z"), pdbids)
  for(j in 2:length(lambdalist[[1]])) {
    for(k in 2:4) {
      data[j-1, k-1] <- lambdalist[[1]][[j]][k]
    }
  }
  return(data)
}

#Get the residue number for each protein
fetchResNum <- function(dataset, username=FALSE) {
  if (username == FALSE) {
    username = myUsername
  }
  sql = paste("SELECT pdb_id, res_num FROM [", username, "@washington.edu].[",dataset, "_1.csv]", sep="")
  resNumlist <- fetchdata(sql)

  pdbids <- array("", length(resNumlist[[1]]) - 1)
  for(i in 2:length(resNumlist[[1]])) {
    pdbids[i - 1] <- resNumlist[[1]][[i]][1]
  }

  data <- empty.df("res_Num",pdbids)
  for(j in 2:length(resNumlist[[1]])) {
    data[j-1, 1] <- resNumlist[[1]][[j]][2]
  }

  return(data)
}
  

#Add error bars to a bargraph
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#Create an empty dataframe
empty.df<- function(cnames, rnames, default=NA){
        df<-data.frame(matrix(rep(default,length(cnames)*length(rnames)), nrow=length(rnames)))
        colnames(df)<-cnames
        rownames(df) <- rnames
        return(df)
}

#Create an empty matrix for use in counting pairs
empty.paircounter <- function(cnames, rnames){

  counter <- matrix(rep(0,length(cnames) * length(rnames)), nrow=length(rnames) )
  rownames(counter) <- rnames
  colnames(counter) <- cnames
           
  return(counter)
}

#load the entire pairlist for a given dataset without using the pairs SQLShare table
loadpairs <- function(dataset, cutoff, ids) {
  pairlist <- vector(mode="list", length=length(ids))

  cat(paste("loading", dataset, "pairs\n"))
  
  for(i in 1:length(ids)) {
    cat(ids[i])
    pairlist[[i]] <- fetchResiduePairs(dataset, cutoff, ids[i], symm=TRUE)
#    cat(paste("\r",i,"/",length(ids), " "))
  }
  cat("\n")

  return(pairlist)
}
