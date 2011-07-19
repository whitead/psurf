#!/usr/bin/env Rscript

source("/home/whitead/Documents/ProteinSurfaces/circlecorr.R")
library(RCurl)

username <- "whitead"
apikey <- "58f90137316aedb538b85a54955173c0"

cutoff <- 0.3
aalist = read.table("/home/whitead/Documents/ProteinSurfaces/aalist", header=T)$res
aalist.sh = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
labels.sh.ord <- c("K", "R", "H", "E", "D", "N", "Q", "G", "S", "T", "A", "I", "L", "M", "P", "V", "C", "F", "W", "Y")
fetchdata <- function(sql) {
  host = "https://sqlshare-rest.cloudapp.net"
  selector = "/REST.svc/v1/db/file"
  query=paste("?sql=",URLencode(sql),sep="")

  data = getURL(paste(host,selector,query,sep=""),
httpheader=c(Authorization =paste("ss_apikey ", username, "@washington.edu :", apikey, sep="")), verbose = FALSE, ssl.verifypeer =
FALSE)

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

fetchLengths <- function(dataset) {
  sql = paste("SELECT res_num FROM [", username, "@washington.edu].[",
    dataset,"_1.csv]",sep="")
  sqlList <- fetchdata(sql)
  lengths <- rep(NA, length(sqlList[[1]]) - 1)
  for(i in 2:length(sqlList[[1]])) {
    lengths[i - 1] <- as.integer(sqlList[[1]][[i]][1])
  }
  return(lengths)
}

fetchChargeAndSA <- function(dataset) {

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
      data[curid,names(rs.table)[j]] <- rs.table[j]
    }
  }
  
  return(data)
}

fetchAllBuriedResidues <- function(dataset, cutoff, normalize=FALSE) {
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

fetchResiduePairs <- function(dataset, cutoff, pdbid, symm=TRUE, counter=NULL) {

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

fetchAllResiduePairs <- function(dataset, pdbIDs=fetchPDBIDs(dataset), symm=TRUE) {


  sql = paste("SELECT * FROM [whitead@washington.edu].[",dataset,"_pairs]" ,sep="")
  print(sql)
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

fetchAllResidueTriplets <- function(dataset, symm=TRUE) {

  sql = paste("SELECT * FROM [whitead@washington.edu].[",dataset,"_triplets]" ,sep="")
  rlist <- fetchdata(sql)
  triples <- rep(NA, length(rlist[[1]]) - 1)
  for(i in 2:length(rlist[[1]])) {
    triples[i - 1] <- paste(as.character(rlist[[1]][[i]][1]), as.character(rlist[[1]][[i]][2]),
            as.character(rlist[[1]][[i]][3]), sep="")
  }
  
  return(triples)
}

getSurfFrac <- function(dataset, cutoff, pdbid) {
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

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

empty.df<- function(cnames, rnames, default=NA){
        df<-data.frame(matrix(rep(default,length(cnames)*length(rnames)), nrow=length(rnames)))
        colnames(df)<-cnames
        rownames(df) <- rnames
        return(df)
}

empty.paircounter <- function(cnames, rnames){

  counter <- matrix(rep(0,length(cnames) * length(rnames)), nrow=length(rnames) )
  rownames(counter) <- rnames
  colnames(counter) <- cnames
           
  return(counter)
}

correlationPicture <- function(matrix, name) {

  cm <- cor(matrix)
  cm[cm == 1] <- 0
  ms <- apply(matrix, MARGIN=2, FUN=mean)
  cairo_pdf(paste(name, "_corr.pdf", sep=""), width=7, height=7, pointsize=8)
  circle.corr(cm, ms, order=F, bg="gray50", col=colorRampPalette(c("blue", "white", "red"))(50))
  graphics.off()
  cairo_pdf(paste(name, "_corr_PCA.pdf", sep=""), width=7, height=7, pointsize=8)
  circle.corr(cm, ms, order=T, bg="gray50", col=colorRampPalette(c("blue", "white", "red"))(50))
  graphics.off()
}

plotpairs <- function(pairslist, fracs, svars, name="pairs", sigNumber=10, bootstrapNumber=20) {

  psum <- pairslist[[1]]
  #find most significant
  for(i in 2:length(pairslist)) {
    psum <- psum + pairslist[[i]]
  }

  counts <- sum(psum)

  
  pmax <- c()
  pmax.which <- c()
  pmax.indices <- matrix(rep(NA, 2 * sigNumber), ncol=2)
  pmax.background <- c()
  pmax.background.var <- c()
  
  for(i in 1:sigNumber) {
    pmax <- c(pmax, max(psum))
    pmax.which <- c(pmax.which, which.max(psum))
    pmax.indices[i,1] <- ceiling(pmax.which[i] / nrow(psum))
    pmax.indices[i,2] <- (pmax.which[i] - 1) %% nrow(psum) + 1

    anames <- sort(c(as.character(aalist.sh[pmax.indices[i,1]]), as.character(aalist.sh[pmax.indices[i,2]])))
    if(i == 1) {
      names(pmax) <- paste(anames[1], anames[2], sep="-")
    } else {
      names(pmax) <- c(names(pmax)[-i], paste(anames[1], anames[2], sep="-"))
    }

    if(pmax.indices[i,1] == pmax.indices[i,2]) {
      pmax.background <- c(pmax.background, fracs[pmax.indices[i,1]] * counts * fracs[pmax.indices[i,2]] / 2)
      pmax.background.var <- c(pmax.background.var, counts ** 2 * (svars[pmax.indices[i,1]] * fracs[pmax.indices[i,1]] ** 2 ))
    } else {

      pmax.background <- c(pmax.background, fracs[pmax.indices[i,1]] * counts *
                           fracs[pmax.indices[i,2]])
      pmax.background.var <- c(pmax.background.var, counts ** 2 * (svars[pmax.indices[i,1]] * fracs[pmax.indices[i,2]] ** 2 + svars[pmax.indices[i,2]] * fracs[pmax.indices[i,1]] ** 2))
      
    }
    #Zero it out
    psum[pmax.which] <- 0
    if(pmax.indices[i,1] != pmax.indices[i,2]) {
      psum[pmax.indices[i,1], pmax.indices[i,2]] <- 0
      psum[pmax.indices[i,2], pmax.indices[i,1]] <- 0
    }
    
    
  }

  pmax.boot <- matrix(rep(0,bootstrapNumber * sigNumber), nrow=bootstrapNumber)
  #bootstrap
  for(i in 1:bootstrapNumber) {
    isample <- sample(length(pairslist), replace=TRUE)
    psample <- lapply(pairslist[isample], FUN=function(x) x[pmax.which])
    for(j in 1:length(pairslist)) {
      for(k in 1:sigNumber) {
        pmax.boot[i,k] <- pmax.boot[i,k] +  psample[[j]][k]
      }
    }
    cat(paste("\r bootsrap", i, "/", bootstrapNumber,"\n"))
  }
  
  #use quantiles for error bar plotting
  eet <- matrix(c(apply(pmax.boot, MARGIN=2, FUN=function(x) {quantile(x,c(0.975))[1] - median(x)}), sqrt(pmax.background.var)), byrow=T, nrow=2)
  eeb <- matrix(c(apply(pmax.boot, MARGIN=2, FUN=function(x) {median(x) - quantile(x,c(0.025))[1]}), sqrt(pmax.background.var)), byrow=T, nrow=2)
  
  yy <- matrix(c(pmax, pmax.background), nrow=2, byrow=T)
  print(yy)
  cairo_pdf(paste(name, ".pdf", sep=""), width=4.5, height=3.42, pointsize=6)
  par(family="LMRoman10", bg="transparent")
  barx <- barplot(yy, axis.lty=1, col=c("gray50", "gray90"), main="", xlab="Amino Acid Pair", ylab="Observations", beside=T, ylim=c(0,max(eet) + max(yy)), names.arg=names(pmax))
  error.bar(barx,yy,lower=eeb,upper=eet, length=0.04)
  legend("topright", col=c("gray50", "gray90"), legend=c("Observed", "Background"), ncol=1, pch=15)
  graphics.off()
}

loadpairs <- function(dataset, cutoff, ids) {
  pairlist <- vector(mode="list", length=length(ids))

  cat(paste("loading", dataset, "pairs\n"))
  
  for(i in 1:length(ids)) {
    pairlist[[i]] <- fetchResiduePairs(dataset, cutoff, ids[i], symm=TRUE)
    cat(paste("\r",i,"/",length(ids), " "))
  }
  cat("\n")

  return(pairlist)
}

#Get ratios

efracs.surf <- fetchAllSurfResidues("eh1", cutoff, normalize=FALSE)
pfracs.surf <- fetchAllSurfResidues("cph1", cutoff, normalize=FALSE)
rfracs.surf <- fetchAllSurfResidues("h1_w", cutoff, normalize=FALSE)

efracs.interior <- fetchAllBuriedResidues("eh1", 0.05, normalize=FALSE)
pfracs.interior <- fetchAllBuriedResidues("cph1", 0.05, normalize=FALSE)
rfracs.interior <- fetchAllBuriedResidues("h1_w", 0.05, normalize=FALSE)

save.image()

rems.surf <- apply(efracs.surf, MARGIN=2, FUN=sum)
rems.interior <- apply(efracs.interior,MARGIN=2, FUN=sum)
rems <- rems.surf / rems.interior
revs <- apply(efracs.surf, MARGIN=2, FUN=var) * 1 / rems.interior + apply(efracs.interior, MARGIN=2, FUN=var) * rems.surf ** 2 * 1 / rems.interior ** 4

rpms.surf <- apply(pfracs.surf, MARGIN=2, FUN=sum)
rpms.interior <- apply(pfracs.interior,MARGIN=2, FUN=sum)
rpms <- rpms.surf / rpms.interior
rpvs <- apply(pfracs.surf, MARGIN=2, FUN=var) * 1 / rpms.interior + apply(pfracs.interior, MARGIN=2, FUN=var) * rpms.surf ** 2 * 1 / rpms.interior ** 4

rrms.surf <- apply(rfracs.surf, MARGIN=2, FUN=sum)
rrms.interior <- apply(rfracs.interior,MARGIN=2, FUN=sum)
rrms <- rrms.surf / rrms.interior
rrvs <- apply(rfracs.surf, MARGIN=2, FUN=var) * 1 / rrms.interior + apply(rfracs.interior, MARGIN=2, FUN=var) * rrms.surf ** 2 * 1 / rrms.interior ** 4





cairo_pdf("ratios.pdf", width=4.5, height=3.42, pointsize=8)
par(family="LMRoman10", bg="transparent", cex=0.8)
yy <- matrix(c(rrms, rpms, rems), byrow=T, nrow=3)
ee <- matrix(c(sqrt(rrvs), sqrt(rpvs), sqrt(revs)), byrow=T, nrow=3)
barx <- barplot(yy, col=c("light gray", "light blue", "dark red"), main="", xlab="Amino Acid", ylab="Density", beside=T, ylim=c(0,max(ee) + max(yy)), names.arg=aalist.sh)
error.bar(barx,yy,ee,length=0.02)
legend("topright", col=c("light gray", "light blue", "dark red"), legend=c("Human", "Cytoplasma", "Extracellular"), ncol=1, pch=15)
graphics.off()


#get fractions

efracs <- fetchAllSurfResidues("eh1", cutoff, normalize=TRUE)
pfracs <- fetchAllSurfResidues("cph1", cutoff, normalize=TRUE)
rfracs <- fetchAllSurfResidues("h1_w", cutoff, normalize=TRUE)

#get pbd ids
hids <- fetchPDBIDs("h1_w")
pids <- fetchPDBIDs("cph1")
eids <- fetchPDBIDs("eh1")


#Make charge histograms
cols <- colorRampPalette(c("red","white","blue"))
datasets <- c("cph1", "eh1", "h1_w")
for(i in 1:length(datasets)) {
  d <- datasets[i]
  cat(paste("loading", d, "charge densities\n"))
  csa <- fetchChargeAndSA(d)
  save.image()
  cdens <- csa$charge[!is.na(csa$charge)] / csa$surface_area[!is.na(csa$charge)]
  cat(paste("Median density for", d, "is", median(cdens), "\n"))
  cat(paste("[Pos] [Neut] [Neg]", as.double(sum(cdens > 0)) / length(cdens), as.double(sum(cdens == 0)) / length(cdens), as.double(sum(cdens < 0)) / length(cdens)))
  cairo_pdf(paste(d,"_charges.pdf",sep=""), width=3.42, height=3.42, pointsize=8)
  par(family="LMRoman10")
  absmax <- max(cdens, -cdens)
  brs <- seq(from = -absmax, to = absmax, length.out=40)
  hist(cdens, breaks=brs, col=cols(40), xlab="Charge Desnsity", family="LMRoman10", main="") 
  graphics.off()
}




#Make PC plot
colors <- rep("gray50", length(hids))
for(i in 1:length(hids)) {
  if(sum(pids == hids[i]) >= 1) {
    colors[i] <- "light blue"
    if(sum(eids == hids[i]) >= 1) {
      colors[i] <- "green"
    }
  } else if(sum(eids == hids[i]) >= 1) {
    colors[i] = "dark red"
  }
  
}
  
save.image()

cairo_pdf("pca.pdf", width=3.42, height=3.42, pointsize=6)
par(family="LMRoman10", fg="light gray")

pc <- princomp(rfracs)

#plot(pc$scores[,1], pc$scores[,2], col=colors, xlab="PC 1", ylab="PC 2")
gsfrac <- apply(rfracs[,c("GLY", "SER")], MARGIN=1, sum)
kefrac <- apply(rfracs[,c("GLU", "LYS")], MARGIN=1, sum)
plot(gsfrac[colors == "gray50"], kefrac[colors=="gray50"], col="gray50", xlab="GLY + SER", ylab="GLU + LYS")
points(gsfrac[colors == "light blue"], kefrac[colors=="light blue"], col="light blue", xlab="GLY + SER", ylab="GLU + LYS")
points(gsfrac[colors == "dark red"], kefrac[colors=="dark red"], col="dark red", xlab="GLY + SER", ylab="GLU + LYS")
points(gsfrac[colors == "green"], kefrac[colors=="green"], col="green", xlab="GLY + SER", ylab="GLU + LYS")
legend("topright", col=c("gray50", "light blue", "dark red", "green"),
       legend=c("Uncategorized", "Cytoplasma", "Extracellular", "Both"), ncol=1,
       pch=1, text.col="black")

graphics.off()


#surface fractions with variance

pms <- apply(pfracs, MARGIN=2, mean)
pvs <- apply(pfracs, MARGIN=2,var)

ems <- apply(efracs, MARGIN=2, mean)
evs <- apply(efracs, MARGIN=2, var)

hms <- apply(rfracs, MARGIN=2, mean)
hvs <- apply(rfracs, MARGIN=2, var)

cairo_pdf("comb.pdf", width=4.5, height=3.42, pointsize=8)
par(family="LMRoman10", bg="transparent", cex=0.8)
yy <- matrix(c(hms, pms, ems), byrow=T, nrow=3)
ee <- matrix(c(sqrt(hvs / nrow(rfracs)),sqrt(pvs / nrow(pfracs)), sqrt(evs / nrow(efracs))), byrow=T, nrow=3)
barx <- barplot(yy, col=c("light gray", "light blue", "dark red"), main="", xlab="Amino Acid", ylab="Density", beside=T, ylim=c(0,ee[which.max(yy)] + max(yy)), names.arg=aalist.sh)
error.bar(barx,yy,ee,length=0.02)
legend("topright", col=c("light gray", "light blue", "dark red"), legend=c("Human", "Cytoplasma", "Extracellular"), ncol=1, pch=15)
graphics.off()

cairo_pdf("comb_black.pdf", width=7, height=3.42, pointsize=8)
par(family="LMRoman10", bg="black", cex=0.8, fg="white", col.axis="white", col.lab="white")
yy <- matrix(c(hms, pms, ems), byrow=T, nrow=3)
ee <- matrix(c(sqrt(hvs / nrow(rfracs)),sqrt(pvs / nrow(pfracs)), sqrt(evs / nrow(efracs))), byrow=T, nrow=3)
barx <- barplot(yy, col=c("blue", "dark gray", "dark red"), main="", xlab="Amino Acid", ylab="Density", beside=T, ylim=c(0,ee[which.max(yy)] + max(yy)), names.arg=aalist.sh)
error.bar(barx,yy,ee,length=0.02)
legend("topright", col=c("blue", "dark gray", "dark red"), legend=c("Human", "Cytoplasma", "Extracellular"), ncol=1, pch=15)
graphics.off()
               
#surface correlation matrices
correlationPicture(rfracs, "h1_w")
correlationPicture(pfracs, "cph1")
correlationPicture(efracs, "eh1")

#load pairs


hpairs <- loadpairs("h1_w", cutoff, hids)
ppairs <- loadpairs("cph1", cutoff, pids)
epairs <- loadpairs("eh1", cutoff, eids)

save.image()

#bootstrap and plot em
plotpairs(epairs, ems, evs / nrow(efracs), name="eh_pairs")
plotpairs(ppairs, pms, pvs / nrow(pfracs), name="cph_pairs")
plotpairs(hpairs, hms, hvs / nrow(rfracs), name="h_pairs")

