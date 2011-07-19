#! /usr/bin/env Rscript

job <- commandArgs()[6]
data <- read.table(job, header=T)
anames <- colnames(data[-c(1,2)])
rnum <- length(anames)

labels.ord <- c("K", "R", "H", "E", "D", "N", "Q", "G", "S", "T", "A", "I", "L", "M", "P", "V", "C", "F", "W", "Y")
labels.one <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
labels.ord.f <- c("K", "R", "H", "E", "D", "N", "Q", "G", "S", "T", "A", "I", "L", "M", "P", "V", "C", "F", "W", "Y", "*")
labels.one.f <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")

matPlot <- function(x) {

  rownames(x) <- labels.one
  colnames(x) <- labels.one.f

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
interaction <- data.frame(matrix(rep(0,bootstrap * length(labels.one)), nrow=bootstrap))
colnames(interaction) <- labels.one

water <- data.frame(matrix(rep(0,bootstrap * length(labels.one)), nrow=bootstrap))
colnames(water) <- labels.one

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

cairo_pdf("contact_plots.pdf", width=3.42,height=4.5, pointsize=8)
par(family="LMRoman10", mfrow=c(2,1), cex.axis=0.6, mar=c(2,4,1,2), cex.main=0.8)
xp <- barplot(wyy[2,], ylim=c(0,max(wyy[3,] + 0.03)), ylab="Water per Contact", col="light blue")
error.bar(xp, wyy[2,], lower=(wyy[2,] - wyy[1,]), upper=(wyy[3,] - wyy[2,]), length=0.03)
par(mar=c(5,4,2,2))
xp <- barplot(iyy[2,], ylim=c(0,max(iyy[3,] + 0.03)), ylab="Proportion Interacting", col="light blue")
error.bar(xp, iyy[2,], lower=(iyy[2,] - iyy[1,]), upper=(iyy[3,] - iyy[2,]), length=0.03)
graphics.off()



#BELOW IS OLD CODE
q()


bpairs.mean <- calculateBackground(counts)
bpairs.var <- empty.df(labels.one.f, labels.one, default=0)

ppairs.mean <- calculatePairs(counts)
ppairs.var <- empty.df(labels.one.f, labels.one, default=0)

logodd.mean <- -log(ppairs.mean / bpairs.mean)
logodd.var <- empty.df(labels.one.f, labels.one, default=0)

delta <- empty.df(labels.one.f, labels.one, default=0)


for(i in 1:bootstrap) {
  counts <- sampleCounts()
  bpairs <- calculateBackground(counts)
  ppairs <- calculatePairs(counts)
  logodd <- -log(as.matrix(ppairs) / as.matrix(bpairs))
  
  delta <- logodd - logodd.mean
  logodd.mean <- logodd.mean + delta  / i
  logodd.var <- logodd.var + delta * (logodd - logodd.mean)
  delta <- bpairs - bpairs.mean
  bpairs.mean <- bpairs.mean + delta  / i
  bpairs.var <- bpairs.var + delta * (bpairs - bpairs.mean)
  delta <- ppairs - ppairs.mean
  ppairs.mean <- ppairs.mean + delta  / i
  ppairs.var <- ppairs.var + delta * (ppairs - ppairs.mean)
  cat(paste("\r bootsrap", i, "/", bootstrap,"\n"))
}

save.image()

write.table(logodd.mean, file=paste(job, "_table_logdd.txt", sep=""))
write.table(sqrt(logodd.var), file=paste(job, "_table_logodd_err.txt", sep=""))
write.table(ppairs.mean, file=paste(job, "_table_ppairs.txt", sep=""))
write.table(sqrt(ppairs.var), file=paste(job, "_table_ppairs_err.txt", sep=""))
write.table(bpairs.mean, file=paste(job, "_table_bpairs.txt", sep=""))
write.table(sqrt(bpairs.var), file=paste(job, "_table_bpairs_err.txt", sep=""))

print(logodd.mean)
print(logodd.var)



cols <- colorRampPalette(c("black", "white", "red"))

pdf(paste(job, "_true.pdf", sep=""), width=6, height=6)
par(cex.axis=0.5)
image(x=1:rnum, y=1:(rnum + 1), z=matPlot(ppairs.mean), axes=F, xlab="AAj", ylab="AAi", col=cols(50), main="Observed Contact Probabilities")
axis(1, at=1:rnum, labels=labels.ord)
axis(2, at=1:(rnum+1), labels=labels.ord.f)

image(x=1:rnum, y=1:(rnum + 1), z=matPlot(bpairs.mean), axes=F, xlab="AAj", ylab="AAi", col=cols(50), main="Expected Contact Probabilities")
axis(1, at=1:rnum, labels=labels.ord)
axis(2, at=1:(rnum + 1), labels=labels.ord.f)


logodd <- log(as.matrix(ppairs.mean) / as.matrix(bpairs.mean))
print(logodd)

image(x=1:rnum, y=1:(rnum + 1), z=matPlot(logodd), axes=F, xlab="AAj", ylab="AAi", col=cols(50), main="Contact Log Odds")
axis(1, at=1:rnum, labels=labels.ord)
axis(2, at=1:(rnum + 1), labels=labels.ord.f)

graphics.off()

sigNumber <- 30
pmax <- c()
pmax.var <- c()
pmax.which <- c()
pmax.indices <- matrix(rep(NA, 2 * sigNumber), ncol=2)
                       
for(i in 1:sigNumber) {
  pmax <- c(pmax, max(logodd))
  pmax.which <- c(pmax.which, which.max(logodd))
  pmax.indices[i,1] <- ceiling(pmax.which[i] / nrow(logodd))
  pmax.indices[i,2] <- (pmax.which[i]  - 1) %% nrow(logodd) + 1
  ri <- pmax.indices[i,1]
  ci <- pmax.indices[i,2]
  pmax.var <- c(pmax.var, bpairs.var[ri,ci] / bpairs.mean[ri,ci] ** 2 + ppairs.var[ri,ci] / ppairs.mean[ri,ci] ** 2)

  anames <- sort(c(labels.one[pmax.indices[i,1]], labels.one.f[pmax.indices[i,2]]))
  names(pmax) <-  c(names(pmax)[-i], paste(anames[1],anames[2], sep="-"))

  logodd[pmax.which] <- min(logodd)
  if(pmax.indices[i,1] != pmax.indices[i,2] && pmax.indices[i,2] != ncol(logodd)) {

    logodd[pmax.indices[i,1], pmax.indices[i,2]] <- min(logodd)
    logodd[pmax.indices[i,2], pmax.indices[i,1]] <- min(logodd)
    
  }  
}

yy <- pmax
ee <- sqrt(pmax.var)

cairo_pdf(paste(job, "_logodds.pdf"), width=18, height=5)
par(family="LMRoman10", bg="transparent")
barx <- barplot(yy, axis.lty=1, col="gray50", main="", xlab="Amino Acid Pair", ylab="Free Energy [-kT]", names.arg=names(pmax), ylim=c(0,max(ee) + max(yy)))
error.bar(barx, yy, ee, length=0.04)
graphics.off()
