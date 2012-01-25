#!/usr/bin/env Rscript

source("SQLShareLib.R")

#Make a picture of a correlation matrix
correlationPicture <- function(matrix, name) {

  cm <- cor(matrix)
  cat(paste("Correlation matrix for:", name))
  print(cm)
  cm[cm == 1] <- 0
  ms <- apply(matrix, MARGIN=2, FUN=mean)
  cairo_pdf(paste(name, "_corr.pdf", sep=""), width=7, height=7, pointsize=12)
  circle.corr(cm, ms, order=F, bg="gray50", col=colorRampPalette(c("blue", "white", "red"))(50))
  graphics.off()
  cairo_pdf(paste(name, "_corr_PCA.pdf", sep=""), width=7, height=7, pointsize=12)
  circle.corr(cm, ms, order=T, bg="gray50", col=colorRampPalette(c("blue", "white", "red"))(50))
  graphics.off()
  
}


plotpairs <- function(pairslist, fracs, svars, name="pairs", sigNumber=10, bootstrapNumber=5000) {

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
  png(paste(name, ".png", sep=""), width=3.3 * 500, height=2.5 * 500, pointsize=10, res=500)
  par(family="LMSans10", mar=c(3,4,0.5,0.5), cex.axis=0.65)
  barx <- barplot(yy, axis.lty=1, col=c("gray50", "gray90"), main="", xlab="Amino Acid Pair", ylab="Frequency", beside=T, ylim=c(0,max(eet) + max(yy)), names.arg=names(pmax))
  error.bar(barx,yy,lower=eeb,upper=eet, length=0.03)
  legend("topright", col=c("gray50", "gray90"), legend=c("Observed", "Background"), ncol=1, pch=15)
  graphics.off()

  cairo_pdf(paste(name, ".pdf", sep=""), width=3.42, height=2.58, pointsize=8)
  par(family="LMSans10", mar=c(3,4,0.5,0.5), cex.axis=0.65)
  barx <- barplot(yy, axis.lty=1, col=c("gray50", "gray90"), main="", xlab="Amino Acid Pair", ylab="Frequency", beside=T, ylim=c(0,max(eet) + max(yy)), names.arg=names(pmax))
  error.bar(barx,yy,lower=eeb,upper=eet, length=0.03)
  legend("topright", col=c("gray50", "gray90"), legend=c("Observed", "Background"), ncol=1, pch=15)
  graphics.off()
}


#Get ratios

efracs.surf <- fetchAllSurfResidues("eh2_w_nogaps", cutoff, normalize=FALSE)
pfracs.surf <- fetchAllSurfResidues("cph2_w_nogaps", cutoff, normalize=FALSE)
rfracs.surf <- fetchAllSurfResidues("h2_w_nogaps", cutoff, normalize=FALSE)

efracs.interior <- fetchAllBuriedResidues("eh2_w_nogaps", 0.05, normalize=FALSE)
pfracs.interior <- fetchAllBuriedResidues("cph2_w_nogaps", 0.05, normalize=FALSE)
rfracs.interior <- fetchAllBuriedResidues("h2_w_nogaps", 0.05, normalize=FALSE)

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





cairo_pdf("ratios.pdf", width=4.5, height=3.3, pointsize=12)
par(family="LMRoman10", bg="transparent", cex=0.8)
yy <- matrix(c(rrms, rpms, rems), byrow=T, nrow=3)
ee <- matrix(c(sqrt(rrvs), sqrt(rpvs), sqrt(revs)), byrow=T, nrow=3)
barx <- barplot(yy, col=c("light gray", "light blue", "dark red"), main="", xlab="Amino Acid", ylab="Fraction", beside=T, ylim=c(0,max(ee) + max(yy)), names.arg=aalist.sh)
error.bar(barx,yy,ee,length=0.02)
legend("topright", col=c("light gray", "light blue", "dark red"), legend=c("Human", "Cytoplasma", "Extracellular"), ncol=1, pch=15)
graphics.off()


#get fractions

efracs <- fetchAllSurfResidues("eh2_w_nogaps", cutoff, normalize=TRUE)
pfracs <- fetchAllSurfResidues("cph2_w_nogaps", cutoff, normalize=TRUE)
rfracs <- fetchAllSurfResidues("h2_w_nogaps", cutoff, normalize=TRUE)

#get pbd ids
hids <- fetchPDBIDs("h2_w_nogaps")
pids <- fetchPDBIDs("cph2_w_nogaps")
eids <- fetchPDBIDs("eh2_w_nogaps")


#Make charge histograms
cols <- colorRampPalette(c("red","white","blue"))
datasets <- c("cph2_w_nogaps", "eh2_w_nogaps", "h2_w_nogaps")
breakNum <- c(20,20,40)
for(i in 1:length(datasets)) {
  d <- datasets[i]
  cat(paste("loading", d, "charge densities\n"))
  csa <- fetchChargeAndSA(d)
  save.image()
  cdens <- csa$charge[!is.na(csa$charge)] / csa$surface_area[!is.na(csa$charge)]
  cat(paste("Median density for", d, "is", median(cdens), "\n"))
  cat(paste("Absolute density median for", d, "is", median(abs(cdens)), "\n"))
  cat(paste("[Pos] [Neut] [Neg]", as.double(sum(cdens > 0)) / length(cdens), as.double(sum(cdens == 0)) / length(cdens), as.double(sum(cdens < 0)) / length(cdens)))
  png(paste(d,"_charges.png",sep=""), width=3.3*500, height=3.3*500, pointsize=10, res=500)
  par(family="LMSans10", cex=0.8)
  absmax <- max(cdens, -cdens)
  brs <- seq(from = -absmax, to = absmax, length.out=breakNum[i])
  hist(cdens, breaks=brs, col=cols(breakNum[i]), xlab="Charge Desnsity", family="LMSans10", main="") 
  graphics.off()

  cairo_pdf(paste(d,"_charges.pdf",sep=""), width=3.42, height=2.58, pointsize=10)
  par(family="LMSans10", cex=0.8)
  absmax <- max(cdens, -cdens)
  brs <- seq(from = -absmax, to = absmax, length.out=breakNum[i])
  hist(cdens, breaks=brs, col=cols(breakNum[i]), xlab="Charge Desnsity", family="LMSans10", main="") 
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

cairo_pdf("pca.pdf", width=10, height=10, pointsize=12)
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

png("comb.png", width=3.3 * 500, height=2.5 * 500, pointsize=8, res=500)
par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
yy <- matrix(c(hms, pms, ems), byrow=T, nrow=3)
ee <- matrix(c(sqrt(hvs / nrow(rfracs)),sqrt(pvs / nrow(pfracs)), sqrt(evs / nrow(efracs))), byrow=T, nrow=3)
barx <- barplot(yy, col=c("light gray", "light blue", "dark red"), main="", xlab="Amino Acid", ylab="Fraction", beside=T, ylim=c(0,ee[which.max(yy)] + max(yy)), names.arg=aalist.sh)
error.bar(barx,yy,ee,length=0.01)
legend("topright", col=c("light gray", "light blue", "dark red"), legend=c("Human", "Cytoplasma", "Extracellular"), ncol=1, pch=15)
graphics.off()

cairo_pdf("comb.pdf", width=3.42, height=2.58, pointsize=8)
par(mar=c(3,4,0.5,0.5), cex.axis=0.65, family="LMSans10")
yy <- matrix(c(hms, pms, ems), byrow=T, nrow=3)
ee <- matrix(c(sqrt(hvs / nrow(rfracs)),sqrt(pvs / nrow(pfracs)), sqrt(evs / nrow(efracs))), byrow=T, nrow=3)
barx <- barplot(yy, col=c("blue", "gray45", "dark red"), main="", xlab="Amino Acid", ylab="Fraction", beside=T, ylim=c(0,ee[which.max(yy)] + max(yy)), names.arg=aalist.sh)
error.bar(barx,yy,ee,length=0.02)
legend("topright", col=c("blue", "gray45", "dark red"), legend=c("Human", "Cytoplasma", "Extracellular"), ncol=1, pch=15)
graphics.off()

cairo_pdf("comb_black.pdf", width=7, height=3.3, pointsize=12)
par(family="LMRoman10", bg="black", cex=0.8, fg="white", col.axis="white", col.lab="white")
yy <- matrix(c(hms, pms, ems), byrow=T, nrow=3)
ee <- matrix(c(sqrt(hvs / nrow(rfracs)),sqrt(pvs / nrow(pfracs)), sqrt(evs / nrow(efracs))), byrow=T, nrow=3)
barx <- barplot(yy, col=c("blue", "dark gray", "dark red"), main="", xlab="Amino Acid", ylab="Fraction", beside=T, ylim=c(0,ee[which.max(yy)] + max(yy)), names.arg=aalist.sh)
error.bar(barx,yy,ee,length=0.02)
legend("topright", col=c("blue", "dark gray", "dark red"), legend=c("Human", "Cytoplasma", "Extracellular"), ncol=1, pch=15)
graphics.off()
               
#surface correlation matrices
correlationPicture(rfracs, "h2_w_nogaps")
correlationPicture(pfracs, "cph2_w_nogaps")
correlationPicture(efracs, "eh2_w_nogaps")

#load pairs


hpairs <- loadpairs("h2_w_nogaps", cutoff, hids)
ppairs <- loadpairs("cph2_w_nogaps", cutoff, pids)
epairs <- loadpairs("eh2_w_nogaps", cutoff, eids)

save.image()

#bootstrap and plot em
plotpairs(epairs, ems, evs, name="eh_pairs")
plotpairs(ppairs, pms, pvs, name="cph_pairs")
plotpairs(hpairs, hms, hvs, name="h_pairs")



