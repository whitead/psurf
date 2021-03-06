source("SQLShareLib.R")
source("ProtLib.R")
source("circlecorr.R")

#Interaction Table
dataset <- "h2"
username <- "whitead"

chis <- getInteractionEnergy(dataset, username, countMatrix=fetchContacts("h1_surface_contacts.csv"))

interactionTable <- chis

surface <- fetchAllSurfResidues(paste(dataset,"w_nogaps", sep="_"), cutoff, normalize=FALSE, username)
interior <- fetchAllBuriedResidues(paste(dataset,"w_nogaps", sep="_"), cutoff, normalize=FALSE, username)
all <- fetchAllSurfResidues(paste(dataset,"w_nogaps", sep="_"), -1, normalize=FALSE, username)

distributions <- list(Proteins=all, Surface=surface, Interior=interior, Diff=(surface - all))



#sample all the contact energies
contacts <- fetchContacts(paste(paste(dataset,"backbone","contacts", sep="_"), "csv" ,sep="."), username)
bootstrap <- 2
chis.list <- vector("list", bootstrap)

for(b in 1:bootstrap) {

  chis <- getInteractionEnergy(countMatrix=contacts, random=TRUE)
  
  for(i in 1:length(distributions)) {

    dist <- distributions[[i]]
    energies <- matrix(0, nrow=nrow(dist), ncol=20)
    
    for(j in 1:nrow(dist)) {
      
      energies[j,] <- as.matrix(dist[j,]) %*% chis[1:20, 1:20]
    
    }
    chis <- rbind(chis, apply(energies, 2, median)) 
  }

  rownames(chis) <- c(rownames(chis)[-sapply(1:length(distributions), function(x) {nrow(chis) - x + 1})], "Proteins", "Surface", "Interior", "Surface - Proteins")
  
  chis.list[[b]] <- chis
}

#Construct a data.frame containing the medians and a data.frame containing the lower and upper

print(chis.list)
print(rownames(chis.list[[1]]))

interactionTable <- empty.df(colnames(chis.list[[1]]), rownames(chis.list[[1]]))
interactionTable.lower <- empty.df(colnames(chis.list[[1]]), rownames(chis.list[[1]]))
interactionTable.upper <- empty.df(colnames(chis.list[[1]]), rownames(chis.list[[1]]))

for(i in 1:nrow(chis.list[[1]])) {
  for(j in 1:ncol(chis.list[[1]])) {

    values <- as.double(lapply(chis.list, function(x) {x[i, j]}))
    qs <- quantile(values, c(0.025, 0.5, 0.975))
    interactionTable[i, j] <- qs[2]
    interactionTable.lower[i,j] <- qs[1]
    interactionTable.upper[i,j] <- qs[3]
  }
}

chis <- interactionTable[-sapply(1:length(distributions), function(x) {nrow(chis) - x + 1}),]
chis <- as.matrix(chis)

print(chis)
print(interactionTable)
print(interactionTable.lower)


#remove glycine
gindex <- which(colnames(interactionTable) == "GLY")

interactionTable <- interactionTable[-gindex, -gindex]
interactionTable.upper <- interactionTable.upper[-gindex, -gindex]
interactionTable.lower <- interactionTable.lower[-gindex, -gindex]

write.table(round(interactionTable, 2), file="Interaction_Table.txt", quote=FALSE, sep=" & ", eol="\\\\\n")


#Plot table
colGrad <- colorRampPalette(c("blue", "red"))
res <- 250
boxsize <- 0.25
png("aa_interaction_table.png", width=res * (4 + 1 + boxsize * nrow(chis)), height=res * (4 + 1 + boxsize * ncol(chis)), res=res)
par(family="LMRoman10", cex.sub=1.2, fg="dark gray")
chis["CYS", "CYS"] <- NA
image(1:nrow(chis), 1:ncol(chis), -chis, col=colGrad(500), xaxt="n", yaxt="n", xlab="Contacting", ylab="Amino Acid")
axis(1, 1:nrow(chis), c(aalist.sh, "Ar", "Po", "Ch", "Hy"))
print(ncol(chis))
axis(2, 1:ncol(chis), aalist.sh)
graphics.off()

colGrad <- colorRampPalette(c("white", "black"))
res <- 250
boxsize <- 0.25
png("aa_interaction_table_bw.png", width=res * (4 + 1 + boxsize * nrow(chis)), height=res * (4 + 1 + boxsize * ncol(chis)), res=res)
par(family="LMRoman10", cex.sub=1.2, fg="dark gray")
chis["CYS", "CYS"] <- NA
image(1:nrow(chis), 1:ncol(chis), -chis, col=colGrad(500), xaxt="n", yaxt="n", xlab="Contacting", ylab="Amino Acid")
axis(1, 1:nrow(chis), c(aalist.sh, "Ar", "Po", "Ch", "Hy"))
print(ncol(chis))
axis(2, 1:ncol(chis), aalist.sh)
graphics.off()

#Plot table with proportion interacting
#get proportion interacting first
surfContacts <- fetchContacts(paste(paste("h1","surface","contacts", sep="_"), "csv" ,sep="."), username)
scont <- sampleContacts(surfContacts, random=FALSE)
propInteract <- sapply(aalist, function(x) {1 - scont["FREE", x] / scont["TOTAL", x]})

#remove glycine
print(propInteract)
propInteract <- propInteract[-gindex]
print(propInteract)
#make the interaction table vary between -1 and 1, just cuz
interactMatrix <- interactionTable[1:19, 1:19]
interactMatrix <- as.matrix(-interactMatrix)
rownames(interactMatrix) <- aalist.sh[-gindex]
colnames(interactMatrix) <- aalist.sh[-gindex]
 
interval <- (max(interactMatrix) - min(interactMatrix))
disp <- max(interactMatrix) - interval / 2
interactMatrix <- 1.7 *(interactMatrix - disp) / interval
cairo_pdf("prop_interact_matrix.pdf", width=7, height=7, pointsize=12)
par(family="LMSans10")
circle.corr(interactMatrix, propInteract, order=F, yTextRot=0, bg="gray50",
            col=colorRampPalette(c("blue", "white", "red"))(200))
graphics.off()


#Make truncated version
colGrad <- colorRampPalette(c("blue", "red"))
trunc.residues <- c("GLU", "GLN", "ASP", "GLN", "LYS", "ARG", "SER", "ALA")
trunc.labs <- c("E", "Q", "D", "N", "K", "R", "S", "A")
trunc.cats <- c("Polar", "Charged", "Hydrophobic")
trunc.index <- sapply(trunc.residues, function(x) {which(rownames(chis) == x)})
trunc.index.cats <- sapply(trunc.cats, function(x) {which(rownames(chis) == x)})
trunc.labs.cats <- c("Po", "Ch", "Hy")


chis <- chis[c(trunc.index, trunc.index.cats), trunc.index]

res <- 300
boxsize <- 0.25
png("aa_trunc_interaction_table.png", width=res * (4 + 1 + boxsize * (length(trunc.cats) + length(trunc.labs))), height=res * (4 + 1 + boxsize * (length(trunc.labs))), res=res)
par(family="LMRoman10", cex.axis=1.0, fg="gray25", par=c(4, 4, 1, 1))
image(1:nrow(chis), 1:ncol(chis), chis, col=colGrad(500), xaxt="n", yaxt="n", xlab="Contacting", ylab="Amino Acid")
box()
axis(1, 1:(length(trunc.cats) + length(trunc.residues)), c(trunc.labs, trunc.labs.cats), xaxs="r")
print(ncol(chis))
axis(2, 1:length(trunc.residues), trunc.labs)
graphics.off()

#Now make bargraph
cairo_pdf("aa_protein_interactions.pdf", width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.6)
yy <- as.matrix(interactionTable[c("Surface", "Interior"),])
print(yy)
print(dim(yy))
print(dim(as.matrix(interactionTable.lower[c("Surface", "Interior"),])))
yy.lower <- yy - as.matrix(interactionTable.lower[c("Surface", "Interior"),])
yy.upper <- as.matrix(interactionTable.upper[c("Surface", "Interior"),]) - yy
yy.index <- order(yy[1,] - yy[2,])
aalist.sh.nog <- aalist.sh[-which(aalist.sh == "G")]
barx <- barplot(yy[,yy.index], col=c("gray25", "gray80"), main="", xlab="Amino Acid", ylab=expression(paste("Interaction Energy, E[", chi, "]", " kJ/mol")), beside=T, names.arg=aalist.sh.nog[yy.index], ylim=c(-30,30), space=c(0,0.4))
error.bar(barx, yy[,yy.index], lower=yy.lower[,yy.index], upper=yy.upper[,yy.index], length=0.02)
legend("top", col=c("gray25", "gray80"), cex=0.8, legend=c("Surface", "Buried"), pch=15)
graphics.off()

                       

sql <- paste("select * FROM [whitead@washington.edu].[GroEL_counts.csv]")
rawData <- fetchdata(sql)



groDist <- empty.df(rawData[[1]][[1]][-1], c(rawData[[1]][[2]][1], rawData[[1]][[3]][1]))
groDist[1:2,] <- matrix(unlist(lapply(rawData[[1]][-1], function(row) sapply(row[-1], as.integer))), nrow=2, byrow=T)

groDist[1, ] <- as.double(groDist["GroEL_Close", ]) / sum(as.double(groDist["GroL_Close", ]))
groDist[2, ] <- as.double(groDist["GroEL_Open", ]) / sum(as.double(groDist["GroE_Open", ]))


b <- array(0,20)

cairo_pdf('groDist_SE.pdf', width=3.42, height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.6)

barx <- barplot(as.matrix(groDist), col=c("gray15","gray75"), main="", xlab="Amino Acid", ylab="", beside=T, names.arg=aalist.sh, ylim = c(0.00,0.20))
legend("topright", col=c("gray15","gray75"), cex=0.8, legend=c("GroEL-GroES (Closed)","GroEL (Open)"), pch=15)
graphics.off()


#plot2 The plot of 5 species

hspRawCount <- fetchHSPResCount("hsp","in",username="whitead")

for (i in 1:nrow(hspRawCount)) {
  hspRawCount[i,] <- hspRawCount[i,] / sum(hspRawCount[i,])
}
hspRawCount <- as.matrix(hspRawCount)

hspDist <- matrix(0,nrow(hspRawCount)+1, ncol(hspRawCount))
rownames(hspDist) <- c("Thermo Thermopilius","E.Coli GroEL","Thermo GroEL","Group II","HSP90","Yeast CCT")
colnames(hspDist) <- colnames(hspRawCount)

hspDist[2,] <- hspRawCount["1SX4",]
hspDist[3,] <- hspRawCount["1WE3",]
hspDist[4,] <- hspRawCount["3KFB",]
hspDist[5,] <- hspRawCount["2CG9",]
hspDist[6,] <- hspRawCount["3P9D",]

#various statistics
print(sum(hspRawCount["1WE3", c("ASP", "GLU", "LYS", "ARG", "HIS")]))
print(sum(hspRawCount["1SX4", c("ASP", "GLU", "LYS", "ARG", "HIS")]))
print(sum(hspRawCount["3KFB", c("ASP", "GLU")]))
print(sum(hspRawCount["3KFB", c("LYS", "HIS", "ARG")]))

cutoff <- 0.3
psurf <- fetchAllSurfResidues("thermus_nogaps", cutoff, normalize=TRUE, "wenjunh")

chargeSum <- rep(0, nrow(psurf))

for (i in 1:nrow(psurf)) {
  chargeSum[i] <- sum(psurf[i,c("LYS","ARG","GLU","ASP","HIS")])
}
quantile(chargeSum, seq(0,1,0.01))

hspDist[1,] <- apply(psurf, MARGIN=2, FUN=median)

up <- matrix(0,6,ncol(psurf))
low <- matrix(0,6,ncol(psurf))

for (i in 1:ncol(psurf)) {
  up[1,i] <- quantile(psurf[,i],0.95) - median(psurf[,i])
  low[1,i] <- median(psurf[,i]) - quantile(psurf[,i],0.05)
}

cairo_pdf('HSP_Protein_Dist_Thermo.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
barx <- barplot(hspDist, beside=TRUE, col=c('blue','gray75','gray60','gray45','gray30','gray15'), ylim=c(0.0,0.3), names.arg=aalist.sh,)
error.bar(barx,hspDist,lower=low,upper=up,length=0.01)
legend("topright", col=c('blue','gray75','gray60','gray45','gray30','gray15'), legend=c("Thermo Thermopilius","E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic CCT"), pch=rep(15,6), cex=0.8)
graphics.off()


cairo_pdf('HSP_Protein_Dist_No_Ecoli.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
shades <- hcl(h=1:5 * (360 / 5), c=60, l=70 )
barx <- barplot(hspDist[2:6,], beside=TRUE, col=c(shades), ylim=c(0.0,0.3), names.arg=aalist.sh,)
error.bar(barx,hspDist[2:6,],lower=low[2:6,],upper=up[2:6,],length=0.01)
legend("topright", col=shades, legend=c("E.Coli GroEl","Thermo GroEl","Group II HSP","HSP90","Eukaryotic CCT"), pch=rep(15,5), cex=0.8)
graphics.off()


#plot3 The plot of E.Coli, Human, and Thermophilius
cutoff <- 0.3
psurf_ecoli <- fetchAllSurfResidues("ecoli_nogaps", cutoff, normalize=TRUE, "wenjunh")
psurf_thermus <- fetchAllSurfResidues("thermus_nogaps", cutoff, normalize=TRUE, "wenjunh")
psurf_human <- fetchAllSurfResidues("h2_w_nogaps", cutoff, normalize=TRUE, "whitead")

protein_dist <- matrix(0, 3, 20)
rownames(protein_dist) <- c("E.Coli","Thermus Thermopilius","Human")
colnames(protein_dist) <- colnames(psurf_ecoli)

protein_dist[1,] <- apply(psurf_ecoli, MARGIN=2, FUN=median)
protein_dist[2,] <- apply(psurf_thermus, MARGIN=2, FUN=median)
protein_dist[3,] <- apply(psurf_human, MARGIN=2, FUN=median)

up <- matrix(0,3,ncol(psurf_ecoli))
low <- matrix(0,3,ncol(psurf_ecoli))

for (i in 1:ncol(psurf_ecoli)) {
  up[1,i] <- quantile(psurf_ecoli[,i],0.95) - median(psurf_ecoli[,i])
  low[1,i] <- median(psurf_ecoli[,i]) - quantile(psurf_ecoli[,i],0.05)
}

for (i in 1:ncol(psurf_thermus)) {
  up[2,i] <- quantile(psurf_thermus[,i],0.95) - median(psurf_thermus[,i])
  low[2,i] <- median(psurf_thermus[,i]) - quantile(psurf_thermus[,i],0.05)
}

for (i in 1:ncol(psurf_human)) {
  up[3,i] <- quantile(psurf_human[,i],0.95) - median(psurf_human[,i])
  low[3,i] <- median(psurf_human[,i]) - quantile(psurf_human[,i],0.05)
}

cairo_pdf('Three_Protein_SurfRes.pdf',width=5.5, height=2.58, pointsize=9)
par(family='LMSans10', cex.axis=0.75)
barx <- barplot(protein_dist, beside=TRUE, col=c('gray75','gray50','gray25'), ylim=c(0.0,0.3), names.arg=aalist.sh,)
error.bar(barx,protein_dist,lower=low,upper=up,length=0.01)
legend("topright", col=c('gray75','gray50','gray25'), legend=c("E.Coli","Thermus Thermopilius","Human"), pch=rep(15,6), cex=0.8)
graphics.off()
