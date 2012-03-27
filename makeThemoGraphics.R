source("SQLShareLib.R")
source("ProtLib.R")

#Interaction Table
dataset <- "thermus"
username <- "wenjunh"

chis <- getInteractionEnergy(dataset, username)

interactionTable <- chis

surface <- fetchAllSurfResidues(paste(dataset,"nogaps", sep="_"), cutoff, normalize=FALSE, username)
interior <- fetchAllBuriedResidues(paste(dataset,"nogaps", sep="_"), cutoff, normalize=FALSE, username)
all <- fetchAllSurfResidues(paste(dataset,"nogaps", sep="_"), -1, normalize=FALSE, username)

distributions <- list(Proteins=all, Surface=surface, Interior=interior, Diff=(surface - all))



#sample all the contact energies
contacts <- fetchContacts(paste(paste(dataset,"backbone","contacts", sep="_"), "csv" ,sep="."), username)
bootstrap <- 3
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
