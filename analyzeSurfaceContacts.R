#! /usr/bin/env Rscript

source("SQLShareLib.R")

data.name <- "h1_surface_contacts.csv"

data <- fetchContacts(data.name, "whitead")

anames <- colnames(data[-c(1,2)])
rnum <- length(anames)

anames.ord <- order(anames)
data <- data[c(1,2,anames.ord + 2)]


bootstrap <- 500
interaction <- data.frame(matrix(rep(0,bootstrap * length(aalist.sh)), nrow=bootstrap))
colnames(interaction) <- aalist.sh

water <- data.frame(matrix(rep(0,bootstrap * length(aalist.sh)), nrow=bootstrap))
colnames(water) <- aalist.sh

for(i in 1:bootstrap) {
  counts <- sampleContacts(data)
  interaction[i,] <- 1 - counts["FREE",] / counts["TOTAL",]
  water[i,] <- counts["WATER",] / apply(counts[c(anames, "WATER"),], 2, sum)
  cat("\r",i,"/",bootstrap)
}
cat("\n")
save.image()

wyy <- apply(water, 2, function(x) {quantile(x,c(0.025,0.5,0.975), na.rm=T)})
iyy <- apply(interaction, 2, function(x) {quantile(x,c(0.025,0.5,0.975), na.rm=F)})

gi <- which(colnames(wyy) == "G")
wyy <- wyy[,-gi]
iyy <- iyy[,-gi]

print("prop interaction")
print(iyy[2,])

cairo_pdf("contact_plots.pdf", width=3.42,height=5.41, pointsize=8)
par(family="LMSans10", mfrow=c(2,1), cex.axis=0.60, mar=c(2,4,1,0.5))

xp <- barplot(iyy[2,], ylim=c(0,max(iyy[3,] + 0.03)), ylab="Proportion Interacting", col="dark gray")
error.bar(xp, iyy[2,], lower=(iyy[2,] - iyy[1,]), upper=(iyy[3,] - iyy[2,]), length=0.03)
text(0,max(iyy),"(a)")

par(mar=c(3,4,2,0.5))

xp <- barplot(wyy[2,], ylim=c(0,max(wyy[3,] + 0.03)), ylab="Water per Contact", col="dark gray", sub="b)")
error.bar(xp, wyy[2,], lower=(wyy[2,] - wyy[1,]), upper=(wyy[3,] - wyy[2,]), length=0.03)
text(0,max(wyy),"(b)")
graphics.off()

cairo_pdf("contact_plots_1.pdf", width=3.42,height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.6)

xp <- barplot(iyy[2,], ylim=c(0,max(iyy[3,] + 0.03)), ylab="Proportion Interacting", col="gray80", xlab="Amino Acid")
error.bar(xp, iyy[2,], lower=(iyy[2,] - iyy[1,]), upper=(iyy[3,] - iyy[2,]), length=0.03)
graphics.off()

cairo_pdf("contact_plots_2.pdf", width=3.42,height=2.58, pointsize=8)
par(family="LMSans10", cex.axis=0.6)
xp <- barplot(wyy[2,], ylim=c(0,max(wyy[3,] + 0.03)), ylab="Water per Contact", col="gray80", xlab="Amino Acid")
error.bar(xp, wyy[2,], lower=(wyy[2,] - wyy[1,]), upper=(wyy[3,] - wyy[2,]), length=0.03)
graphics.off()

