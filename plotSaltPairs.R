source("SQLShareLib.R")
kt <- 2.4943
bootstrap <- 100
conts <- fetchContacts("h1_surface_contacts.csv")

result <- matrix(rep(0,bootstrap * 9), nrow=bootstrap)

for(i in 1:bootstrap) {
  
  cmat <- sampleContacts(conts)

  pbackground <- cmat["TOTAL",] / sum(cmat["TOTAL",])

  #hydration
  hyd1 <- (cmat["WATER",] / t(apply(cmat[1:21,], MARGIN=2, sum))) * (1 - cmat["FREE",] / cmat["TOTAL",]) + cmat["FREE",] / cmat["TOTAL",]
  hyd2 <- cmat["WATER",] / t(apply(cmat[1:21,], MARGIN=2, sum))
  hyd2 <- hyd2 / (sum(cmat["WATER",]) / sum(cmat[1:21,]))
  hyd3 <- cmat["WATER",] / cmat["TOTAL",]

  #interactions
  pre <- cmat["GLU", "ARG"] / sum(cmat[1:21, "ARG"])
  per <- cmat["ARG", "GLU"] / sum(cmat[1:21, "GLU"])
  re.score <- pre / pbackground["GLU"]
  er.score <- per / pbackground["ARG"]
  
  prd <- cmat["ASP", "ARG"] / sum(cmat[1:21, "ARG"])
  pdr <- cmat["ARG", "ASP"] / sum(cmat[1:21, "ASP"])
  rd.score <- pre / pbackground["ASP"]
  dr.score <- per / pbackground["ARG"]
  
  pke <- cmat["GLU", "LYS"] / sum(cmat[1:21, "LYS"])
  pek <- cmat["LYS", "GLU"] / sum(cmat[1:21, "GLU"])
  ke.score <- pke / pbackground["GLU"]
  ek.score <- pek / pbackground["LYS"]
  
  pkd <- cmat["ASP", "LYS"] / sum(cmat[1:21, "LYS"])
  pdk <- cmat["LYS", "ASP"] / sum(cmat[1:21, "ASP"])
  kd.score <- pkd / pbackground["ASP"]
  dk.score <- pdk / pbackground["LYS"]

  pkk <- cmat["LYS", "LYS"] / sum(cmat[1:21, "LYS"])
  kk.score <- pdk / pbackground["LYS"]

  scores <- as.numeric(c(re.score + er.score, rd.score + dr.score, ke.score + ek.score, kd.score + dk.score, kk.score * 2))
  scores <- -kt * log(scores / 2)

  result[i,] <- c(as.numeric(-kt*log(hyd2[c("ARG", "LYS", "GLU", "ASP")])), scores)
}

save.image()

cairo_pdf("Interaction.pdf", width=4.5, height=3.3, pointsize=12)
par(family="LMRoman10", bg="transparent", cex=0.8, fg="dark gray")
yy <-apply(result[,c(5:9,1:4)], 2, mean)
ee <- sqrt(apply(result[,c(5:9,1:4)], 2, var))
barx <- barplot(yy, main="", xlab="Pair", ylab="Free Energy of Interaction",  ylim=c(min(yy) - max(ee),max(ee) + max(yy)), names.arg=c("RE", "RD", "KE", "KD", "KK", "RW", "KW", "DW", "EW"), col=c(rep("dark red", 4), "light gray", rep("blue", 4)), fg="black")
error.bar(barx,yy,ee,length=0.02)
legend("topright", col=c("dark red", "light gray", "blue"),  legend=c("Salt Bridge", "Known Unfavored", "Crystallographic Water"), text.col="black", cex=0.7, pch=15)
graphics.off()
