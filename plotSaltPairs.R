source("SQLShareLib.R")
kt <- 2.4943
bootstrap <- 3
conts <- fetchContacts("h1_surface_contacts.csv")

result <- matrix(rep(0,bootstrap * 9), nrow=bootstrap)

for(i in 1:bootstrap) {
  
  cmat <- sampleContacts(conts)

  pbackground <- apply(cmat[1:20, 1:20], 1, sum) / sum(cmat[1:20,1:20])

  #hydration
  hyd1 <- (cmat["WATER",] / t(apply(cmat[1:21,], MARGIN=2, sum))) * (1 - cmat["FREE",] / cmat["TOTAL",]) + cmat["FREE",] / cmat["TOTAL",]
  hyd2 <- cmat["WATER",] / t(apply(cmat[1:21,], MARGIN=2, sum))
  hyd3 <- cmat["WATER",] / cmat["TOTAL",]

  #interactions
  pre <- cmat["GLU", "ARG"] / cmat["TOTAL", "ARG"]
  per <- cmat["ARG", "GLU"] / cmat["TOTAL", "GLU"]
  re.score <- pre / pbackground["GLU"]
  er.score <- per / pbackground["ARG"]
  
  prd <- cmat["ASP", "ARG"] / cmat["TOTAL", "ARG"]
  pdr <- cmat["ARG", "ASP"] / cmat["TOTAL", "ASP"]
  rd.score <- pre / pbackground["ASP"]
  dr.score <- per / pbackground["ARG"]
  
  pke <- cmat["GLU", "LYS"] / cmat["TOTAL", "LYS"]
  pek <- cmat["LYS", "GLU"] / cmat["TOTAL", "GLU"]
  ke.score <- pke / pbackground["GLU"]
  ek.score <- pek / pbackground["LYS"]
  
  pkd <- cmat["ASP", "LYS"] / cmat["TOTAL", "LYS"]
  pdk <- cmat["LYS", "ASP"] / cmat["TOTAL", "ASP"]
  kd.score <- pkd / pbackground["ASP"]
  dk.score <- pdk / pbackground["LYS"]

  pkk <- cmat["LYS", "LYS"] / cmat["TOTAL", "LYS"]
  kk.score <- pdk / pbackground["LYS"]

  scores <- as.numeric(c(re.score + er.score, rd.score + dr.score, ke.score + ek.score, kd.score + dk.score, kk.score * 2))
  scores <- kt * log(scores / 2)

  result[i,] <- c(as.numeric(hyd3[c("ARG", "LYS", "GLU", "ASP")]), scores)
}

save.image()

cairo_pdf("Interaction.pdf", width=4.5, height=3.3, pointsize=12)
par(family="LMRoman10", bg="transparent", cex=0.8)
yy <-apply(result[,5:9], 2, mean)
ee <- sqrt(apply(result[,5:9], 2, var))
barx <- barplot(yy, main="", xlab="Pair", ylab="Free Energy of Interaction", ylim=c(0,max(ee) + max(yy)), names.arg=c("RE", "RD", "KE", "KD", "KK"))
error.bar(barx,yy,ee,length=0.02)
graphics.off()
