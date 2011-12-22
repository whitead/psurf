#!/usr/bin/r


#read in the points
pointFile = argv[1]
data <- read.table(pointFile)
dmatrix <- as.matrix(data)
dmatrix <- dmatrix - matrix(apply(dmatrix, 2, mean), byrow=T, nrow=nrow(dmatrix), ncol=3)
smatrix <- matrix(rep(0,9), nrow=3)
#transform it into the radius of gyration tensor
for(i in 1:nrow(data)) {
  smatrix <- smatrix + dmatrix[i,] %*% t(dmatrix[i,])
}
smatrix <- smatrix / (nrow(data) + 1)

#get the eigenvalues
eigenVecs <- eigen(smatrix)

e.values <- eigenVecs$values
cat(e.values)
cat('\n')
