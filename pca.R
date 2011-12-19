#!/usr/bin/r


pointFile = argv[1]
data <- read.table(pointFile)

matrix <- as.matrix(data)
smatrix <- matrix(rep(0,9), nrow=3)

for(i in 1:nrow(data)) {
  smatrix <- smatrix + matrix[i,] %*% t(matrix[i,])
}
smatrix <- smatrix / (nrow(data) + 1)

eigenVecs <- eigen(smatrix)

e.values <- eigenVecs$values
cat(e.values)
cat('\n')
