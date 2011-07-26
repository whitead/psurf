#' Represents Correlation circles
#'
#' @author Taiyun Wei
#' @param corr Correlation matrix to represent
#' @param col vector the fill color of circles from 1 to -1
#'        the length of it may not be 2, eg rainbow(50)
#' @param bg background color of graph
#' @param cex numeric, for the variable names
#' @param order whether reorder the variables using principal 
#'         component analysis of the correlation matrix
#' @param title title of the graph
#' @param ... extra parameters, currenlty ignored
circle.corr <- function(corr, bars=NULL, col=c("black","white"), bg = "white", 
	cex = 1, order = FALSE, title = "", ...){
	
     rootbg <- bg    

    if (is.null(corr)) 
        return(invisible())
    if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 
        6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1)) 
        stop("Need a correlation matrix!")
    n <- nrow(corr)
    m <- ncol(corr)
    
    ## reorder the variables using principal component analysis
    if (order) {
    	if(!n==m){
    		stop("The matrix must be squre if order is TRUE!")
    	}
      x.eigen <- eigen(corr)$vectors[, 1:2]
      e1 <- x.eigen[, 1]
      e2 <- x.eigen[, 2]
      alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
      corr <- corr[order(alpha), order(alpha)]
      if(!is.null(bars)) {
        bars <- bars[order(alpha)]
      }
    }
    
    ## set up variable names
    rname <- rownames(corr)
    cname <- colnames(corr)
    if (is.null(rname)) 
        rname <- 1:n
    if (is.null(cname)) 
        cname <- 1:m
    rname <- as.character(rname)
    cname <- as.character(cname)

    ## calculate label-text width approximately
    par(mar = c(0, 0, 2, 0), bg = "white")
    plot.new()
    plot.window(c(0, m), c(0, n), asp = 1)
    xlabwidth <- max(strwidth(rname, cex = cex))
    ylabwidth <- max(strwidth(cname, cex = cex))
    
    barlength <- m / 4
    
    if(is.null(bars)) {
      barlength <- 0
    } 

    ## set up an empty plot with the appropriate dimensions
    plot.window(c(-xlabwidth + 0.5, m + barlength + 0.5), c(-barlength, n + 1 + ylabwidth),
                asp = 1, xlab="", ylab="")
    rect(0.5, 0.5, m + 0.5, n + 0.5, col = bg)	##background color

    ## add variable names and title
    text(rep(-xlabwidth/2, n), n:1, rname, col = "black", cex = cex)
    text(1:m, rep(n + 1 + ylabwidth/2, m), cname, srt = 90, col = "black", 
        cex = cex )
    title(title)

    ## add grid
    segments(rep(0.5, n + 1), 0.5 + 0:n, rep(m + 0.5, n + 1), 
        0.5 + 0:n, col = "gray")
    segments(0.5 + 0:m, rep(0.5, m + 1), 0.5 + 0:m, rep(n + 0.5, 
        m), col = "gray")

    ## assign circles' fill color
    nc <- length(col)
    if(nc==1)
        	bg <- rep(col, n*m)
    else{
        ff <- seq(-1,1, length=nc+1)
        bg2 = rep(0, n * m)
        for (i in 1:(n * m)){
          bg2[i] <- rank(c(ff[2:nc], as.vector(corr)[i]), 
                            ties.method = "random")[nc]
        }
        bg <- (col[nc:1])[bg2]
    }

    ## plot n*m circles using vector language, suggested by Yihui Xie
    ## the area of circles denotes the absolute value of coefficient
    symbols(rep(1:m, each = n), rep(n:1, m), add = TRUE, inches = F, 
        circles = as.vector(sqrt(abs(corr))/2), bg = as.vector(bg))
	
    ##Cover up diagonal
    if(m == n) {
        symbols(1:m, n:1, add=T, inches=F, circles=rep(0.25, m), bg=rootbg, fg=rootbg)
    }

    ##Add bars
    if(!is.null(bars) & n == m) {
      barmarg  <- 0
      col <- rgb(0.4, 0.4, 0.8)
      for(i in 1:m) {
        rect(m + 0.5, 0.5 + (m - i) + barmarg, 0.5 + m + barlength * bars[i] / max(bars), 0.5 + (m - i) + 1 - barmarg, col=col, border="gray")
	rect(-0.5 + i + barmarg, 0.5, -0.5 + i + 1 - barmarg, 0.5 - barlength * bars[i] / max(bars), col=col, border="gray")
      }
    }
}
## examples
#data(mtcars)
#circle.corr( cor(mtcars), order = TRUE, bg = "gray50", 
#	col = colorRampPalette(c("blue","white","red"))(100) )


