
## ======================================================================
## Functions for the areal modeling in R.
##
## Spatial Statistics M - Spring 2013
## Created: Peter Craigmile, peter.craigmile@glasgow.ac.uk, Feb 2013
## Contact: pfc@stat.osu.edu
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



Euclidean.dist.matrix <- function (sites) {
  ## ======================================================================
  ## Purpose: The function calculate the Euclidean distances among sites.
  ## Assumes: 'sites' are matrices or data frames.
  ## ======================================================================

  dd <- as.matrix(dist(sites, upper=TRUE, diag=TRUE))
  dimnames(dd) <- NULL
  dd
}



Moran.I <- function (y, W) {
  ## ======================================================================
  ## Purpose: Calculate Moran's I for data 'y' and proximity matrix 'W'.
  ## ======================================================================

  z <- y - mean(y)
  
  length(z) * sum(W * outer(z,z)) / (sum(W) * sum(z^2))
}





Moran.I.perm.test <- function (y, W, show.perms=TRUE, num.perms=1000) {
  ## ======================================================================
  ## Purpose: Perform a two sided hypothesis test for spatial
  ##          dependence in areal data using Moran's I for
  ##          data 'y' and proximity matrix 'W'.
  ## ======================================================================
  
  obs.I <- Moran.I(y, W)

  ## Permute the data 1000 times, and calculate the I statistics.
  perm.stats <- sapply(1:num.perms, function (k) Moran.I(sample(y), W))
  
  ## Produce a histogram of these statistics
  if (show.perms) {
    hist(perm.stats, breaks=seq(-1, 1, 0.05), main="",
         xlab="Permutation statistics for Moran's I")
    abline(v=obs.I, lwd=2)
  }
  
  ## Calculate the two sided Monte Carlo p-value
  p.value <- 2*sum(perm.stats > abs(obs.I)) / (length(perm.stats)+1)
  
  list(statistic=obs.I,
       parameter="I",
       p.value=p.value)
}




plot.poly <- function (xx.polylist, aux.var, intervals,
                       legend.x, legend.y, ...) {
  ## ======================================================================
  ## Plotting spatial polygons, with an auxiliary variable 'aux.var',
  ## broken down by the 'intervals'.
  ## ======================================================================

  if (missing(aux.var)) {
    plot(xx.polylist, ...)
         ##forcefill=FALSE, ...)
  }
  else {
    cols <- grey(seq(0.2, 0.8, length=length(intervals)))

    the.cols <- cols[findInterval(aux.var, intervals)]
    plot(xx.polylist, col=the.cols, ...)

    ys <- sort(seq(legend.y[1], legend.y[2], len=length(intervals)))
    
    image(range(legend.x), ys,
          rbind(intervals, intervals), col=cols, add=T)
    
    text(min(legend.x), ys, intervals, pos=2, cex=0.9)
  }

  invisible()
}

