
## ======================================================================
## Functions for the estimation of CAR processes using R.
##
## Spatial Statistics M - Spring 2013
## Created: Peter Craigmile, peter.craigmile@glasgow.ac.uk, Feb 2013
## Contact: pfc@stat.osu.edu
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


CAR.precision <- function (alpha, tau2, W, wts = rep(1, nrow(W))) {
  ## ======================================================================
  ## Calculate the precision matrix for a CAR model
  ## with parameters 'alpha', 'tau2',
  ## and proximity matrix 'W', and weights 'wts'.
  ## ======================================================================

  diag(wts) %*% (diag(nrow(W)) - alpha * W) / tau2
}



CAR.alpha.range <- function (W, wts = rep(1, nrow(W))) {
  ## ======================================================================
  ## Calculate the range of a possible alpha values for a CAR model
  ## with proximity matrix 'W' and weights 'wts'.
  ## ======================================================================
  
  1/range(eigen(W)$values)
}



CAR.sim <- function (n, mu, alpha, tau2, W, wts = rep(1, nrow(W))) {
  ## ======================================================================
  ## Simulate 'n' realizations of a CAR model with mean mu,
  ## autocorrelation parameter 'alpha' and variance 'tau2'.
  ## In this model 'W' is the proximity matrix, and weights 'wts'.
  ## ======================================================================

  require(MASS)
  mvrnorm(n, mu, solve(CAR.precision(alpha, tau2, W, wts)))
}



CAR.mle.beta <- function (alpha, y, W, wts = rep(1, nrow(W)),
                          X=cbind(rep(1, length(y)))) {
  ## ======================================================================
  ## Calculate the MLE of beta
  ## ======================================================================

  A <- crossprod(X, CAR.precision(alpha, 1, W, wts))
  
  solve(A %*% X) %*% A %*% y
}



CAR.mle.tau2 <- function (alpha, y, W, wts = rep(1, nrow(W)),
                          X=cbind(rep(1, length(y)))) {
  ## ======================================================================
  ## Calculate the MLE of tau
  ## ======================================================================
  
  prec   <- CAR.precision(alpha, 1, W, wts)
  A      <- crossprod(X, prec)
  mu.hat <- X %*% solve(A %*% X) %*% A %*% y
  z      <- y - as.numeric(mu.hat)
  
  drop( (crossprod(z, prec) %*% z)/length(y) )
}



CAR.m2l.alpha <- function (alpha, y, W, wts = rep(1, nrow(W)),
                           X=cbind(rep(1, length(y))), penalty=0,
                           debug=FALSE) {
  ## ======================================================================
  ## Calculate minus two times the profile log likelihood plus a penalty
  ## for alpha
  ## ======================================================================

  if (debug)
    cat("alpha: ", alpha, "\n")

  n      <- length(y)
  prec   <- CAR.precision(alpha, 1, W, wts)
  A      <- crossprod(X, prec)
  z      <- y - X %*% solve(A %*% X) %*% A %*% y

  drop(n * log(2*pi) + n * log((crossprod(z, prec) %*% z)/n) -
       log(det(prec)) + penalty)
}



CAR.mle.alpha <- function (y, W, wts = rep(1, nrow(W)),
                           X=cbind(rep(1, length(y))),
                           debug=FALSE) {
  ## ======================================================================
  ## Calculate the MLE of alpha.
  ## ======================================================================
  
  optimize(CAR.m2l.alpha, CAR.alpha.range(W, wts),
           y=y, W=W, X=X, wts=wts, debug=debug)$min
}


CAR.cov.mle.beta <- function (alpha, tau2, W, wts = rep(1, nrow(W)),
                              X=cbind(rep(1, nrow(W)))) {
  ## ======================================================================
  ## Calculate the estimated covariance matrix for the MLE of beta.
  ## ======================================================================

  solve(t(X) %*% CAR.precision(alpha, tau2, W, wts) %*% X)
}




CAR.mle <- function (y, W, wts = rep(1, nrow(W)),
                     X=cbind(rep(1, nrow(W)))) {
  ## ======================================================================
  ## Estimate the parameters of the CAR model based on the data 'y',
  ## proximity matrix 'W', weights 'wts', and design matrix 'X'.
  ## ======================================================================

  
  alpha.hat <- CAR.mle.alpha(y, W, wts, X=X)
  tau2.hat  <- CAR.mle.tau2(alpha.hat, y, W, wts, X=X)
  beta.hat  <- CAR.mle.beta(alpha.hat, y, W, wts, X=X)

  cov.beta.hat <- CAR.cov.mle.beta(alpha.hat, tau2.hat, W, wts, X=X)

  npars <- length(beta.hat) + 2

  the.loglik <- -0.5*CAR.m2l.alpha(alpha.hat, y, W, wts, X=X)

  structure(list(y=y, W=W, wts=wts, X=X,
                 alpha.hat = alpha.hat,
                 tau2.hat  = tau2.hat,
                 beta.hat  = beta.hat,
                 cov.beta.hat = cov.beta.hat,
                 loglik = the.loglik,
                 AIC = -2*the.loglik + 2*npars),
            class="CAR")
}



summary.CAR <- function (x) {
  ## ======================================================================
  ## Summarize this CAR model
  ## ======================================================================

  est  <- x$beta.hat
  se   <- sqrt(diag(x$cov.beta.hat))
  tval <- est/se
  error.df <- length(x$y) - length(est)

  ans <- list()
  ans$coef <- cbind(est, se, tval, 2 * pt(abs(tval),
                                          error.df, lower.tail = FALSE))
  dimnames(ans$coef)[[2]] <-
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  ans$alpha.hat <- x$alpha.hat
  ans$tau2.hat  <- x$tau2.hat
  ans$error.df <- error.df
  ans$loglik    <- x$loglik
  ans$AIC       <- x$AIC
  ans
}


