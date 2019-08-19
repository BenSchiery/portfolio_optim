rm(list = ls())
while(length(dev.list()) > 0){dev.off()}
library(quadprog) # need the solve.QP() function
library(Matrix) # need the nearPD() function

# a function to find the range of expected returns that are possible given individual assets' expected 
# returns and weight constraints
rx.range <- function(rx, min.allo = 0, max.allo = 1, weight.sum = 1){
  # to get the minimum expected return of the portfolio, apply as much weight as is allowable to the 
  # assets with the lowest individual return rates
  rx <- sort(rx) # lowest to highest individual asset expected return values
  x <- weight.sum - length(rx) * min.allo # how much of the total weight we have left to apply to assets
  rx.min <- sum(rx * min.allo) # start by applying minimum weight to eash asset
  for(i in rx){
    if(x == 0){break} # once we have no more weight to apply, move on
    y <- min(max.allo, x) # find the max allowable weight, unless we have less than that much left
    rx.min <- rx.min + i * y # apply that weight
    x <- x - y # decrement the amount of weight we have left to apply
  }
  
  rx <- sort(rx, decreasing = T)
  x <- weight.sum - length(rx) * min.allo
  rx.max <- sum(rx * min.allo)
  for(i in rx){
    if(x == 0){break}
    y <- min(max.allo, x)
    rx.max <- rx.max + i * y
    x <- x - y
  }
  
  c(rx.min, rx.max)
}

# a function to give a series of weight vectors constituting portfolios on the efficient frontier
# subject to weight constraints.
# minimizes portfolio standard deviation at a sequence of portfolio expected return values
# `contraction` shrinks the range of exp. ret. values checked by a small amount because
# trouble tends to arrise when the solver hits the extremes of that range
eff.frontier <- function(time.series, min.allo = 0, max.allo = 1, 
                         weight.sum = 1, n.step = 100, contraction = 0.99, 
                         include.lower.branch = FALSE){
  n.asset <- dim(time.series)[2]
  if(is.null(colnames(time.series))){
    colnames(time.series) <- paste("A", 1:n.asset, sep = "") # name the assets if they aren't already
  }
  if(min.allo * n.asset > weight.sum){stop("Minimum allocation too high for given weight sum constraint.")}
  if(max.allo * n.asset < weight.sum){stop("Maximum allocation too low for given weight sum constraint.")}
  
  Dmat <- matrix(nearPD(cov(time.series))$mat, ncol = n.asset) # quadratic term
  cov.mat.norm <- norm(Dmat)
  dvec <- colMeans(time.series) # linear term
  Amat <- cbind("sum.col" = rep(1, times = n.asset),
                "rx.col" = dvec,
                "min.allo" = diag(n.asset),
                "max.allo" = -diag(n.asset)) # constraint matrix
  # 2 equality constraints; sum(weights) = weight.sum AND 
  # sum(weights * asset expected returns) = desired portfolio return (which are sapply'ed over)
  meq <- 2
  # range of possible expected return values for portfolios of the given assets
  rx.rng <- rx.range(rx = dvec,
                     min.allo = min.allo,
                     max.allo = max.allo,
                     weight.sum = weight.sum)  
  
  rx.rng <- (rx.rng - mean(rx.rng)) * contraction + mean(rx.rng) # shrink the range a little bit
  # generate a sequence of exp. ret. values at which to minimize portfolio stdev
  rx.seq <- seq(rx.rng[1], rx.rng[2], length.out = n.step) 
  
  ftr <- t(sapply(X = rx.seq,
                  FUN = function(rx){
                    bvec <- c(weight.sum, rx, 
                              rep(min.allo, times = n.asset), 
                              -rep(max.allo, times = n.asset))
                    w <- solve.QP(Dmat = Dmat / cov.mat.norm,
                                  dvec = dvec / cov.mat.norm, 
                                  Amat = Amat, 
                                  bvec = bvec, 
                                  meq = meq)$solution
                    stdev <- sqrt(t(w) %*% Dmat %*% w)
                    c(w, "sd" = stdev, "rx" = rx)
                  })) # sapply solve.QP over all the rx values
  colnames(ftr) <- c(colnames(time.series), "sd", "rx") # give names
  if(include.lower.branch){
    ftr
  }else{
    lower.rx <- ftr[,"rx"][which.min(ftr[,"sd"])]
    ftr[ftr[,"rx"] >= lower.rx,]
  }
}

folio.optim <- function(rx.goal, time.series, min.allo = 0, max.allo = 1, weight.sum = 1){
  n.asset <- dim(time.series)[2]
  if(is.null(colnames(time.series))){
    colnames(time.series) <- paste("A", 1:n.asset, sep = "")
  }
  if(min.allo * n.asset > weight.sum){stop("Minimum allocation too high for given weight sum constraint.")}
  if(max.allo * n.asset < weight.sum){stop("Maximum allocation too low for given weight sum constraint.")}
  
  Dmat <- matrix(nearPD(cov(time.series))$mat, ncol = n.asset)
  dvec <- colMeans(time.series)
  Amat <- cbind("sum.col" = rep(1, times = n.asset),
                "rx.col" = dvec,
                "min.allo" = diag(n.asset),
                "max.allo" = -diag(n.asset))
  meq <- 2
  rx.rng <- rx.range(rx = dvec,
                     min.allo = min.allo,
                     max.allo = max.allo,
                     weight.sum = weight.sum)
  
  if(rx.goal < rx.rng[1]){
    stop(paste("Expected return goal too low. Try a goal between", 
               round(rx.rng[1], digits = 2), "and", round(rx.rng[2], digits = 2),"."))
  }else if(rx.goal > rx.rng[2]){
    stop(paste("Expected return goal too high. Try a goal between", 
               round(rx.rng[1], digits = 2), "and", round(rx.rng[2], digits = 2),"."))
  }
  
  bvec <- c(weight.sum, rx.goal, 
            rep(min.allo, times = n.asset), 
            -rep(max.allo, times = n.asset))
  soln <- solve.QP(Dmat = Dmat,
                   dvec = dvec, 
                   Amat = Amat, 
                   bvec = bvec, 
                   meq = meq)$solution
  names(soln) <- colnames(time.series)
  soln
}

folio <- function(time.series){
  E <- new.env()
  E$rx <- colMeans(time.series)
  E$cov.mat <- matrix(nearPD(cov(time.series))$mat, ncol = dim(time.series)[2])
  function(weights){
    c("sd" = as.numeric(sqrt(t(weights) %*% E$cov.mat %*% weights)),
      "rx" = sum(E$rx * weights))
  }
}

################################
#### Generate a time series ####
################################

n.asset <- 10

time.series <- sapply(X = 1:n.asset,
                      FUN = function(i){
                        cumsum(rnorm(n = n.asset^2, sd = 10))
                      })
asset.names <- paste("A", 1:n.asset, sep = "")
colnames(time.series) <- asset.names

#####################################
#### Plot the efficient frontier ####
#####################################

ftr <- eff.frontier(time.series = time.series, min.allo = 0, max.allo = 1, weight.sum = 1)

psd <- ftr[,"sd"] # portfolio sd
prx <- ftr[,"rx"] # portfolio exp. ret.
x.min <- min(0, psd)
x.max <- max(0, psd)
y.min <- min(0, prx)
y.max <- max(0, prx)
plot(x = psd, y = prx, type = "l", lwd = 2, 
     xlim = c(x.min, x.max), ylim = c(y.min, y.max),
     xlab = "Portfolio SD", ylab = "Portfolio Expected Return",
     main = "Efficient Frontier")
points(x = 0, y = 0, pch = 16, col = "red") # add origin for reference

########################
#### Test optimizer ####
########################

rx.goal <- mean(prx[prx > 0])
w <- folio.optim(rx.goal = rx.goal, time.series = time.series)
pt <- folio(time.series)(weights = w)
points(x = pt[1], y = pt[2], pch = 16, col = "cyan", cex = 2)
