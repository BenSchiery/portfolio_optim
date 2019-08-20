rm(list = ls())
while(length(dev.list()) > 0){dev.off()}
library(quadprog) # need the solve.QP() function
library(Matrix) # need the nearPD() function

# a function to find the range of expected returns that are possible given individual assets' expected 
# returns (rx) and weight constraints
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
  
  # similarly, apply as much weight as is allowable to the assets with the highest individual return rates
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
eff.frontier <- function(returns, min.allo = 0, max.allo = 1, 
                         weight.sum = 1, n.step = 100, contraction = 0.99, 
                         include.lower.branch = FALSE){
  n.asset <- dim(returns)[2]
  if(is.null(colnames(returns))){
    colnames(returns) <- paste("A", 1:n.asset, sep = "") # name the assets if they aren't already
  }
  if(min.allo * n.asset > weight.sum){stop("Minimum allocation too high for given weight sum constraint.")}
  if(max.allo * n.asset < weight.sum){stop("Maximum allocation too low for given weight sum constraint.")}
  
  Dmat <- matrix(nearPD(cov(returns))$mat, ncol = n.asset) # quadratic term
  cov.mat.norm <- norm(Dmat)
  
  # if not including the lower branch, we find the weights which minimize 
  # portfolio SD without regard for portfolio returns, then find the return 
  # rate of a portfolio having those weights; this will be the lower bound 
  # on the sequence of return rates for which we will minimize portfolio SD
  if(!include.lower.branch){
    dvec <- matrix(0, nrow = n.asset, ncol = 1) # a dvec of zeroes makes us disregard returns
    Amat <- cbind("sum.col" = rep(1, times = n.asset),
                  "min.allo" = diag(n.asset),
                  "max.allo" = -diag(n.asset)) # constraint matrix without a return constraint
    bvec <- c(weight.sum,
              rep(min.allo, times = n.asset), 
              -rep(max.allo, times = n.asset)) # bvec includes no goal return rate 
    meq <- 1 # only the first constraint demands equality, sum(weights) = weight.sum
    w <- solve.QP(Dmat = Dmat / cov.mat.norm,
                  dvec = dvec, 
                  Amat = Amat, 
                  bvec = bvec, 
                  meq = meq)$solution # the weights which minimize portfolio SD given our constraints
    rx.lower <- sum(w * colMeans(returns)) # soon-to-be lower bound on expected returns range
  }
  
  dvec <- colMeans(returns, na.rm = T) # linear term
  Amat <- cbind("sum.col" = rep(1, times = n.asset),
                "rx.col" = dvec,
                "min.allo" = diag(n.asset),
                "max.allo" = -diag(n.asset)) # constraint matrix _with_ a return constraint
  # 2 equality constraints; sum(weights) = weight.sum AND 
  # sum(weights * asset expected returns) = desired portfolio return (which are sapply'ed over)
  meq <- 2
  # range of possible expected return values for portfolios of the given assets
  rx.rng <- rx.range(rx = dvec,
                     min.allo = min.allo,
                     max.allo = max.allo,
                     weight.sum = weight.sum) 
  
  if(!include.lower.branch){
    rx.rng[1] <- rx.lower
  } # replace the lower bound if needed
  
  rx.rng <- (rx.rng - mean(rx.rng)) * contraction + mean(rx.rng) # shrink the range a little bit
  # generate a sequence of exp. ret. values at which to minimize portfolio stdev
  rx.seq <- seq(rx.rng[1], rx.rng[2], length.out = n.step) 
  
  ftr <- t(sapply(X = rx.seq,
                  FUN = function(rx){
                    # change bvec to reflect changing goal for return rate
                    bvec <- c(weight.sum, rx, 
                              rep(min.allo, times = n.asset), 
                              -rep(max.allo, times = n.asset))
                    w <- solve.QP(Dmat = Dmat / cov.mat.norm,
                                  dvec = dvec / cov.mat.norm, 
                                  Amat = Amat, 
                                  bvec = bvec, 
                                  meq = meq)$solution # optimal weights
                    stdev <- sqrt(t(w) %*% Dmat %*% w) # portfolio SD at optimal weights
                    c(w, "sd" = stdev, "rx" = rx)
                  })) # sapply solve.QP over all the rx values
  colnames(ftr) <- c(colnames(returns), "sd", "rx") # give names
  
  ftr
}

folio.optim <- function(rx.goal, returns, min.allo = 0, max.allo = 1, weight.sum = 1){
  n.asset <- dim(returns)[2]
  if(is.null(colnames(returns))){
    colnames(returns) <- paste("A", 1:n.asset, sep = "")
  }
  if(min.allo * n.asset > weight.sum){stop("Minimum allocation too high for given weight sum constraint.")}
  if(max.allo * n.asset < weight.sum){stop("Maximum allocation too low for given weight sum constraint.")}
  
  Dmat <- matrix(nearPD(cov(returns))$mat, ncol = n.asset)
  cov.mat.norm <- norm(Dmat)
  dvec <- colMeans(returns)
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
  soln <- solve.QP(Dmat = Dmat / cov.mat.norm,
                   dvec = dvec / cov.mat.norm, 
                   Amat = Amat, 
                   bvec = bvec, 
                   meq = meq)$solution
  names(soln) <- colnames(returns)
  soln
}

folio <- function(returns){
  E <- new.env()
  E$rx <- colMeans(returns)
  E$cov.mat <- matrix(nearPD(cov(returns))$mat, ncol = dim(returns)[2])
  function(weights){
    c("sd" = as.numeric(sqrt(t(weights) %*% E$cov.mat %*% weights)),
      "rx" = sum(E$rx * weights))
  }
}

###########################################
#### Generate a time series of returns ####
###########################################

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

ftr <- eff.frontier(returns = time.series, 
                    min.allo = 0, 
                    max.allo = 1, 
                    weight.sum = 1, 
                    include.lower.branch = F)

psd <- ftr[,"sd"] # portfolio sd
prx <- ftr[,"rx"] # portfolio exp. ret.
bbox <- c("xmin" = min(0, psd),
          "xmax" = max(0, psd),
          "ymin" = min(0, prx),
          "ymax" = max(0, prx)) # bounding box for graph
plot(x = psd, y = prx, type = "l", lwd = 2, 
     xlim = bbox[1:2], ylim = bbox[3:4],
     xlab = "Portfolio SD", ylab = "Portfolio Expected Return",
     main = "Efficient Frontier")
points(x = 0, y = 0, pch = 16, col = "red") # add origin for reference

########################
#### Test optimizer ####
########################

rx.goal <- mean(prx[prx > 0])
w <- folio.optim(rx.goal = rx.goal, returns = time.series)
pt <- folio(time.series)(weights = w)
points(x = pt[1], y = pt[2], pch = 16, col = "cyan", cex = 2)
