
mycent<-function (obj, xvar = NULL, 
                  cent = c(1,3,5,15,25,50,75,85,95,97,99), legend = TRUE,
                  ylab = "y", xlab = "x", main = NULL, 
                  main.gsub = "@", xleg = min(xvar), 
                  yleg = max(obj$y), xlim = range(xvar), 
                  ylim = range(obj$y), save = FALSE, plot = TRUE, points = TRUE, 
                  pch = "+", col = "blue", col.centiles = 1:length(cent) + 
                    2, lty.centiles = 1, lwd.centiles = 1, ...) 
{
  if (!is.gamlss(obj)) 
    stop(paste("This is not an gamlss object", "\n", ""))
  if (is.null(xvar)) 
    stop(paste("The xvar argument is not specified", "\n", 
               ""))
  fname <- obj$family[1]
  qfun <- paste("q", fname, sep = "")
  Title <- paste("Centile curves using", fname, sep = " ")
  main <- if (is.null(main)) 
    paste("Centile curves using", fname, sep = " ")
  else gsub(main.gsub, Title, main)
  oxvar <- xvar[order(xvar)]
  oyvar <- obj$y[order(xvar)]
  if (is.matrix(obj$y)) {
    oyvar <- obj$y[, 1][order(xvar)]
    ylim <- range(obj$y[, 1])
    yleg = max(obj$y[, 1])
  }
  if (plot) {
    lty.centiles <- rep(lty.centiles, length(cent))
    lwd.centiles <- rep(lwd.centiles, length(cent))
    col.centiles <- rep(col.centiles, length(cent))
    if (points == TRUE) {
      plot(oxvar, oyvar, type = "p", col = col, pch = pch, 
           xlab = xlab, ylab = ylab, xlim = xlim, ylim, 
           ...)
    }
    else {
      plot(oxvar, oyvar, type = "n", col = col, pch = pch, 
           xlab = xlab, ylab = ylab, xlim = xlim, ylim, 
           ...)
    }
    title(main)
  }
  col <- 3
  lpar <- length(obj$parameters)
  ii <- 0
  per <- rep(0, length(cent))
  for (var in cent) {
    if (lpar == 1) {
      newcall <- call(qfun, var/100,
                      mu = fitted(obj, "mu")[order(xvar)])
    }
    else if (lpar == 2) {
      newcall <- call(qfun, var/100, mu = fitted(obj, "mu")[order(xvar)], 
                      sigma = fitted(obj, "sigma")[order(xvar)])
    }
    else if (lpar == 3) {
      newcall <- call(qfun, var/100, mu = fitted(obj, "mu")[order(xvar)], 
                      sigma = fitted(obj,
                                     "sigma")[order(xvar)], nu = fitted(obj,"nu")[order(xvar)])
    }
    else {
      newcall <- call(qfun, var/100, mu = fitted(obj, "mu")[order(xvar)], 
                      sigma = fitted(obj, 
                                     "sigma")[order(xvar)], nu = fitted(obj, 
                                                                        "nu")[order(xvar)], tau = fitted(obj, "tau")[order(xvar)])
    }
    ii <- ii + 1
    ll <- eval(newcall)
    if (plot) {
      lines(oxvar, ll, col = col.centiles[ii], lty = lty.centiles[ii], 
            lwd = lwd.centiles[ii], ...)
    }
    per[ii] <- (1 - sum(oyvar > ll)/length(oyvar)) * 100
    if (!save) 
      cat("% of cases below ", var, "centile is ", per[ii], 
          "\n")
  }
  if (plot) {
    if (legend == TRUE) 
      legend(list(x = xleg, y = yleg), legend = cent, col = col.centiles, 
             lty = lty.centiles, lwd = lwd.centiles, ncol = 1, 
             ...)
  }
  if (save) {
    #return(cbind(cent, per))
     L=nu;Median=mu;S=sigma
    return(cbind(L,Median,S,cent, per))
  }
}

mycent.pred<-function(
  obj, type = c("centiles", "z-scores", "standard-centiles"), 
  xname = NULL, xvalues = NULL, power = NULL, yval = NULL, 
  cent = c(1,3,5,15,25,50,75,85,95,97,99), 
  dev = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
  plot = FALSE, legend = TRUE, 
  ...) 
{
  calc.cent <- function(xvar, cent) {
    o <- order(xvar)
    mat <- xvar[o]
    cent <- cent
    for (var in cent) {
      if (lpar == 1) {
        newcall <- call(qfun, var/100, mu = mu[o])
      }
      else if (lpar == 2) {
        newcall <- call(qfun, var/100, mu = mu[o], sigma = sigma[o])
      }
      else if (lpar == 3) {
        newcall <- call(qfun, var/100, mu = mu[o], sigma = sigma[o], 
                        nu = nu[o])
      }
      else {
        newcall <- call(qfun, var/100, mu = mu[o], sigma = sigma[o], 
                        nu = nu[o], tau = tau[o])
      }
      ll <- eval(newcall)
      mat <- cbind(mat, ll)
    }
    mat <- as.data.frame(mat)
    nnn <- paste("C", as.character(cent), sep = "")
    names(mat) <- c(xname, nnn)
    return(mat)
  }
  plot.mat <- function(mat, cent, legend, ...) {
    lcent <- dim(mat)[2]
    xleg <- min(mat[, 1])
    yleg <- max(mat[, 2:lcent])
    plot(mat[, 1], mat[, 2], type = "n", ...)
    for (i in 2:lcent) lines(mat[, 1], mat[, i], col = i)
    if (legend) 
      legend(list(x = xleg, y = yleg), 
             legend = cent, col = c(2,3, 4, 5, 6, 7, 8, 9, 10,
                                    11, 12), lty = 1, ncol = 1, bg = "white")
    invisible()
  }
  if (!is.gamlss(obj)) 
    stop(paste("This is not an gamlss object", "\n", ""))
  if (is.null(xvalues)) 
    stop(paste("The xvalues  argument is not specified", 
               "\n", ""))
  if (is.null(xname)) 
    stop(paste("The xname argument is not specified", "\n", 
               ""))
  if (!is.character(xname)) 
    stop(paste("The xname argument is not a character", "\n", 
               ""))
  xvar <- if (!is.null(power)) 
    xvar <- xvalues^power
  else xvalues
  newx <- data.frame(xvar)
  colnames(newx) <- xname
  lpar <- length(obj$parameters)
  if ("mu" %in% obj$parameters) {
    if (is.null(obj$mu.fix)) 
      mu <- predict(obj, what = "mu", newdata = newx, type = "response", 
                    ...)
    else if (obj$mu.fix == TRUE) 
      mu <- rep(fitted(obj, "mu")[1], length(xvar))
  }
  if ("sigma" %in% obj$parameters) {
    if (is.null(obj$sigma.fix)) 
      sigma <- predict(obj, what = "sigma", newdata = newx, 
                       type = "response", ...)
    else if (obj$sigma.fix == TRUE) 
      sigma <- rep(fitted(obj, "sigma")[1], length(xvar))
  }
  if ("nu" %in% obj$parameters) {
    if (is.null(obj$nu.fix)) 
      nu <- predict(obj, what = "nu", newdata = newx, type = "response", 
                    ...)
    else if (obj$nu.fix == TRUE) 
      nu <- rep(fitted(obj, "nu")[1], length(xvar))
  }
  if ("tau" %in% obj$parameters) {
    if (is.null(obj$tau.fix)) 
      tau <- predict(obj, what = "tau", newdata = newx, 
                     type = "response", ...)
    else if (obj$tau.fix == TRUE) 
      tau <- rep(fitted(obj, "tau")[1], length(xvar))
  }
  type <- match.arg(type)
  if (type == "centiles") {
    fname <- obj$family[1]
    qfun <- paste("q", fname, sep = "")
    xvar <- xvalues
    mat <- calc.cent(xvar = xvar, cent = cent)
    if (plot) 
      plot.mat(mat, cent, legend, ...)
    #return(mat)
   L=nu;Median=mu;S=sigma;Age<-mat[,1]
    return(cbind(Age,L,Median,S,mat[,-1]))
  }
  if (type == "z-scores") {
    if (is.null(yval)) 
      stop("the y values should be set if type=z-scores is used")
    if (length(yval) != length(xvalues)) 
      stop("length of xvalues and yval is not the same")
    fname <- obj$family[1]
    qfun <- paste("p", fname, sep = "")
    if (lpar == 1) {
      newcall <- call(qfun, yval, mu = mu)
    }
    else if (lpar == 2) {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma)
    }
    else if (lpar == 3) {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma, 
                      nu = nu)
    }
    else {
      newcall <- call(qfun, yval, mu = mu, sigma = sigma, 
                      nu = nu, tau = tau)
    }
    cdf <- eval(newcall)
    rqres <- qnorm(cdf)
    return(rqres)
  }
  if (type == "standard-centiles") {
    cent <- pnorm(dev) * 100
    fname <- obj$family[1]
    qfun <- paste("q", fname, sep = "")
    xvar <- xvalues
    mat <- calc.cent(xvar = xvar, cent = cent)
    nnn <- paste(as.character(dev), sep = "")
    names(mat) <- c(xname, nnn)
    if (plot) 
      plot.mat(mat, dev, legend, ...)
    #return(mat)
  L=nu;Median=mu;S=sigma;Age<-mat[,1]
    return(cbind(Age,L,Median,S,mat[,-1]))
  }
}