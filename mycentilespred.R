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
    #TAU=tau;L=nu;Median=mu;S=sigma;Age<-mat[,1]
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
    # TAU=tau;L=nu;Median=mu;S=sigma;Age<-mat[,1]
    L=nu;Median=mu;S=sigma;Age<-mat[,1]
        return(cbind(Age,L,Median,S,mat[,-1]))
  }
}
