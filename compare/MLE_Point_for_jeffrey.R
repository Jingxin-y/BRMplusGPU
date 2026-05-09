#' Penalized Maximum‐Likelihood Estimation
#'
#' Penalizing the log-likelihood function with the Jeffry's prior, but with the prior directly applied to p0,p1
#'
#' Compute Determinant of Fisher Information for RR/RD Model
#'
#' @param param Character scalar, either \code{"RR"} or \code{"RD"}.
#' @param x Binary exposure indicator (0/1).
#' @param alpha.ml Numeric vector of length \(p_a\).  Fitted \(\alpha\) parameters.
#' @param beta.ml Numeric vector of length \(p_b\).  Fitted \(\beta\) parameters.
#' @param va Numeric matrix \(n\times p_a\).  Design matrix for the \(\alpha\) component.
#' @param vb Numeric matrix \(n\times p_b\).  Design matrix for the \(\beta\) component.
#' @param weight Numeric vector of length \(n\).  Observation weight (not used in this simple approximation).
#'
### augmentation calculation
fisher.detf = function(param, x, alpha.ml, beta.ml, va, vb, weight) {

  getProb = if (param == "RR") getProbRR else getProbRD

  p0p1 = getProb(va %*% alpha.ml, vb %*% beta.ml)
  p0 = p0p1[x == 0, 1]
  p1 = p0p1[x == 1, 2]

  fisher.det = sum(1/(p0*(1-p0)))*sum(1/(p1*(1-p1)))
  return(fisher.det)
}


#' Penalized Maximum‐Likelihood Estimation
#'
#' Penalizing the log-likelihood function with the Jeffry's prior, but with the prior directly applied to p0,p1
#'
#' @param param Character scalar, either \code{"RR"} or \code{"RD"}.
#' @param y Numeric vector of length \(n\).  Binary outcomes (0/1).
#' @param x Numeric vector of length \(n\).  Binary exposure indicator (0/1).
#' @param va Numeric matrix \(n\times p_a\).
#' @param vb Numeric matrix \(n\times p_b\).
#' @param alpha.start Numeric vector of length \(p_a\).  Initial values for
#'   the \eqn{\alpha} parameters.
#' @param beta.start Numeric vector of length \eqn{p_b}.  Initial values for
#'   the \eqn{\beta} parameters.
#' @param weight Numeric vector of length \(n\).  Observation weights.
#' @param max.step Integer.  Maximum number of alternating updates.
#' @param thres Numeric.  Convergence threshold on relative parameter change.
#' @param pa Integer.  Number of \(\alpha\) parameters (\(p_a\)).
#' @param pb Integer.  Number of \(\beta\) parameters (\(p_b\)).
#'

# The difference between this file and "MLE_Point_of_estimator_for_jeffrey.R" lies in the function used to compute the Fisher information.
# We can merge the two files by adding a conditional statement based on the value of argument "method".
max.likelihood.jeffrey.direct = function(param, y, x, va, vb, alpha.start, beta.start, weight,
                          max.step, thres, pa, pb) {

  startpars = c(alpha.start, beta.start)

  getProb = if (param == "RR") getProbRR else getProbRD

  ## negative log likelihood function
  neg.log.likelihood = function(pars) {
    alpha = pars[1:pa]
    beta = pars[(pa + 1):(pa + pb)]
    p0p1 = getProb(va %*% alpha, vb %*% beta)
    p0 = p0p1[, 1];   p1 = p0p1[, 2]

    fisher.det  = fisher.detf(param, x, alpha.start, beta.start, va, vb, weight)

    return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                  (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) - sum((1 - y[x ==
                                                                                  1]) * log(1 - p1[x == 1]) * weight[x == 1] + (y[x == 1]) * log(p1[x ==
                                                                                                                                                       1]) * weight[x == 1])-
             log(fisher.det)/2)

  }

  neg.log.likelihood.alpha = function(alpha){
    p0p1 = getProb(va %*% alpha, vb %*% beta)
    p0    = p0p1[,1];  p1 = p0p1[,2]

    fisher.det  = fisher.detf(param, x, alpha.start, beta.start, va, vb, weight)

    return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                  (y[x==0])*log(p0[x==0])*weight[x==0]) -
             sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                   (y[x==1])*log(p1[x==1])*weight[x==1])-
             log(fisher.det)/2)
  }

  neg.log.likelihood.beta = function(beta){
    p0p1 = getProb(va %*% alpha, vb %*% beta)
    p0    = p0p1[,1];  p1 = p0p1[,2]

    fisher.det  = fisher.detf(param, x, alpha.start, beta.start, va, vb, weight)

    return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                  (y[x==0])*log(p0[x==0])*weight[x==0]) -
             sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                   (y[x==1])*log(p1[x==1])*weight[x==1])-
             log(fisher.det)/2)
  }


  ## Optimization

  Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
  alpha = alpha.start; beta = beta.start
  diff = thres + 1; step = 0
  while(diff > thres & step < max.step){
    step = step + 1
    opt1 = stats::optim(alpha,neg.log.likelihood.alpha,control=list(maxit=max(100,max.step/10)))
    diff1 = Diff(opt1$par,alpha)
    alpha = opt1$par
    opt2 = stats::optim(beta,neg.log.likelihood.beta,control=list(maxit=max(100,max.step/10)))
    diff  = max(diff1,Diff(opt2$par,beta))
    beta = opt2$par
  }

  opt = list(par = c(alpha,beta), convergence = (step < max.step),
             value = neg.log.likelihood(c(alpha,beta)), step = step)

  return(opt)
}


