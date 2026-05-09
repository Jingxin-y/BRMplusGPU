#' Penalized Maximum‐Likelihood Estimation
#'
#' Penalizing the log-likelihood function with the Jeffry's prior
#'
#' Alternating coordinate‐descent updates are used to optimize over
#' \eqn{\alpha} and \eqn{\beta} in turn.
#'
#' @param param Character scalar, either \code{"RR"} or \code{"RD"}.
#' @param y Numeric vector of length \eqn{n}.  Binary outcomes (0/1).
#' @param x Numeric vector of length \eqn{n}.  Binary exposure indicator (0/1).
#' @param va Numeric matrix \eqn{n\times p_a}.
#' @param vb Numeric matrix \eqn{n\times p_b}.
#' @param alpha.start Numeric vector of length \eqn{p_a}.  Initial values for
#'   the \eqn{\alpha} parameters.
#' @param beta.start Numeric vector of length \eqn{p_b}.  Initial values for
#'   the \eqn{\beta} parameters.
#' @param weight Numeric vector of length \eqn{n}.  Observation weight.
#' @param max.step Integer.  Maximum number of alternating coordinate‐descent
#'   iterations.
#' @param thres Numeric.  Convergence threshold on relative parameter change.
#' @param pa Integer.  Number of \eqn{\alpha} parameters (\eqn{p_a}).
#' @param pb Integer.  Number of \eqn{\beta} parameters (\eqn{p_b}).
#'
max.likelihood.jeffrey = function(param, y, x, va, vb, alpha.start, beta.start, weight,
                          max.step, thres, pa, pb) {

  startpars = c(alpha.start, beta.start)

  getProb = if (param == "RR") getProbRR else getProbRD

  ## negative log likelihood function
  neg.log.likelihood = function(pars) {
    alpha = pars[1:pa]
    beta = pars[(pa + 1):(pa + pb)]
    p0p1 = getProb(va %*% alpha, vb %*% beta)
    p0 = p0p1[, 1];   p1 = p0p1[, 2]

   if (param == "RR") fisher  = var.mle.rr (x, alpha.start, beta.start, va, vb, weight) else  fisher  = var.mle.rd (x, alpha.start, beta.start, va, vb, weight)

    return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                  (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) - sum((1 - y[x ==
                                                                                  1]) * log(1 - p1[x == 1]) * weight[x == 1] + (y[x == 1]) * log(p1[x ==
                                                                                                                                                       1]) * weight[x == 1])+
             log(det(fisher))/2) ### add the Jeffrey's prior

  }

  neg.log.likelihood.alpha = function(alpha){
    p0p1 = getProb(va %*% alpha, vb %*% beta)
    p0    = p0p1[,1];  p1 = p0p1[,2]

    if (param == "RR") fisher  = var.mle.rr (x, alpha.start, beta.start, va, vb, weight) else  fisher  = var.mle.rd (x, alpha.start, beta.start, va, vb, weight)

    return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                  (y[x==0])*log(p0[x==0])*weight[x==0]) -
             sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                   (y[x==1])*log(p1[x==1])*weight[x==1])+
             log(det(fisher))/2)  ### add the Jeffrey's prior
  }

  neg.log.likelihood.beta = function(beta){
    p0p1 = getProb(va %*% alpha, vb %*% beta)
    p0    = p0p1[,1];  p1 = p0p1[,2]

    if (param == "RR") fisher  = var.mle.rr (x, alpha.start, beta.start, va, vb, weight) else  fisher  = var.mle.rd (x, alpha.start, beta.start, va, vb, weight)

    return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                  (y[x==0])*log(p0[x==0])*weight[x==0]) -
             sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                   (y[x==1])*log(p1[x==1])*weight[x==1])+
             log(det(fisher))/2)  ### add the Jeefrey's prior
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
