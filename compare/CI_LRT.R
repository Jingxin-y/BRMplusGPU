#' Construct Likelihood-Ratio Confidence Intervals via Profiling (for \eqn{\alpha})
#'
#' Computes **profile-likelihood** confidence intervals and LRT p-values for each
#' component of \eqn{\alpha} in a binary-outcome model with treatment indicator
#' \eqn{x \in \{0,1\}}. For each \eqn{\alpha_j}, the function fixes \eqn{\alpha_j}
#' at a grid of values, maximizes the log-likelihood over the remaining
#' \eqn{\alpha_{-j}} and all \eqn{\beta}, evaluates the likelihood-ratio statistic,
#' and finds where it falls below the \eqn{\chi^2_1} cutoff to form a
#' \eqn{95\%} profile CI.
#'
#' @param param Character. Model family switch: use \code{"RR"} for a relative-risk
#'   parametrization (via \code{getProbRR}), otherwise \code{"RD"} for a
#'   risk-difference parametrization (via \code{getProbRD}).
#' @param y Numeric vector of length \eqn{n}. Binary outcomes (0/1).
#' @param x Numeric vector of length \eqn{n}. Binary exposure indicator (0/1).
#' @param va Numeric matrix \eqn{n \times p_a}. Design for \eqn{\alpha}.
#' @param vb Numeric matrix \eqn{n \times p_b}. Design for \eqn{\beta}.
#' @param weight Numeric vector of length \eqn{n}. Observation weights.
#' @param max.step Integer. Maximum number of alternating (block) optimization
#'   iterations used when profiling each \eqn{\alpha_j}.
#' @param thres Numeric. Convergence tolerance; the inner alternation stops when
#'   the relative parameter change is below this value.
#' @param pars Numeric vector of length \eqn{p_a + p_b}. Concatenated MLEs
#'   \code{c(alpha.ml, beta.ml)} used to center search ranges and for the null
#'   (unconstrained) fit in the LRT.
#' @param se Numeric vector (typically of length \eqn{p_a}) giving marginal
#'   standard errors for \eqn{\alpha} at the MLE; used to build the profiling
#'   grid \eqn{\alpha_j \in [\alpha_{j,ml} \pm 3\,se_j]} (truncated to \code{[-12, 12]}).
#' @param pa Integer. Number of \eqn{\alpha} parameters (\eqn{p_a}).
#' @param pb Integer. Number of \eqn{\beta} parameters (\eqn{p_b}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{low}}{Numeric vector of length \eqn{p_a}. Lower 95\% profile CI
#'         for each \eqn{\alpha_j}.}
#'   \item{\code{up}}{Numeric vector of length \eqn{p_a}. Upper 95\% profile CI
#'         for each \eqn{\alpha_j}.}
#'   \item{\code{p}}{Numeric vector of length \eqn{p_a}. LRT p-values for testing
#'         \eqn{H_0:\ \alpha_j=0} vs \eqn{H_1:\ \alpha_j \neq 0}.}
#' }
#'

profile <- function(param,y, x, va, vb, weight, max.step, thres, pars, se, pa, pb){
  ## real data
  getProb = if (param == "RR") getProbRR else getProbRD
  alpha.ml = pars[1:pa]
  beta.ml = pars[(pa + 1):(pa + pb)]
  p0p1 = getProb(va %*% alpha.ml, vb %*% beta.ml)
  p0.ml = p0p1[, 1];   p1.ml = p0p1[, 2]
  ## profile

  alpha.start <- rep(0,pa)
  beta.start <- rep(0,pb)

  optm.beta <- function(alphaj,j){

    neg.log.likelihood = function(pars) {
      alpha = pars[1:pa]
      beta = pars[(pa + 1):(pa + pb)]
      p0p1 = getProb(va %*% alpha, vb %*% beta)
      p0 = p0p1[, 1];   p1 = p0p1[, 2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)

      return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                    (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) - sum((1 - y[x ==
                                                                                    1]) * log(1 - p1[x == 1]) * weight[x == 1] + (y[x == 1]) * log(p1[x ==
                                                                                                                                                         1]) * weight[x == 1]))
    }

    neg.log.likelihood.alpha = function(alpha){
      p0p1 = getProb(va %*% alpha, vb %*% beta)
      p0    = p0p1[,1];  p1 = p0p1[,2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)

      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))
    }

    neg.log.likelihood.beta = function(beta){
      p0p1 = getProb(va %*% alpha, vb %*% beta)
      p0    = p0p1[,1];  p1 = p0p1[,2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)


      return(-sum((1-y[x==0])*log(1-p0[x==0])*weight[x==0] +
                    (y[x==0])*log(p0[x==0])*weight[x==0]) -
               sum((1-y[x==1])*log(1-p1[x==1])*weight[x==1] +
                     (y[x==1])*log(p1[x==1])*weight[x==1]))
    }

    Diff = function(x,y) sum((x-y)^2)/sum(x^2+thres)
    alpha = alpha.start
    alpha[j] = alphaj
    beta = beta.start
    diff = thres + 1; step = 0
    while(diff > thres & step < max.step){
      step = step + 1
      opt1 = stats::optim(alpha,neg.log.likelihood.alpha,control=list(maxit=max(100,max.step/10)))
      diff1 = Diff(opt1$par,alpha)
      alpha = opt1$par
      alpha[j] = alphaj
      opt2 = stats::optim(beta,neg.log.likelihood.beta,control=list(maxit=max(100,max.step/10)))
      diff  = max(diff1,Diff(opt2$par,beta))
      beta = opt2$par
    }
    return(neg.log.likelihood(c(alpha,beta)))
  }

  LRT.alpha <- function(alpha,j){
    return(2*optm.beta(alpha.ml[j],j)-2*optm.beta(alpha[j],j))
  }

  get.lrt <- function(alpha){
    lrt <- rep(0,length(alpha))
    for (j in 1:length(alpha)) {
      lrt[j] <- LRT.alpha(alpha,j)
    }
    return(lrt)
  }

  chi.th <- qchisq(0.95, df = 1)

  alpha.seq <- lapply(1:pa, function(j){
    seq(max(alpha.ml[j] - 3*se[j], -12),
        min(alpha.ml[j] + 3*se[j],  12),
        length.out = 40)
  })
  alpha.mat <- do.call(cbind, alpha.seq)

  result.lrt <- apply(alpha.mat, 1, get.lrt)

  lrt.mat <- if(pa>1) {t(result.lrt)} else{as.matrix(result.lrt,ncol = 1)}

  alpha.up = rep(0,pa)
  for (j in 1:pa) {
    alpha.up[j] <- max(alpha.mat[which(lrt.mat[,j] <= chi.th),j])
  }
  alpha.low = rep(0,pa)
  for (j in 1:pa) {
    alpha.low[j] <- min(alpha.mat[which(lrt.mat[,j] <= chi.th),j])
  }

  p.values <- pchisq(get.lrt(alpha.start), df = 1, lower.tail = FALSE)
  return(list(low = alpha.low,
              up = alpha.up,
              p = p.values))
}


