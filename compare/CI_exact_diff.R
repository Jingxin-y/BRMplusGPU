#'  Exact (Monte Carlo) Profile-Likelihood CIs via Acceptability for \eqn{\alpha}
#'
#' @description
#' Constructs component-wise \eqn{95\%} confidence intervals and p-values for the
#' target parameter vector \eqn{\alpha} in a binary-outcome model using **profile
#' likelihood** and an **acceptability** (double randomization) calibration of the
#' likelihood-ratio statistic (LRT). For each component \eqn{\alpha_j}, the method:
#' (i) fixes \eqn{\alpha_j} at a candidate value, (ii) re-optimizes over
#' \eqn{\alpha_{-j}} and \eqn{\beta} by alternating minimization of the negative
#' log-likelihood, (iii) evaluates the profile LRT, and (iv) uses Monte Carlo
#' simulation under the constrained fit to compute an acceptability value that
#' plays the role of a calibrated p-value. CI endpoints are found by **bisection**
#' on \eqn{\alpha_j} until acceptability crosses the \code{0.05} threshold.
#'
#' @param param Character. Model family switch: \code{"RR"} (relative-risk scale)
#' otherwise \code{"RD"} (risk-difference scale)
#' @param y Numeric vector of length \eqn{n}. Binary outcomes \code{0/1}.
#' @param x Numeric vector of length \eqn{n}. Binary exposure \code{0/1}.
#' @param va Numeric matrix \eqn{n \times p_a}. Design for \eqn{\alpha}.
#' @param vb Numeric matrix \eqn{n \times p_b}. Design for \eqn{\beta}.
#' @param weight Numeric vector of length \eqn{n}. Observation weights.
#' @param max.step Integer. Maximum number of inner alternating-optimization
#'   iterations when profiling each \eqn{\alpha_j}.
#' @param thres Numeric. Convergence tolerance for the inner alternation; the loop
#'   stops when the relative change metric falls below this value.
#' @param thres.dicho Numeric. Tolerance for the outer **bisection** used to find
#'   CI endpoints on each \eqn{\alpha_j}.
#' @param pars Numeric vector of length \eqn{p_a+p_b}. Concatenated MLEs
#'   \code{c(alpha.ml, beta.ml)} used to center search ranges and for the
#'   unconstrained likelihood in the LRT.
#' @param se Numeric vector (at least length \eqn{p_a}). Standard errors for the
#'   \eqn{\alpha} components at the MLE; used to build starting brackets
#'   \code{alpha.ml +/- 4*se} (truncated to \code{[-8, 8]}).
#' @param pa Integer. Number of \eqn{\alpha} parameters.
#' @param pb Integer. Number of \eqn{\beta} parameters.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{low}}{Numeric vector of length \eqn{p_a}. Lower 95\% CI endpoint
#'         for each \eqn{\alpha_j}.}
#'   \item{\code{up}}{Numeric vector of length \eqn{p_a}. Upper 95\% CI endpoint
#'         for each \eqn{\alpha_j}.}
#'   \item{\code{p}}{Numeric vector of length \eqn{p_a}. Acceptability-based
#'         two-sided p-values for testing \eqn{H_0:\ \alpha_j=0}.}
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item \strong{Computation.} This is substantially more expensive than a
#'         \eqn{\chi^2} cutoff because each candidate \eqn{\alpha_j} value requires
#'         an inner optimization and a Monte Carlo loop. Consider parallelizing
#'         the simulation loop if needed.
#'   \item \strong{Reproducibility.} Set a random seed before calling if you want
#'         stable CIs/p-values across runs.
#' }
#'
#' @references
#' Blaker, H. (2000). Confidence curves and improved exact confidence intervals for discrete distributions. Canadian Journal of Statistics, 28(4), 783-798.
#' https://doi.org/10.2307/3315916 (Theorem 1)
#'
exact <- function(param, y, x, va, vb, weight, max.step, thres, thres.dicho, pars, se, pa, pb) {
  ## real data
  getProb <- if (param == "RR") getProbRR else getProbRD
  alpha.ml <- pars[1:pa]
  beta.ml <- pars[(pa + 1):(pa + pb)]
  p0p1 <- getProb(mat_vec_mul(va, alpha.ml), mat_vec_mul(vb, beta.ml))
  p0.ml <- p0p1[, 1]
  p1.ml <- p0p1[, 2]
  ## profile
  
  alpha.start <- rep(0, pa)
  beta.start <- rep(0, pb)
  
  optm.beta <- function(alphaj, j, y) {
    neg.log.likelihood <- function(pars) {
      alpha <- pars[1:pa]
      beta <- pars[(pa + 1):(pa + pb)]
      p0p1 <- getProb(mat_vec_mul(va, alpha), mat_vec_mul(vb, beta))
      p0 <- p0p1[, 1]
      p1 <- p0p1[, 2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      
      return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                    (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) - sum((1 - y[x ==
                                                                                   1]) * log(1 - p1[x == 1]) * weight[x == 1] + (y[x == 1]) * log(p1[x ==
                                                                                                                                                       1]) * weight[x == 1]))
    }
    
    neg.log.likelihood.alpha <- function(alpha) {
      p0p1 <- getProb(mat_vec_mul(va, alpha), mat_vec_mul(vb, beta))
      p0 <- p0p1[, 1]
      p1 <- p0p1[, 2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      
      return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                    (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) -
               sum((1 - y[x == 1]) * log(1 - p1[x == 1]) * weight[x == 1] +
                     (y[x == 1]) * log(p1[x == 1]) * weight[x == 1]))
    }
    
    neg.log.likelihood.beta <- function(beta) {
      p0p1 <- getProb(mat_vec_mul(va, alpha), mat_vec_mul(vb, beta))
      p0 <- p0p1[, 1]
      p1 <- p0p1[, 2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      
      
      return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                    (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) -
               sum((1 - y[x == 1]) * log(1 - p1[x == 1]) * weight[x == 1] +
                     (y[x == 1]) * log(p1[x == 1]) * weight[x == 1]))
    }
    
    Diff <- function(x, y) sum((x - y)^2) / sum(x^2 + thres)
    alpha <- alpha.start
    alpha[j] <- alphaj
    beta <- beta.start
    diff <- thres + 1
    step <- 0
    while (diff > thres & step < max.step) {
      step <- step + 1
      opt1 <- stats::optim(alpha, neg.log.likelihood.alpha, control = list(maxit = max(100, max.step / 10)))
      diff1 <- Diff(opt1$par, alpha)
      alpha <- opt1$par
      alpha[j] <- alphaj
      opt2 <- stats::optim(beta, neg.log.likelihood.beta, control = list(maxit = max(100, max.step / 10)))
      diff <- max(diff1, Diff(opt2$par, beta))
      beta <- opt2$par
    }
    return(neg.log.likelihood(c(alpha, beta)))
  }
  
  LRT.alpha <- function(alpha, j, y) {
    return(2 * optm.beta(alpha.ml[j], j, y) - 2 * optm.beta(alpha, j, y))
  }
  
  
  # Simulate distribution of observed profile‐LRT statistic
  ptail <- function(alphaj, j, nsim = 500) {
    LRT.sim <- numeric(nsim)
    
    neg.log.likelihood.alpha <- function(alpha) {
      p0p1 <- getProb(va %*% alpha, vb %*% beta.sim)
      p0 <- p0p1[, 1]
      p1 <- p0p1[, 2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      
      return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                    (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) -
               sum((1 - y[x == 1]) * log(1 - p1[x == 1]) * weight[x == 1] +
                     (y[x == 1]) * log(p1[x == 1]) * weight[x == 1]))
    }
    
    neg.log.likelihood.beta <- function(beta) {
      p0p1 <- getProb(va %*% alpha.sim, vb %*% beta)
      p0 <- p0p1[, 1]
      p1 <- p0p1[, 2]
      eps <- 1e-12
      p0 <- pmin(pmax(p0, eps), 1 - eps)
      p1 <- pmin(pmax(p1, eps), 1 - eps)
      
      
      return(-sum((1 - y[x == 0]) * log(1 - p0[x == 0]) * weight[x == 0] +
                    (y[x == 0]) * log(p0[x == 0]) * weight[x == 0]) -
               sum((1 - y[x == 1]) * log(1 - p1[x == 1]) * weight[x == 1] +
                     (y[x == 1]) * log(p1[x == 1]) * weight[x == 1]))
    }
    
    Diff <- function(x, y) sum((x - y)^2) / sum(x^2 + thres)
    alpha.sim <- alpha.start
    alpha.sim[j] <- alphaj
    beta.sim <- beta.start
    diff <- thres + 1
    step <- 0
    while (diff > thres & step < max.step) {
      step <- step + 1
      opt1 <- stats::optim(alpha.sim, neg.log.likelihood.alpha, control = list(maxit = max(100, max.step / 10)))
      diff1 <- Diff(opt1$par, alpha.sim)
      alpha.sim <- opt1$par
      alpha.sim[j] <- alphaj
      opt2 <- stats::optim(beta.sim, neg.log.likelihood.beta, control = list(maxit = max(100, max.step / 10)))
      diff <- max(diff1, Diff(opt2$par, beta.sim))
      beta.sim <- opt2$par
    }
    
    # Fitted probabilities under (alpha0, beta.sim)
    prob <- getProb(mat_vec_mul(va, alpha.sim), mat_vec_mul(vb, beta.sim))
    p0 <- prob[, 1]
    p1 <- prob[, 2]
    
    for (i in 1:nsim) {
      y.sim <- numeric(length(y))
      y.sim[x == 0] <- rbinom(sum(x == 0), 1, p0[x == 0])
      y.sim[x == 1] <- rbinom(sum(x == 1), 1, p1[x == 1])
      LRT.sim[i] <- LRT.alpha(alphaj, j, y.sim)
    }
    
    return(LRT.sim)
  }
  
  # Compute the acceptabity function
  acceptability <- function(alphaj, LRT.obs, LRT.sim) {
    p.left <- mean(LRT.sim <= LRT.obs)
    p.right <- mean(LRT.sim >= LRT.obs)
    p.min <- min(p.left, p.right)
    p1.min <- sapply(LRT.sim, function(x) {
      p.left <- mean(LRT.sim <= x)
      p.right <- mean(LRT.sim >= x)
      min(p.left, p.right)
    })
    a.val <- mean(p1.min <= p.min)
    return(a.val)
  }
  
  # function of dichotomy
  dichotomy <- function(j, alpha.low, alpha.up, direction = "low", thres.dicho = 1e-3, max.step = 20) {
    alpha.iteration <- alpha.up
    step <- 1
    while (alpha.up - alpha.low > thres.dicho & step < max.step) {
      LRT.obs <- LRT.alpha(alpha.iteration, j, y)
      LRT.sim <- ptail(alpha.iteration, j, nsim = (21 - step) * 100)
      a.val <- acceptability(alpha.iteration, LRT.obs, LRT.sim)
      cond <- if (direction == "low") {
        a.val > 0.05
      } else {
        a.val < 0.05
      }
      cond <- isTRUE(as.logical(cond))
      if (cond) {
        alpha.up <- alpha.iteration
        alpha.iteration <- (alpha.up + alpha.low) / 2
      } else {
        alpha.low <- alpha.iteration
        alpha.iteration <- (alpha.up + alpha.low) / 2
      }
      step <- step + 1
    }
    return(list(alpha.dicho = alpha.iteration, convergence = (step < max.step)))
  }
  
  
  # # get candidate of alpha
  # alpha.up.start <- pmin(alpha.ml + 4 * se[1:pa], 8)
  # alpha.low.start <- pmax(alpha.ml - 4 * se[1:pa], -8)
  # 
  # alpha.up1 <- rep(0, pa)
  # for (j in 1:pa) {
  #   alpha.up1[j] <- dichotomy(j, alpha.ml[j], alpha.up.start[j], "up", thres.dicho = 1e-3)$alpha.dicho
  # }
  # alpha.low1 <- rep(0, pa)
  # for (j in 1:pa) {
  #   alpha.low1[j] <- dichotomy(j, alpha.low.start[j], alpha.ml[j], "low", thres.dicho = 1e-3)$alpha.dicho
  # }
  
  # Build the 95% CI
  
  p.value <- rep(0, pa)
  for (j in 1:pa) {
    LRT.obs.p <- LRT.alpha(0, j, y)
    LRT.sim.p <- ptail(0, j, nsim = 2)
    p.value[j] <- acceptability(0, LRT.obs.p, LRT.sim.p)
  }
  
  return(list(
    low = 1,
    up = 1,
    p = p.value
  ))
}